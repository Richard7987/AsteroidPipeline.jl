using AsteroidPipeline
using Aqua
using TypedTables
using HTTP
using WCS
using FITSIO
using Random
using Test
using Statistics
using Reproject
using Distributed
using Dates: DateTime, datetime2julian

# Loads .env (gitignored) for local runs, e.g. ASTROMETRY_API_KEY for the live
# plate_solve test; variables already set in the environment win.
let env_file = joinpath(pkgdir(AsteroidPipeline), ".env")
    if isfile(env_file)
        for line in eachline(env_file)
            line = strip(line)
            (isempty(line) || startswith(line, "#")) && continue
            k, v = split(line, "="; limit=2)
            k = strip(k)
            v = strip(strip(v), ['"', '\''])
            haskey(ENV, k) || (ENV[k] = v)
        end
    end
end

@testset "AsteroidPipeline.jl" begin

    @testset "package quality (Aqua)" begin
        Aqua.test_all(AsteroidPipeline)
    end

    @testset "detect_sources" begin
        Random.seed!(1)
        image = 100.0 .+ 5.0 .* randn(64, 64)

        x0, y0, amp, sigma = 32, 40, 500.0, 2.0
        for j in axes(image, 1), i in axes(image, 2)
            r2 = (i - x0)^2 + (j - y0)^2
            image[j, i] += amp * exp(-r2 / (2 * sigma^2))
        end

        sources = detect_sources(image; threshold=5.0)
        @test length(sources) >= 1
        best = sources[argmax(sources.peak)]
        @test best.x == x0
        @test best.y == y0
        @test best.flux > 0
        @test best.flux_err > 0

        # a finite gain adds the source's shot noise: it may only widen flux_err,
        # and only where there is a source
        with_gain = detect_sources(image; threshold=5.0, gain=2.0)
        best_gain = with_gain[argmax(with_gain.peak)]
        @test best_gain.flux_err > best.flux_err

        noise_only = 100.0 .+ 5.0 .* randn(32, 32)
        @test length(detect_sources(noise_only; threshold=100.0)) == 0
    end

    @testset "detect_sources xcen/ycen and mask; fill_value_mask" begin
        Random.seed!(11)
        image = 100.0 .+ 5.0 .* randn(80, 60)   # (y, x)
        x0, y0 = 20.3, 40.7                     # deliberately off-grid
        for j in axes(image, 2), i in axes(image, 1)
            image[i, j] += 800.0 * exp(-((j - x0)^2 + (i - y0)^2) / (2 * 1.5^2))
        end
        # a detector-gap-style constant strip, as in real IASC/PS1 frames
        image[:, 40:45] .= 300.0

        mask = fill_value_mask(image)
        @test all(mask[:, 40:45])
        @test !any(mask[:, 1:30])            # real noisy sky is never a plateau

        masked = detect_sources(image; threshold=5.0, mask)
        @test !any(d -> 38 <= d.x <= 47, masked)
        star = masked[argmax(masked.peak)]
        @test star.xcen ≈ x0 atol=0.2                     # sub-pixel, unlike x/y
        @test star.ycen ≈ y0 atol=0.2
        @test star.x == round(Int, x0) && star.y == round(Int, y0)
        @test_throws DimensionMismatch detect_sources(image; threshold=5.0, mask=falses(3, 3))

        # a cosmic-ray hit: one pixel, nothing in its neighbours
        image[15, 50] += 800.0
        all_peaks = detect_sources(image; threshold=5.0, mask)
        cr = only(filter(d -> d.x == 50 && d.y == 15, all_peaks))
        @test abs(cr.sharpness) < 0.1
        @test star.sharpness > 0.5
        cleaned = detect_sources(image; threshold=5.0, mask, min_sharpness=0.3)
        @test !any(d -> d.x == 50 && d.y == 15, cleaned)
        @test any(d -> d.x == star.x && d.y == star.y, cleaned)
    end

    @testset "stack_difference" begin
        Random.seed!(12)
        nx, ny = 120, 100
        truth = WCSTransform(2; crpix=[60.0, 50.0], crval=[150.0, 20.0],
                             cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
        stars_xy = [(30.0, 25.0), (90.0, 70.0), (55.0, 80.0), (100.0, 20.0)]
        dithers = [(0, 0), (3, -2), (-4, 1), (2, 5)]       # whole-pixel pointing offsets
        transparency = [1.0, 1.0, 1.0, 0.7]               # one hazier frame
        mover = [(40.0 + 8k, 60.0 - 5k) for k in 0:3]      # frame-1 grid
        images = Matrix{Float64}[]; wcss = WCSTransform[]
        for k in 1:4
            dx, dy = dithers[k]
            w = WCSTransform(2; crpix=truth.crpix .+ [dx, dy], crval=truth.crval,
                             cdelt=truth.cdelt, ctype=truth.ctype)
            img = 100.0 .+ 4.0 .* randn(ny, nx)
            add!(x, y, amp) = for j in 1:nx, i in 1:ny
                img[i, j] += amp * exp(-((j - x)^2 + (i - y)^2) / (2 * 1.5^2))
            end
            for (sx, sy) in stars_xy
                add!(sx + dx, sy + dy, 3000.0 * transparency[k])
            end
            add!(mover[k][1] + dx, mover[k][2] + dy, 300.0 * transparency[k])
            push!(images, img); push!(wcss, w)
        end

        diffs, masks = stack_difference(images, wcss; bright_star_sigma=Inf)
        for k in 1:4
            dx, dy = dithers[k]
            noise = std(diffs[k][.!masks[k]])
            mx, my = round(Int, mover[k][1] + dx), round(Int, mover[k][2] + dy)
            @test diffs[k][my, mx] > 20 * noise            # the mover survives...
            for (sx, sy) in stars_xy                       # ...static stars don't
                @test abs(diffs[k][round(Int, sy + dy), round(Int, sx + dx)]) < 6 * noise
            end
        end
        # default bright-star masking covers each star's core and halo
        _, masked = stack_difference(images, wcss)
        @test all(masked[k][round(Int, sy + dithers[k][2]), round(Int, sx + dithers[k][1])]
                  for k in 1:4, (sx, sy) in stars_xy)
        @test_throws ArgumentError stack_difference(images[1:2], wcss[1:2])

        # sub-pixel alignment: same contract, mover still found
        sp, spm = stack_difference(images, wcss; bright_star_sigma=Inf, subpixel=true)
        for k in 1:4
            mx, my = round(Int, mover[k][1] + dithers[k][1]), round(Int, mover[k][2] + dithers[k][2])
            @test sp[k][my, mx] > 10 * std(sp[k][.!spm[k]])
        end
        B = AsteroidPipeline._bilinear
        A = [1.0 2.0; 3.0 4.0]
        @test B(A, nothing, 1.5, 1.5) ≈ 2.5 && B(A, nothing, 1.0, 2.0) ≈ 2.0
        @test isnan(B(A, nothing, 2.5, 1.0)) && isnan(B(A, BitMatrix([0 0; 0 1]), 1.5, 1.5))
    end

    @testset "refine_wcs" begin
        Random.seed!(13)
        nx, ny = 300, 260
        truth = WCSTransform(2; crpix=[150.0, 130.0], crval=[3.05, -10.74],
                             cdelt=[-0.257 / 3600, 0.257 / 3600], ctype=["RA---TAN", "DEC--TAN"])
        image = 100.0 .+ 5.0 .* randn(ny, nx)
        xy = [(20 + 260rand(), 20 + 220rand()) for _ in 1:40]
        for (x, y) in xy, j in max(1, floor(Int, x) - 8):min(nx, ceil(Int, x) + 8),
            i in max(1, floor(Int, y) - 8):min(ny, ceil(Int, y) + 8)
            image[i, j] += 2000.0 * exp(-((j - x)^2 + (i - y)^2) / (2 * 1.6^2))
        end
        sky = [pix_to_world(truth, [x, y]) for (x, y) in xy]
        stars = Table(ra=[s[1] for s in sky], dec=[s[2] for s in sky])
        # a header solution off by ~26 px in y, like the real PS1 headers
        header = WCSTransform(2; crpix=truth.crpix .+ [-1.0, 26.5], crval=truth.crval,
                              cdelt=truth.cdelt, ctype=truth.ctype)

        r = refine_wcs(image, header, stars)
        @test r.refined
        @test r.n_matches >= 30
        @test r.offset_px[1] ≈ 1.0 atol=0.5          # detected minus header-predicted
        @test r.offset_px[2] ≈ -26.5 atol=0.5
        @test r.rms_arcsec < 0.05
        for (x, y) in xy[1:5]
            a, b = pix_to_world(r.wcs, [x, y]), pix_to_world(truth, [x, y])
            @test hypot((a[1] - b[1]) * cosd(b[2]), a[2] - b[2]) * 3600 < 0.05
        end

        too_few = refine_wcs(image, header, stars[1:3])
        @test !too_few.refined
        @test too_few.wcs === header

        # a damaged frame: thousands of sources, none of them the catalogue's stars
        damaged = 100.0 .+ 5.0 .* randn(ny, nx)
        for _ in 1:2000
            x, y = 1 + (nx - 1)rand(), 1 + (ny - 1)rand()
            for j in max(1, floor(Int, x) - 4):min(nx, ceil(Int, x) + 4), i in max(1, floor(Int, y) - 4):min(ny, ceil(Int, y) + 4)
                damaged[i, j] += 2000.0 * exp(-((j - x)^2 + (i - y)^2) / (2 * 1.0^2))
            end
        end
        sky = [pix_to_world(truth, [1 + (nx - 1)rand(), 1 + (ny - 1)rand()]) for _ in 1:400]
        many = Table(ra=[s[1] for s in sky], dec=[s[2] for s in sky])
        r = refine_wcs(damaged, header, many)
        @test !r.refined
        @test r.wcs === header
        @test refine_wcs(damaged, header, many; min_significance=0).refined   # what it accepted before
    end

    @testset "refine_wcs_set" begin
        Random.seed!(17)
        nx, ny = 300, 260
        base = WCSTransform(2; crpix=[150.0, 130.0], crval=[3.05, -10.74],
                            cdelt=[-0.257 / 3600, 0.257 / 3600], ctype=["RA---TAN", "DEC--TAN"])
        sky = [pix_to_world(base, [20 + 260rand(), 20 + 220rand()]) for _ in 1:40]
        stars = Table(ra=[s[1] for s in sky], dec=[s[2] for s in sky])
        # three dithered frames, 5 different catalogue stars bright enough in each
        # -- too few for any one frame's own fit
        truths = [WCSTransform(2; crpix=base.crpix .+ [3.0k, -2.0k], crval=base.crval, cdelt=base.cdelt,
                               ctype=base.ctype) for k in 1:3]
        images = map(1:3) do k
            image = 100.0 .+ 5.0 .* randn(ny, nx)
            for s in sky[5k-4:5k]
                x, y = world_to_pix(truths[k], [s...])
                for j in max(1, floor(Int, x) - 8):min(nx, ceil(Int, x) + 8),
                    i in max(1, floor(Int, y) - 8):min(ny, ceil(Int, y) + 8)
                    image[i, j] += 2000.0 * exp(-((j - x)^2 + (i - y)^2) / (2 * 1.6^2))
                end
            end
            image
        end
        off(k, d) = WCSTransform(2; crpix=truths[k].crpix .+ d, crval=base.crval, cdelt=base.cdelt, ctype=base.ctype)
        error_arcsec(w, k) = maximum(sky[1:15]) do s
            a = pix_to_world(w, world_to_pix(truths[k], [s...]))
            hypot((a[1] - s[1]) * cosd(s[2]), a[2] - s[2]) * 3600
        end

        # headers 80 px off (beyond refine_wcs's default search), all by the same amount
        headers = [off(k, [-80.0, 10.0]) for k in 1:3]
        @test all(k -> !refine_wcs(images[k], headers[k], stars).refined, 1:3)
        r = refine_wcs_set(images, headers, stars)
        @test all(x -> x.refined && x.joint && x.n_matches == 15, r)
        @test all(k -> error_arcsec(r[k].wcs, k) < 0.1, 1:3)

        # one header disagreeing with the others by 8 px: no joint fit, and each
        # frame's own offset corrected instead
        headers[3] = off(3, [-72.0, 10.0])
        r = refine_wcs_set(images, headers, stars)
        @test !any(x -> x.joint, r)
        @test all(x -> !x.refined && x.shift_only, r)
        @test all(k -> error_arcsec(r[k].wcs, k) < 0.3, 1:3)
    end

    @testset "_filter_tracklets" begin
        dets = [Table(x=[10, 50], y=[10, 50], flux=[100.0, 100.0], xcen=[10.0, 50.0], ycen=[10.0, 50.0]),
                Table(x=[12, 50], y=[10, 50], flux=[90.0, 100.0], xcen=[12.0, 50.0], ycen=[10.0, 50.0]),
                Table(x=[14, 50], y=[10, 50], flux=[400.0, 100.0], xcen=[14.0, 50.0], ycen=[10.0, 50.0])]
        ts = [0.0, 1.0, 2.0]
        moving = [(frame=k, x=Float64(8 + 2k), y=10.0) for k in 1:3]       # 2 px/day, flux 100->400
        static = [(frame=k, x=50.0, y=50.0) for k in 1:3]
        tracklets = [moving, static]
        F = AsteroidPipeline._filter_tracklets
        @test F(tracklets, dets, ts) === tracklets
        @test F(tracklets, dets, ts; min_speed=1.0) == [moving]
        @test F(tracklets, dets, ts; max_flux_ratio=2.5) == [static]  # moving varies 4.4x
        @test isempty(F(tracklets, dets, ts; min_speed=1.0, max_flux_ratio=2.5))
        bent = [(frame=1, x=50.0, y=50.0), (frame=2, x=52.0, y=55.0), (frame=3, x=54.0, y=50.0)]
        @test AsteroidPipeline._line_residual(moving, ts) < 1e-9                  # exactly straight
        # a fourth frame both miss: both partial, so the bent one (~2.4 px RMS) goes
        dets4, ts4 = vcat(dets, [dets[1]]), [ts; 3.0]
        @test F([moving, bent], dets4, ts4; max_residual=1.0) == [moving]
        @test F([bent], dets, ts; max_residual=1.0) == [bent]                    # complete: left alone

        # four frames, one of them faint (fluxes 121, 114, 47, 95: 2019 UA185's):
        # kept without that point, as a partial tracklet that must still be straight
        flux4 = [121.0, 114.0, 47.0, 95.0]
        dets_ua = [Table(x=[10 + 2k], y=[10], flux=[flux4[k]], xcen=[10.0 + 2k], ycen=[10.0]) for k in 1:4]
        ua = [(frame=k, x=10.0 + 2k, y=10.0) for k in 1:4]
        @test F([ua], dets_ua, ts4; max_flux_ratio=2.5, max_residual=1.0) == [ua[[1, 2, 4]]]
        @test isempty(F([ua], dets_ua, ts4; max_flux_ratio=2.5, max_residual=1.0, min_frames=4))
        # not when the rest isn't straight
        ua_bent = [ua[1], ua[2], ua[3], (frame=4, x=18.0, y=20.0)]
        dets_bent = vcat(dets_ua[1:3], [Table(x=[18], y=[20], flux=[95.0], xcen=[18.0], ycen=[20.0])])
        @test isempty(F([ua_bent], dets_bent, ts4; max_flux_ratio=2.5, max_residual=1.0))
    end

    @testset "link_candidates" begin
        frame1 = Table(x=[10.0, 50.0], y=[10.0, 20.0])
        frame2 = Table(x=[12.0, 55.0], y=[11.0, 20.0])
        frame3 = Table(x=[14.0, 60.0], y=[12.0, 20.0])
        timestamps = [0.0, 1.0, 2.0]

        tracklets = link_candidates([frame1, frame2, frame3], timestamps; match_radius=0.5)
        @test length(tracklets) == 2
        endpoints = sort([t[end].x for t in tracklets])
        @test endpoints ≈ [14.0, 60.0]

        @test isempty(link_candidates([frame1, frame2, frame3], timestamps; max_speed=0.1))
        @test_throws ArgumentError link_candidates([frame1, frame2], [0.0, 0.0])

        # velocity refit: frame 2's noisy position makes the 2-point seed miss
        # frame 5 by 2.33 px; only the least-squares refit recovers it.
        noisy1 = Table(x=[10.0], y=[10.0])
        noisy2 = Table(x=[15.5], y=[6.7])
        noisy3 = Table(x=[20.0], y=[4.0])
        noisy4 = Table(x=[25.0], y=[1.0])
        noisy5 = Table(x=[30.0], y=[-2.0])
        noisy_ts = [0.0, 1.0, 2.0, 3.0, 4.0]

        refit_tracklets = link_candidates([noisy1, noisy2, noisy3, noisy4, noisy5], noisy_ts;
                                           match_radius=2.0, min_frames=5)
        @test length(refit_tracklets) == 1
        @test length(refit_tracklets[1]) == 5
        @test refit_tracklets[1][end].frame == 5
        @test [p.frame for p in refit_tracklets[1]] == 1:5  # sorted by frame

        # If the refit's velocity exceeds max_speed, the tracklet keeps its seed
        # points instead of being dropped (seed ~4.66, refit ~5.83, max 5.0).
        under1 = Table(x=[10.0], y=[10.0])
        under2 = Table(x=[14.0], y=[7.6])
        under3 = Table(x=[20.0], y=[4.0])
        under_tracklets = link_candidates([under1, under2, under3], [0.0, 1.0, 2.0];
                                           match_radius=3.0, min_frames=1, max_speed=5.0)
        @test length(under_tracklets) == 1
        @test length(under_tracklets[1]) == 3
        @test under_tracklets[1][2].x == 14.0 && under_tracklets[1][2].y == 7.6  # unrefit seed point kept
    end

    @testset "link_candidates seeds from any frame pair; _merge_clusters" begin
        # an object missing from frame 1 (detector gap, masked over a star):
        # seeding only from frames 1-2 could never link it, even at min_frames=3
        ts = [0.0, 1.0, 2.0, 3.0]
        frames = [Table(x=[500], y=[500]),                      # frame 1: just an unrelated source
                  Table(x=[12], y=[20]), Table(x=[14], y=[20]), Table(x=[16], y=[20])]
        @test isempty(link_candidates(frames, ts; match_radius=0.5))           # all 4 required
        found = link_candidates(frames, ts; match_radius=0.5, min_frames=3)
        @test length(found) == 1
        @test [p.frame for p in only(found)] == [2, 3, 4]

        # a fast object's trail, fragmented into three peaks, becomes one
        # detection at their flux-weighted centre; a satellite-length chain
        # is dropped; a lone detection passes untouched
        dets = Table(x=[100, 108, 116, 400, 412, 424, 436, 448, 460, 472, 900],
                     y=fill(50, 11), peak=fill(10.0, 11),
                     flux=[100.0, 200.0, 100.0, fill(50.0, 7)..., 300.0], flux_err=fill(5.0, 11),
                     xcen=Float64[100, 108, 116, 400, 412, 424, 436, 448, 460, 472, 900], ycen=fill(50.0, 11),
                     sharpness=fill(0.6, 11))
        merged = AsteroidPipeline._merge_clusters(dets, 15.0, 60.0)
        @test length(merged) == 2
        trail = only(filter(d -> d.flux == 400.0, collect(merged)))
        @test trail.xcen ≈ 108.0 && trail.flux_err ≈ sqrt(3) * 5.0
        @test any(d -> d.xcen == 900.0 && d.flux == 300.0, merged)            # the lone one
        @test !any(d -> 390 <= d.xcen <= 480, merged)                          # 72 px chain: dropped
    end

    @testset "find_variable_sources" begin
        # columns: a constant star, a variable star, and a mover, which must
        # not be picked up as variable
        frame1 = Table(x=[50.0, 80.0, 10.0], y=[50.0, 80.0, 10.0],
                        flux=[1000.0, 1000.0, 1500.0], flux_err=[10.0, 10.0, 10.0])
        frame2 = Table(x=[50.0, 80.0, 15.0], y=[50.0, 80.0, 10.0],
                        flux=[1000.0, 2000.0, 1500.0], flux_err=[10.0, 10.0, 10.0])
        frame3 = Table(x=[50.0, 80.0, 20.0], y=[50.0, 80.0, 10.0],
                        flux=[1000.0, 800.0, 1500.0], flux_err=[10.0, 10.0, 10.0])
        frame4 = Table(x=[50.0, 80.0, 25.0], y=[50.0, 80.0, 10.0],
                        flux=[1000.0, 2200.0, 1500.0], flux_err=[10.0, 10.0, 10.0])
        frame5 = Table(x=[50.0, 80.0, 30.0], y=[50.0, 80.0, 10.0],
                        flux=[1000.0, 900.0, 1500.0], flux_err=[10.0, 10.0, 10.0])
        timestamps = [0.0, 1.0, 2.0, 3.0, 4.0]

        groups = find_variable_sources([frame1, frame2, frame3, frame4, frame5], timestamps;
                                        position_tolerance=1.0)
        @test length(groups) == 1
        @test length(groups[1]) == 5
        @test all(p.x == 80.0 && p.y == 80.0 for p in groups[1])

        constant_chi2, constant_dof = variability_chi2(fill(1000.0, 5), fill(10.0, 5))
        @test constant_chi2 / constant_dof < 3.0

        variable_chi2, variable_dof = variability_chi2([p.flux for p in groups[1]],
                                                         [p.flux_err for p in groups[1]])
        @test variable_chi2 / variable_dof > 3.0

        wcs = WCSTransform(2; crpix=[50.0, 50.0], crval=[150.0, 20.0],
                            cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
        calibrated = astrometric_calibrate(groups, fill(wcs, 5), timestamps)
        @test length(calibrated) == 5
        @test all(==(1), calibrated.id)
    end

    @testset "photometric_scale / normalized find_variable_sources" begin
        # 5 constant stars at 0.8x flux in frame 2 (a transparency change, as seen
        # on real ZTF frames) and one variable that doesn't follow the scaling
        star_x = [10.0, 20.0, 30.0, 40.0, 50.0]
        star_y = copy(star_x)
        frame1 = Table(x=vcat(star_x, 70.0), y=vcat(star_y, 70.0),
                        flux=vcat(fill(1000.0, 5), 2000.0), flux_err=vcat(fill(5.0, 5), 10.0))
        frame2 = Table(x=vcat(star_x, 70.0), y=vcat(star_y, 70.0),
                        flux=vcat(fill(800.0, 5), 2000.0), flux_err=vcat(fill(4.0, 5), 10.0))
        frame3 = Table(x=vcat(star_x, 70.0), y=vcat(star_y, 70.0),
                        flux=vcat(fill(1000.0, 5), 2000.0), flux_err=vcat(fill(5.0, 5), 10.0))
        timestamps = [0.0, 1.0, 2.0]

        scales = photometric_scale([frame1, frame2, frame3])
        @test scales[1] == 1.0
        @test scales[2]≈1.25 rtol=1e-6   # 1000/800
        @test scales[3]≈1.0 rtol=1e-6

        groups = find_variable_sources([frame1, frame2, frame3], timestamps)
        @test length(groups) == 1
        @test all(p.x == 70.0 && p.y == 70.0 for p in groups[1])

        # without normalization the 0.8x mismatch can flag rescaled constant stars
        # too, so disabling it must not find fewer candidates
        raw_groups = find_variable_sources([frame1, frame2, frame3], timestamps; normalize=false)
        @test length(raw_groups) >= length(groups)
    end

    @testset "photometric_scale min_stars guard" begin
        # only 2 matched stars per frame (below min_stars=5): the scale stays 1.0
        frame1 = Table(x=[10.0, 20.0], y=[10.0, 20.0], flux=[1000.0, 1000.0], flux_err=[5.0, 5.0])
        frame2 = Table(x=[10.0, 20.0], y=[10.0, 20.0], flux=[800.0, 800.0], flux_err=[4.0, 4.0])

        scales = @test_logs (:warn, r"too few") match_mode = :any photometric_scale([frame1, frame2])
        @test scales == [1.0, 1.0]
    end

    @testset "astrometry" begin
        wcs = WCSTransform(2; crpix=[500.0, 500.0], crval=[150.0, 20.0],
                            cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

        ra, dec = pix_to_sky(wcs, 500.0, 500.0)
        @test ra ≈ 150.0 atol=1e-9
        @test dec ≈ 20.0 atol=1e-9

        loaded = load_wcs(WCS.to_header(wcs))
        ra2, dec2 = pix_to_sky(loaded, 500.0, 500.0)
        @test ra2 ≈ ra atol=1e-9
        @test dec2 ≈ dec atol=1e-9

        # CNPIX1/CNPIX2 (a legacy IRAF keyword pair, seen on real PS1 headers) makes
        # wcslib reject the whole header as singular; see load_wcs.
        cnpix_cards = rpad("CNPIX1  =                    0", 80) * rpad("CNPIX2  =                    0", 80)
        header_with_cnpix = WCS.to_header(wcs) * cnpix_cards
        @test_throws "singular" WCS.from_header(header_with_cnpix)
        loaded_cnpix = load_wcs(header_with_cnpix)
        ra3, dec3 = pix_to_sky(loaded_cnpix, 500.0, 500.0)
        @test ra3 ≈ ra atol=1e-9
        @test dec3 ≈ dec atol=1e-9

        tracklets = [[(frame=1, x=500.0, y=500.0), (frame=2, x=501.0, y=500.0)]]
        timestamps = [2460000.5, 2460000.51]
        candidates = astrometric_calibrate(tracklets, [wcs, wcs], timestamps)

        @test length(candidates) == 2
        @test all(==(1), candidates.id)
        @test candidates[1].ra ≈ 150.0 atol=1e-9
        @test candidates[2].epoch == timestamps[2]
        @test candidates[2].ra != candidates[1].ra
    end

    @testset "reprojection round-trip" begin
        # Grid A -> B -> A should preserve a source's position to well
        # under a pixel, since B differs from A only by a few pixels of
        # shift (like real same-night ZTF dither, not a large rotation).
        nx, ny = 60, 60
        wcs_a = WCSTransform(2; crpix=[30.0, 30.0], crval=[150.0, 20.0],
                              cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
        wcs_b = WCSTransform(2; crpix=[33.5, 27.5], crval=[150.0, 20.0],
                              cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

        img = fill(10.0, nx, ny)
        img[30, 30] = 5000.0

        onto_b, _ = reproject((img, wcs_a), wcs_b; shape_out=(nx, ny))
        back_to_a, mask = reproject((onto_b, wcs_b), wcs_a; shape_out=(nx, ny))

        # NaN marks out-of-footprint pixels, and argmax over NaN returns the first NaN
        peak = argmax(ifelse.(mask, back_to_a, -Inf))
        @test peak == CartesianIndex(30, 30)
        @test mask[30, 30]
    end

    @testset "build_reference" begin
        nx, ny = 50, 50
        wcs = WCSTransform(2; crpix=[25.0, 25.0], crval=[150.0, 20.0],
                            cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

        Random.seed!(11)
        frames = []
        for k in 1:7
            img = 100.0 .+ 2.0 .* randn(nx, ny)
            img[25, 25] += 400.0            # static source, present in every frame
            img[10 + k, 40 - k] += 3000.0   # "moving" outlier at a different pixel each frame
            push!(frames, (image=img, wcs=wcs, magzp=25.0, sigma=2.0))
        end

        image, sigma, mask = build_reference(frames, wcs, (nx, ny))

        @test image[25, 25] > 300.0          # static source survives the median
        @test all(mask)
        # none of the per-frame outlier pixels should show up at ~3000
        # above background in the median-combined reference
        outlier_pixels = [image[10+k, 40-k] for k in 1:7]
        @test all(v -> v < 200.0, outlier_pixels)
    end

    @testset "build_reference with workers" begin
        # The `workers` path must reproduce the sequential result exactly.
        nx, ny = 50, 50
        wcs = WCSTransform(2; crpix=[25.0, 25.0], crval=[150.0, 20.0],
                            cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

        Random.seed!(11)
        frames = []
        for k in 1:7
            img = 100.0 .+ 2.0 .* randn(nx, ny)
            img[25, 25] += 400.0
            img[10 + k, 40 - k] += 3000.0
            push!(frames, (image=img, wcs=wcs, magzp=25.0, sigma=2.0))
        end

        sequential = build_reference(frames, wcs, (nx, ny))

        new_workers = addprocs(2; exeflags="--project=$(Base.active_project())")
        try
            @everywhere new_workers using AsteroidPipeline
            distributed = build_reference(frames, wcs, (nx, ny); workers=new_workers)
            @test distributed[1] ≈ sequential[1]
            @test distributed[2] ≈ sequential[2]
            @test distributed[3] == sequential[3]
        finally
            rmprocs(new_workers)
        end
    end

    @testset "estimate_psf" begin
        Random.seed!(12)
        nx, ny = 120, 120
        img = 100.0 .+ 3.0 .* randn(nx, ny)
        true_sigma = 1.8
        centers = [(30, 40), (90, 30), (60, 90), (100, 100)]
        for (cx, cy) in centers, i in 1:nx, j in 1:ny
            r2 = (i - cx)^2 + (j - cy)^2
            img[i, j] += 4000.0 * exp(-r2 / (2 * true_sigma^2))
        end

        psf = estimate_psf(img; threshold=15.0, min_separation=20.0)
        @test sum(psf) ≈ 1.0 atol=1e-9

        c = size(psf, 1) ÷ 2 + 1
        @test argmax(psf) == CartesianIndex(c, c)
        # empirical FWHM should match the injected sigma * 2sqrt(2ln2); findfirst
        # over 0:(c-1) returns a position, one more than the radius
        half_max_pos = findfirst(r -> psf[c+r, c] < psf[c, c] / 2, 0:(c - 1))
        @test half_max_pos !== nothing
        half_max_r = half_max_pos - 1
        @test half_max_r - 1 <= true_sigma * 2sqrt(2log(2)) / 2 <= half_max_r + 1
    end

    @testset "fit_moffat_psf / estimate_psf fallback" begin
        Random.seed!(13)
        nx, ny = 120, 120
        img = 100.0 .+ 3.0 .* randn(nx, ny)
        true_alpha, true_beta = 3.0, 2.5
        # centers ~25 px apart, under min_separation=40: every star fails isolation,
        # forcing the analytic fallback (relaxation_attempts=0 keeps it that way)
        centers = [(50, 50), (70, 55), (55, 70), (75, 75)]
        for (cx, cy) in centers, i in 1:nx, j in 1:ny
            r2 = (i - cx)^2 + (j - cy)^2
            img[i, j] += 3000.0 / (1 + r2 / true_alpha^2)^true_beta
        end

        @test_throws ErrorException estimate_psf(img; threshold=15.0, fallback=false,
                                                   relaxation_attempts=0)

        psf = estimate_psf(img; threshold=15.0, relaxation_attempts=0)  # fallback=true by default
        @test sum(psf) ≈ 1.0 atol=1e-6

        c = size(psf, 1) ÷ 2 + 1
        @test argmax(psf) == CartesianIndex(c, c)

        true_fwhm = 2 * true_alpha * sqrt(2^(1 / true_beta) - 1)
        half_max_pos = findfirst(r -> psf[c+r, c] < psf[c, c] / 2, 0:(c - 1))
        @test half_max_pos !== nothing
        half_max_r = half_max_pos - 1
        @test half_max_r - 1 <= true_fwhm / 2 <= half_max_r + 1

        empty_image = fill(100.0, 60, 60)
        @test_throws ErrorException fit_moffat_psf(empty_image; threshold=15.0)
    end

    @testset "estimate_psf min_separation relaxation" begin
        Random.seed!(14)
        nx, ny = 120, 120
        img = 100.0 .+ 3.0 .* randn(nx, ny)
        true_sigma = 1.8
        # 25 px apart: fails min_separation=40 but passes the first relaxed attempt
        # (20), so a real empirical PSF is recovered
        centers = [(40, 60), (65, 60)]
        for (cx, cy) in centers, i in 1:nx, j in 1:ny
            r2 = (i - cx)^2 + (j - cy)^2
            img[i, j] += 4000.0 * exp(-r2 / (2 * true_sigma^2))
        end

        # relaxation_attempts=0 must not relax at all: with fallback=false
        # (so this can't succeed via fit_moffat_psf either), it has to
        # fail cleanly at the original min_separation=40.
        @test_throws ErrorException estimate_psf(img; threshold=15.0, fallback=false,
                                                   min_separation=40.0, relaxation_attempts=0)

        # With relaxation (the default, relaxation_attempts=2), the same
        # data succeeds via the empirical path instead of falling back —
        # verified by the recovered FWHM matching the injected Gaussian's.
        psf_relaxed = estimate_psf(img; threshold=15.0, min_separation=40.0)
        @test sum(psf_relaxed) ≈ 1.0 atol=1e-9
        c = size(psf_relaxed, 1) ÷ 2 + 1
        @test argmax(psf_relaxed) == CartesianIndex(c, c)
        half_max_pos = findfirst(r -> psf_relaxed[c+r, c] < psf_relaxed[c, c] / 2, 0:(c - 1))
        @test half_max_pos !== nothing
        half_max_r = half_max_pos - 1
        @test half_max_r - 1 <= true_sigma * 2sqrt(2log(2)) / 2 <= half_max_r + 1
    end

    @testset "zogy_subtract" begin
        psf = zeros(9, 9)
        for i in 1:9, j in 1:9
            r2 = (i - 5)^2 + (j - 5)^2
            psf[i, j] = exp(-r2 / (2 * 1.5^2))
        end
        psf ./= sum(psf)
        sz = (96, 96)

        # N == R, identical PSF, no noise -> D ~ 0 to machine precision
        img = fill(100.0, sz)
        img[60:68, 60:68] .+= 5000.0 .* psf
        D, _ = zogy_subtract(img, img; psf_n=psf, psf_r=psf, sigma_n=1.0, sigma_r=1.0)
        @test maximum(abs.(D)) < 1e-8

        # pure noise -> std(S_corr) ~ 1; the high gain keeps _variance_map's Poisson
        # term, which noise excursions would inflate, out of this test
        Random.seed!(13)
        sigma = 5.0
        n_noise = 100.0 .+ sigma .* randn(sz)
        r_noise = 100.0 .+ sigma .* randn(sz)
        _, s_corr_noise = zogy_subtract(n_noise, r_noise; psf_n=psf, psf_r=psf,
                                         sigma_n=sigma, sigma_r=sigma, gain_n=1e8, gain_r=1e8)
        @test isapprox(std(s_corr_noise), 1.0; atol=0.15)
        @test isapprox(mean(s_corr_noise), 0.0; atol=0.15)

        # injected point source of known flux -> peak S_corr matches the
        # analytic matched-filter SNR, F0 * sqrt(sum(psf^2) / (sigma_n^2 + sigma_r^2))
        Random.seed!(14)
        flux = 2000.0
        n_img = 100.0 .+ sigma .* randn(sz)
        r_img = 100.0 .+ sigma .* randn(sz)
        n_img[60:68, 60:68] .+= flux .* psf
        _, s_corr_src = zogy_subtract(n_img, r_img; psf_n=psf, psf_r=psf,
                                       sigma_n=sigma, sigma_r=sigma, gain_n=1e8, gain_r=1e8)
        predicted_snr = flux * sqrt(sum(psf .^ 2) / (2 * sigma^2))
        @test isapprox(maximum(s_corr_src), predicted_snr; rtol=0.05)
        @test argmax(s_corr_src) == CartesianIndex(64, 64)

        # asymmetric PSFs and noise levels: with psf_n == psf_r most of the general
        # terms in _k_n/_k_r cancel, so the symmetric checks miss them
        gaussian_psf(s) = (p = [exp(-((i - 5)^2 + (j - 5)^2) / (2s^2)) for i in 1:9, j in 1:9]; p ./ sum(p))
        psf_n2, psf_r2 = gaussian_psf(1.5), gaussian_psf(2.2)
        sigma_n2, sigma_r2 = 5.0, 2.0
        Random.seed!(15)
        n_noise2 = 100.0 .+ sigma_n2 .* randn(sz)
        r_noise2 = 100.0 .+ sigma_r2 .* randn(sz)
        _, s_corr_asym = zogy_subtract(n_noise2, r_noise2; psf_n=psf_n2, psf_r=psf_r2,
                                        sigma_n=sigma_n2, sigma_r=sigma_r2, gain_n=1e8, gain_r=1e8)
        @test isapprox(std(s_corr_asym), 1.0; atol=0.15)

        # mismatched sky background between N and R: without subtracting each image's
        # own background, the FFT's DC term swamped S_corr (seen on real ZTF data)
        Random.seed!(16)
        n_bg_noise = 400.0 .+ sigma .* randn(sz)   # science: brighter sky
        r_bg_noise = 100.0 .+ sigma .* randn(sz)   # reference: darker sky
        _, s_corr_bg = zogy_subtract(n_bg_noise, r_bg_noise; psf_n=psf, psf_r=psf,
                                      sigma_n=sigma, sigma_r=sigma, gain_n=1e8, gain_r=1e8)
        @test isapprox(mean(s_corr_bg), 0.0; atol=0.2)
        @test isapprox(std(s_corr_bg), 1.0; atol=0.15)

        # regression test: omitting n_sources/r_sources (silently dropping
        # V_ast) used to be undiscoverable without reading the docstring —
        # now warns explicitly, one call with, one without.
        sources_n = Table(x=[64.0], y=[64.0])
        sources_r = Table(x=[64.2], y=[64.1])
        @test_logs (:warn, r"V_ast") match_mode=:any zogy_subtract(
            img, img; psf_n=psf, psf_r=psf, sigma_n=1.0, sigma_r=1.0)
        @test_logs zogy_subtract(img, img; psf_n=psf, psf_r=psf, sigma_n=1.0, sigma_r=1.0,
                                  n_sources=sources_n, r_sources=sources_r)
    end

    @testset "run_pipeline with reference" begin
        mktempdir() do dir
            nx, ny = 100, 100
            wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
            header(mjd) = FITSHeader(
                ["MJD-OBS", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                 "CTYPE1", "CTYPE2", "GAIN", "MAGZP"],
                Any[mjd, wcs.crpix[1], wcs.crpix[2], wcs.crval[1], wcs.crval[2],
                    wcs.cdelt[1], wcs.cdelt[2], wcs.ctype[1], wcs.ctype[2], 1.0, 25.0],
                fill("", 11),
            )
            star!(img, cx, cy, flux; sigma=1.8) = for i in 1:size(img, 1), j in 1:size(img, 2)
                r2 = (i - cx)^2 + (j - cy)^2
                img[i, j] += flux * exp(-r2 / (2 * sigma^2))
            end

            Random.seed!(15)
            # deep reference: 6 frames, only ever the bright static star,
            # so the moving source is only ever seen against a clean stack
            refpaths = String[]
            for k in 1:6
                img = 100.0 .+ 3.0 .* randn(nx, ny)
                star!(img, 20.0, 80.0, 2000.0)
                path = joinpath(dir, "ref$k.fits")
                FITS(path, "w") do f
                    write(f, img; header=header(59000.0 + k))
                end
                push!(refpaths, path)
            end

            # science triplet: the static star plus a mover clear of it, each frame
            # dithered a few pixels, so reprojection leaves a NaN border to sanitize
            dither(mjd, dx, dy) = FITSHeader(
                ["MJD-OBS", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                 "CTYPE1", "CTYPE2", "GAIN", "MAGZP"],
                Any[mjd, wcs.crpix[1] + dx, wcs.crpix[2] + dy, wcs.crval[1], wcs.crval[2],
                    wcs.cdelt[1], wcs.cdelt[2], wcs.ctype[1], wcs.ctype[2], 1.0, 25.0],
                fill("", 11),
            )
            x0, y0, dx, dy = 70.0, 20.0, 6.0, -4.0
            mjd0 = 60000.0
            scipaths = String[]
            for k in 0:2
                img = 100.0 .+ 3.0 .* randn(nx, ny)
                star!(img, 20.0, 80.0, 2000.0)
                star!(img, x0 + k * dx, y0 + k * dy, 600.0)
                path = joinpath(dir, "sci$k.fits")
                FITS(path, "w") do f
                    write(f, img; header=dither(mjd0 + k * 1e-3, 4.0k, -3.0k))
                end
                push!(scipaths, path)
            end

            sci_wcs = load_frame(scipaths[1]).wcs
            refs = [load_frame(p) for p in refpaths]
            image, sigma, mask = build_reference(refs, sci_wcs, (nx, ny))
            psf = estimate_psf(image; threshold=20.0, min_separation=15.0)
            reference = (image=image, sigma=sigma, mask=mask, psf=psf, wcs=sci_wcs)

            # quality_max_std=Inf: on a 100x100 image the gate's plain std is skewed by
            # the injected source (see "run_pipeline quality gate")
            candidates = run_pipeline(scipaths; threshold=6.0, match_radius=2.0,
                                       reference=reference,
                                       psf_threshold=20.0, psf_min_separation=15.0,
                                       quality_max_std=Inf)

            @test length(unique(candidates.id)) >= 1
            @test !any(isnan, candidates.x) && !any(isnan, candidates.ra)
            by_frame = Dict(row.frame => row for row in candidates if row.id == candidates[1].id)
            @test haskey(by_frame, 1) && haskey(by_frame, 3)
            if haskey(by_frame, 1) && haskey(by_frame, 3)
                # expected positions on reference.wcs's grid, corrected for each frame's
                # CRPIX dither
                @test by_frame[1].x ≈ x0 atol=1.5
                @test by_frame[1].y ≈ y0 atol=1.5
                @test by_frame[3].x ≈ x0 + 2dx - 4.0 * 2 atol=1.5
                @test by_frame[3].y ≈ y0 + 2dy - (-3.0 * 2) atol=1.5
            end
        end
    end

    @testset "run_pipeline quality gate" begin
        mktempdir() do dir
            # 600x600: the gate's plain std over a small image is skewed by any bright
            # source (a ~60 sigma one pushed a clean 100x100 frame to ~4, here ~1.1)
            nx, ny = 600, 600
            wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
            header(mjd) = FITSHeader(
                ["MJD-OBS", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                 "CTYPE1", "CTYPE2", "GAIN", "MAGZP"],
                Any[mjd, wcs.crpix[1], wcs.crpix[2], wcs.crval[1], wcs.crval[2],
                    wcs.cdelt[1], wcs.cdelt[2], wcs.ctype[1], wcs.ctype[2], 1.0, 25.0],
                fill("", 11),
            )
            star!(img, cx, cy, flux; sigma=1.8) = for i in 1:size(img, 1), j in 1:size(img, 2)
                r2 = (i - cx)^2 + (j - cy)^2
                img[i, j] += flux * exp(-r2 / (2 * sigma^2))
            end

            Random.seed!(17)
            refpaths = String[]
            for k in 1:6
                img = 100.0 .+ 3.0 .* randn(nx, ny)
                star!(img, 300.0, 300.0, 2000.0)
                path = joinpath(dir, "ref$k.fits")
                FITS(path, "w") do f
                    write(f, img; header=header(59000.0 + k))
                end
                push!(refpaths, path)
            end

            # Frames 1-2: clean, plus a bright transient that must not trip the gate.
            # Frame 3: a broad low-contrast glow, which sigma_n doesn't correct for.
            x0, y0 = 450.0, 150.0
            scipaths = String[]
            for k in 1:3
                img = 100.0 .+ 3.0 .* randn(nx, ny)
                star!(img, 300.0, 300.0, 2000.0)
                star!(img, x0, y0, 250.0)
                k == 3 && star!(img, 150.0, 450.0, 60.0; sigma=20.0)
                path = joinpath(dir, "sci$k.fits")
                FITS(path, "w") do f
                    write(f, img; header=header(60000.0 + (k - 1) * 1e-3))
                end
                push!(scipaths, path)
            end

            sci_wcs = load_frame(scipaths[1]).wcs
            refs = [load_frame(p) for p in refpaths]
            image, sigma, mask = build_reference(refs, sci_wcs, (nx, ny))
            psf = estimate_psf(image; threshold=20.0, min_separation=15.0)
            reference = (image=image, sigma=sigma, mask=mask, psf=psf, wcs=sci_wcs)

            candidates = @test_logs (:warn, r"quality_max_std") match_mode = :any run_pipeline(
                scipaths; threshold=6.0, match_radius=5.0, min_frames=1,
                reference=reference, psf_threshold=20.0, psf_min_separation=15.0)

            @test !isempty(candidates)
            @test !(3 in candidates.frame)  # the glow-patch frame contributed nothing
            @test 1 in candidates.frame     # the real transient still recovered elsewhere

            # min_frames omitted: its default subtracts the gated frame, giving 2, and
            # must match the explicit min_frames=1 result
            auto_candidates = @test_logs (:warn, r"quality_max_std") match_mode = :any run_pipeline(
                scipaths; threshold=6.0, match_radius=5.0,
                reference=reference, psf_threshold=20.0, psf_min_separation=15.0)
            @test !isempty(auto_candidates)
            @test !(3 in auto_candidates.frame)
            @test 1 in auto_candidates.frame

            # lowering the threshold further should also drop the clean frames
            @test_logs((:warn, r"quality_max_std"), (:warn, r"quality_max_std"),
                       (:warn, r"quality_max_std"), match_mode = :any,
                       run_pipeline(scipaths; threshold=6.0, match_radius=5.0, min_frames=1,
                                    reference=reference, psf_threshold=20.0, psf_min_separation=15.0,
                                    quality_max_std=0.0))
        end
    end

    @testset "run_pipeline" begin
        mktempdir() do dir
            nx, ny = 80, 60  # deliberately non-square, to catch axis-order bugs
            wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

            # true asteroid track, in FITS pixel coordinates (x, y differ so a
            # transposed axis order would be caught)
            x0, y0, dx, dy = 20.0, 45.0, 5.0, -3.0
            mjd0 = 60000.0

            Random.seed!(2)
            paths = String[]
            for k in 0:2
                raw = 100.0 .+ 5.0 .* randn(nx, ny)  # FITS-native (x, y) layout
                xk, yk = x0 + k * dx, y0 + k * dy
                for i in 1:nx, j in 1:ny
                    r2 = (i - xk)^2 + (j - yk)^2
                    raw[i, j] += 500.0 * exp(-r2 / (2 * 2.0^2))
                end

                header = FITSHeader(
                    ["MJD-OBS", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                     "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"],
                    Any[mjd0 + k * 1e-3, wcs.crpix[1], wcs.crpix[2], wcs.crval[1], wcs.crval[2],
                        wcs.cdelt[1], wcs.cdelt[2], wcs.ctype[1], wcs.ctype[2]],
                    fill("", 9),
                )

                path = joinpath(dir, "frame$k.fits")
                FITS(path, "w") do f
                    write(f, raw; header=header)
                end
                push!(paths, path)
            end

            candidates = run_pipeline(paths; threshold=5.0, match_radius=1.0)

            @test length(candidates) == 3
            @test length(unique(candidates.id)) == 1
            @test sort(candidates.frame) == [1, 2, 3]

            by_frame = Dict(row.frame => row for row in candidates)
            @test by_frame[1].x ≈ x0 atol=1.0
            @test by_frame[1].y ≈ y0 atol=1.0
            @test by_frame[3].x ≈ x0 + 2dx atol=1.0
            @test by_frame[3].y ≈ y0 + 2dy atol=1.0
        end
    end

    @testset "run_pipeline links across dithered frames" begin
        # Pointing offsets between frames (up to 9 px, over match_radius) must not
        # bend a straight track: frames are linked on a common grid.
        mktempdir() do dir
            nx, ny = 120, 90
            crval = [150.0, 20.0]; cdelt = [-1 / 3600, 1 / 3600]
            truth = WCSTransform(2; crpix=[60.0, 45.0], crval, cdelt, ctype=["RA---TAN", "DEC--TAN"])
            dithers = [(0, 0), (6, 9), (-5, 4), (3, -8)]
            track = [(30.0 + 12k, 60.0 - 7k) for k in 0:3]     # straight line, truth grid
            mjd0 = 60000.0
            Random.seed!(14)
            paths = String[]
            for k in 1:4
                dx, dy = dithers[k]
                raw = 100.0 .+ 5.0 .* randn(nx, ny)            # FITS-native (x, y)
                xk, yk = track[k][1] + dx, track[k][2] + dy
                for i in 1:nx, j in 1:ny
                    raw[i, j] += 600.0 * exp(-((i - xk)^2 + (j - yk)^2) / (2 * 1.8^2))
                end
                header = FITSHeader(
                    ["MJD-OBS", "EXPTIME", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"],
                    Any[mjd0 + (k - 1) * 1e-2, 60.0, truth.crpix[1] + dx, truth.crpix[2] + dy, crval[1], crval[2],
                        cdelt[1], cdelt[2], "RA---TAN", "DEC--TAN"],
                    fill("", 10))
                path = joinpath(dir, "frame$k.fits")
                FITS(f -> write(f, raw; header=header), path, "w")
                push!(paths, path)
            end

            candidates = run_pipeline(paths; threshold=5.0, match_radius=3.0)
            @test length(unique(candidates.id)) == 1
            @test sort(candidates.frame) == [1, 2, 3, 4]
            for row in candidates
                ra, dec = pix_to_world(truth, collect(track[row.frame]))
                @test hypot((row.ra - ra) * cosd(dec), row.dec - dec) * 3600 < 0.5
                # epochs are mid-exposure (EXPTIME=60 s -> +30 s), as MPC reports need
                @test row.epoch ≈ mjd0 + (row.frame - 1) * 1e-2 + 30 / 86400 + 2400000.5 atol=1e-7
            end
            starts = run_pipeline(paths; threshold=5.0, match_radius=3.0, exptime_key=nothing)
            @test first(starts).epoch ≈ first(candidates).epoch - 30 / 86400 atol=1e-7
        end
    end

    @testset "search_field" begin
        mktempdir() do dir
            # 300x225: on a smaller canvas a variable this bright corrupts the
            # background estimate (hundreds of spurious detections at 80x60)
            nx, ny = 300, 225
            wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

            # true asteroid track (as in the "run_pipeline" testset above)
            x0, y0, dx, dy = 20.0, 45.0, 5.0, -3.0
            # stationary variable star clear of the track, with a flux swing far beyond
            # photon noise and S/N above find_variable_sources's floor
            xv, yv = 220.0, 150.0
            var_amps = [3000.0, 8000.0, 3000.0]
            mjd0 = 60000.0

            Random.seed!(3)
            paths = String[]
            for k in 0:2
                raw = 100.0 .+ 5.0 .* randn(nx, ny)
                xk, yk = x0 + k * dx, y0 + k * dy
                for i in 1:nx, j in 1:ny
                    r2 = (i - xk)^2 + (j - yk)^2
                    raw[i, j] += 500.0 * exp(-r2 / (2 * 2.0^2))
                    r2v = (i - xv)^2 + (j - yv)^2
                    raw[i, j] += var_amps[k+1] * exp(-r2v / (2 * 2.0^2))
                end

                header = FITSHeader(
                    ["MJD-OBS", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                     "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"],
                    Any[mjd0 + k * 1e-3, wcs.crpix[1], wcs.crpix[2], wcs.crval[1], wcs.crval[2],
                        wcs.cdelt[1], wcs.cdelt[2], wcs.ctype[1], wcs.ctype[2]],
                    fill("", 9),
                )

                path = joinpath(dir, "sf_frame$k.fits")
                FITS(path, "w") do f
                    write(f, raw; header=header)
                end
                push!(paths, path)
            end

            result = search_field(paths; threshold=5.0, match_radius=1.0,
                                   variability_position_tolerance=1.0)
            baseline = run_pipeline(paths; threshold=5.0, match_radius=1.0)

            # one shared _detect_all_frames pass must produce exactly the
            # same movers table run_pipeline computes on its own
            @test result.movers.id == baseline.id
            @test result.movers.frame == baseline.frame
            @test result.movers.x == baseline.x
            @test result.movers.y == baseline.y
            @test result.movers.ra == baseline.ra
            @test result.movers.dec == baseline.dec
            @test result.movers.epoch == baseline.epoch

            mover_ids = [row.id for row in result.movers
                         if row.frame == 1 && isapprox(row.x, x0; atol=1.0) && isapprox(row.y, y0; atol=1.0)]
            mover_id = only(mover_ids)
            mover_rows = Dict(row.frame => row for row in result.movers if row.id == mover_id)
            @test mover_rows[3].x ≈ x0 + 2dx atol=1.0
            @test mover_rows[3].y ≈ y0 + 2dy atol=1.0

            @test length(result.variables) == 3
            @test all(row -> isapprox(row.x, xv; atol=1.0) && isapprox(row.y, yv; atol=1.0), result.variables)
        end
    end

    @testset "light_curve / recover_rotation_period" begin
        mktempdir() do dir
            nx, ny = 60, 60
            wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
            ra, dec = pix_to_sky(wcs, 30.0, 30.0)

            true_period = 0.2  # days
            F0, amp = 3000.0, 0.4
            Random.seed!(21)

            paths = String[]
            mjd0 = 60000.0
            n = 40
            for k in 0:n-1
                t = mjd0 + k * 0.011 + 0.003 * rand()  # uneven sampling over ~2 periods
                flux_k = F0 * (1 + amp * sin(2pi * t / true_period))
                img = 100.0 .+ 3.0 .* randn(nx, ny)
                for i in 1:nx, j in 1:ny
                    r2 = (i - 30)^2 + (j - 30)^2
                    img[i, j] += flux_k * exp(-r2 / (2 * 1.8^2)) / (2pi * 1.8^2)
                end
                header = FITSHeader(
                    ["MJD-OBS", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                     "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"],
                    Any[t, wcs.crpix[1], wcs.crpix[2], wcs.crval[1], wcs.crval[2],
                        wcs.cdelt[1], wcs.cdelt[2], wcs.ctype[1], wcs.ctype[2]],
                    fill("", 9),
                )
                path = joinpath(dir, "f$k.fits")
                FITS(path, "w") do f
                    write(f, img; header=header)
                end
                push!(paths, path)
            end

            times, flux, flux_err = light_curve(paths, ra, dec; aperture_radius=6.0)
            @test length(times) == n
            @test all(>(0), flux_err)
            @test maximum(flux) > minimum(flux)  # the injected variability is present

            result = recover_rotation_period(times, flux; minimum_period=0.05, maximum_period=1.0)
            @test result.period ≈ true_period rtol=0.05
            @test result.false_alarm_probability < 1e-3

            @test_throws ArgumentError recover_rotation_period(times, flux;
                                                                 minimum_period=1.0, maximum_period=0.5)
        end
    end

    @testset "SkyBoT response status / retry" begin
        # No network: response shapes captured from the live service. Flag -1 is a
        # server-side crash behind HTTP 200, which must not parse as no matches.
        header = "# Flag: 1\n# Ticket: 1\n# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec) \n"
        found = header * " - | 2018 LT | 00 12 10.1460 | -10 46 56.759 | Mars-Crosser | 20.0 | 2.974 | 140.155 \n"
        none = "# Flag: 0\n# Ticket: 2\nNo solar system object was found in the requested FOV: ...\n"
        crashed = "# Flag: -1\n# Ticket: 3\nSkyBoT asteroid conesearch -> numIntAstromJ2000_para[1]: computeEphemeris: \n" *
                  "Program received signal SIGBUS: Access to an undefined portion of a memory object.\n"

        @test only(AsteroidPipeline._parse_skybot(found)).name == "2018 LT"
        @test isempty(AsteroidPipeline._parse_skybot(none))
        @test_throws SkyBoTServiceError AsteroidPipeline._parse_skybot(crashed)
        @test_throws SkyBoTServiceError AsteroidPipeline._parse_skybot("")
        @test_throws SkyBoTServiceError AsteroidPipeline._parse_skybot("<html>502 Bad Gateway</html>")

        c = (id=1, ra=3.0492, dec=-10.7441, epoch=2460934.9349)
        # transient: two crashes, then a real answer — retried, not dropped
        responses = [crashed, crashed, found]
        calls = Ref(0)
        flaky = _ -> (calls[] += 1; responses[calls[]])
        matches = AsteroidPipeline._skybot_matches(c, 0.12; fetch=flaky, backoff_seconds=0)
        @test calls[] == 3
        @test only(matches).name == "2018 LT"

        # persistent: gives up after _SKYBOT_ATTEMPTS and throws, never empty
        calls[] = 0
        always_down = _ -> (calls[] += 1; crashed)
        @test_throws SkyBoTServiceError AsteroidPipeline._skybot_matches(c, 0.12; fetch=always_down, backoff_seconds=0)
        @test calls[] == AsteroidPipeline._SKYBOT_ATTEMPTS

        # observatory -> topocentric ephemerides (SkyBoT's -loc); unset -> geocentric
        sent = Ref{Dict}()
        capture = q -> (sent[] = q; found)
        AsteroidPipeline._skybot_matches(c, 0.12; observatory="F51", fetch=capture)
        @test sent[]["-observer"] == "F51"
        AsteroidPipeline._skybot_matches(c, 0.12; fetch=capture)
        @test !haskey(sent[], "-observer")
        # an empty field is HTTP 204 with no body; an empty 200 is still a failure
        @test isempty(AsteroidPipeline._parse_skybot(AsteroidPipeline._skybot_text(HTTP.Response(204))))
        @test_throws SkyBoTServiceError AsteroidPipeline._parse_skybot(AsteroidPipeline._skybot_text(HTTP.Response(200, "")))
        # a small cone is asked of SkyBoT as 2′ (it loses objects in tiny ones),
        # and what lies beyond the radius asked for is dropped here
        @test isempty(AsteroidPipeline._skybot_matches(c, 5 / 3600; fetch=capture))      # 2018 LT is 140″ off
        @test parse(Float64, sent[]["-rd"]) ≈ 2 / 60
        @test only(AsteroidPipeline._skybot_matches(c, 150 / 3600; fetch=capture)).name == "2018 LT"
        # a request failing inside @sync/@async comes out as itself, which
        # callers treat as SkyBoT being unavailable, not wrapped
        wrapped = try
            @sync @async throw(SkyBoTServiceError("down"))
        catch e
            e
        end
        @test wrapped isa CompositeException
        @test AsteroidPipeline._task_root(wrapped) isa SkyBoTServiceError
        @test_throws ArgumentError crossmatch_catalog([c], :vsx; radius=5.0, observatory="F51")
        @test_throws ArgumentError crossmatch_catalog([c], :skybot; radius=5.0, observatory="Palomar")
    end

    @testset "JPL Small-Body Identification, while SkyBoT is unavailable" begin
        A = AsteroidPipeline
        @test A._jpl_sexagesimal(3.0416 / 15, "-") == "00-12-09.98"
        @test A._jpl_sexagesimal(-10.7814, "M") == "M10-46-53.04"
        @test A._jpl_name("610588 (2006 CX81)") == "2006 CX81"
        @test A._jpl_name("(2018 LT)") == "2018 LT"
        @test A._jpl_name("433 Eros (A898 PA)") == "Eros"
        @test A._jpl_name("C/2020 F3 (NEOWISE)") == "C/2020 F3 (NEOWISE)"

        # the answer's shape, from the live service (2025-09-16 10:02 UTC, F51)
        answer = """{"fields_second":["Object name","Astrometric RA (hh:mm:ss)","Astrometric Dec (dd mm'ss\\")",
            "Dist. from center RA (\\")","Dist. from center Dec (\\")","Dist. from center Norm (\\")",
            "Visual magnitude (V)","RA rate (\\"/h)","Dec rate (\\"/h)"],"n_second_pass":2,
            "data_second_pass":[["610588 (2006 CX81)","00:12:33.19","-10 49'06.9\\"","342.","-134.","367.0","21.4","-2.507E+01","-1.950E+01"],
                                ["(2018 LT)","00:12:09.98","-10 46'52.9\\"","-.340",".111",".3579","20.0","-2.249E+01","-3.664E+01"]]}"""
        objects = A._parse_jpl_sbident(answer)
        @test [o.name for o in objects] == ["2006 CX81", "2018 LT"]
        lt = objects[2]
        @test lt.ra ≈ 15 * (12 / 60 + 9.98 / 3600) && lt.dec ≈ -(10 + 46 / 60 + 52.9 / 3600)
        @test (lt.mv, lt.ra_rate, lt.dec_rate) == (20.0, -22.49, -36.64)
        @test isempty(A._parse_jpl_sbident("""{"signature":{},"n_second_pass":0}"""))
        @test_throws ErrorException A._parse_jpl_sbident("""{"message":"invalid MPC code"}""")

        # with SkyBoT marked down, a search is answered from a fetched field,
        # its objects moved along their rates: 2018 LT 48 min later, per Horizons
        t0 = 2460934.918246
        down, fields = A._SKYBOT_DOWN_UNTIL[], copy(A._JPL_FIELDS)
        try
            A._SKYBOT_DOWN_UNTIL[] = time() + 600
            empty!(A._JPL_FIELDS)
            push!(A._JPL_FIELDS, (ra=lt.ra, dec=lt.dec, half_width=0.5, epoch=t0, observatory="F51", objects=objects))
            later = (ra=3.036479058, dec=-10.789470924, epoch=2460934.951521)
            m = only(A._known_objects_near(later, 5 / 3600; observatory="F51"))
            @test m.name == "2018 LT" && m.source == "JPL"
            @test m.distance_arcsec < 0.5
        finally
            A._SKYBOT_DOWN_UNTIL[] = down
            empty!(A._JPL_FIELDS); append!(A._JPL_FIELDS, fields)
        end
    end

    @testset "gaia_reference_stars" begin
        @test_throws ArgumentError gaia_reference_stars(3.05, -10.74, 10.0; epoch_jd=2460934.9)
        try
            # The real IASC practice field (XY54_p10) refine_wcs was
            # validated on: VizieR returned 144 Gaia DR3 stars in 0.13 deg.
            stars = gaia_reference_stars(3.0492, -10.7441, 0.05; epoch_jd=2460934.9)
            @test length(stars) > 5
            @test all(s -> hypot((s.ra - 3.0492) * cosd(-10.7441), s.dec + 10.7441) < 0.051, stars)
            @test all(s -> 12.0 <= s.gmag <= 21.0, stars)
            # proper motion is applied: the same stars a decade earlier sit elsewhere
            then = gaia_reference_stars(3.0492, -10.7441, 0.05; epoch_jd=2457389.0)
            @test length(then) == length(stars)
            @test any(i -> stars.ra[i] != then.ra[i] || stars.dec[i] != then.dec[i], eachindex(stars))
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end
    end

    @testset "frame_epoch" begin
        # Pan-STARRS1-shaped header: DATE-OBS is UTC, MJD-OBS the same instant in TAI (+37 s)
        utc_start = 2460934.5 + (10 * 3600 + 1 * 60 + 53.966548) / 86400     # 2025-09-16T10:01:53.966548
        h = FITSHeader(["DATE-OBS", "MJD-OBS", "EXPTIME"],
                       Any["2025-09-16T10:01:53.966548", utc_start - 2400000.5 + 37 / 86400, 45.0],
                       fill("", 3))
        @test frame_epoch(h) ≈ utc_start + 22.5 / 86400 atol=1e-8        # DATE-OBS, mid-exposure
        @test frame_epoch(h; exptime_key=nothing) ≈ utc_start atol=1e-8
        @test frame_epoch(h; timestamp_key="MJD-OBS") ≈ utc_start + (37 + 22.5) / 86400 atol=1e-8
        no_date = FITSHeader(["MJD-OBS"], Any[60000.25], [""])
        @test frame_epoch(no_date) ≈ 2460000.75                           # falls back to MJD-OBS
        @test frame_epoch(FITSHeader(["DATE-OBS"], Any["2020-01-01"], [""]); exptime_key=nothing) ≈ 2458849.5
        @test_throws ArgumentError frame_epoch(FITSHeader(["DATE-OBS"], Any["yesterday"], [""]))
        # other telescopes' conventions, found without a timestamp_key
        ztf = FITSHeader(["OBSJD", "OBSMJD", "EXPTIME"], Any[2458779.6955903, 58779.1955903, 30.0], fill("", 3))
        @test frame_epoch(ztf) ≈ 2458779.6955903 + 15 / 86400 atol=1e-8          # OBSMJD
        @test frame_epoch(FITSHeader(["JD"], Any[2458779.5], [""])) ≈ 2458779.5
        @test frame_epoch(FITSHeader(["DATE-OBS", "TIME-OBS"], Any["2020-01-01", "06:00:00"], fill("", 2))) ≈ 2458849.75
        @test frame_epoch(FITSHeader(["DATE-OBS", "UT"], Any["2020-01-01", "06:00:00.5"], fill("", 2))) ≈
              2458849.75 + 0.5 / 86400 atol=1e-8
        @test_throws ArgumentError frame_epoch(FITSHeader(["EXPTIME"], Any[30.0], [""]))
    end

    @testset "calibrated photometry / inject_movers / injection_recovery" begin
        # three synthetic frames of one star field: 40 catalog stars rendered
        # with a zero point of exactly 25, a known WCS, real-looking headers
        Random.seed!(15)
        nx, ny = 320, 280
        truth_wcs = WCSTransform(2; crpix=[160.0, 140.0], crval=[3.05, -10.74],
                                 cdelt=[-0.257 / 3600, 0.257 / 3600], ctype=["RA---TAN", "DEC--TAN"])
        star_xy = [(15 + 290rand(), 15 + 250rand()) for _ in 1:40]
        gmag = 16.0 .+ 3.0 .* rand(40)
        sky = [pix_to_world(truth_wcs, [x, y]) for (x, y) in star_xy]
        catalog = Table(ra=first.(sky), dec=last.(sky), gmag=gmag)
        zp, σ = 25.0, 1.5
        mktempdir() do dir
            paths = String[]
            for k in 1:3
                raw = 100.0 .+ 5.0 .* randn(nx, ny)                         # FITS-native (x, y)
                for ((x, y), g) in zip(star_xy, gmag), i in max(1, floor(Int, x) - 8):min(nx, ceil(Int, x) + 8),
                    j in max(1, floor(Int, y) - 8):min(ny, ceil(Int, y) + 8)
                    raw[i, j] += 10^(-0.4 * (g - zp)) / (2π * σ^2) * exp(-((i - x)^2 + (j - y)^2) / (2σ^2))
                end
                mjd = 60000.0 + (k - 1) * 0.01
                header = FITSHeader(
                    ["MJD-OBS", "EXPTIME", "GAIN", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"],
                    Any[mjd, 45.0, 1.0, 160.0, 140.0, 3.05, -10.74, -0.257 / 3600, 0.257 / 3600, "RA---TAN", "DEC--TAN"],
                    fill("", 11))
                path = joinpath(dir, "f$k.fits")
                FITS(f -> write(f, raw; header=header), path, "w")
                push!(paths, path)
            end

            image = permutedims(FITS(f -> Float64.(read(f[1])), paths[1]))
            z = photometric_zeropoint(image, truth_wcs, catalog)
            @test z.n_stars >= 20
            @test z.zeropoint ≈ zp atol=0.05      # 5 px aperture holds ~99.6% of a σ=1.5 px star

            injected, truth = inject_movers(paths, joinpath(dir, "inj"); n=4, mag_range=(17.0, 17.5),
                                            speed_range=(20.0, 40.0), rng=Xoshiro(3),
                                            reference_stars=catalog, margin=60)
            @test length(truth) == 4 * 3
            @test all(isfile, injected)
            @test truth.epoch[truth.frame .== 1][1] ≈ 2400000.5 + 60000.0 + 22.5 / 86400 atol=1e-8

            candidates = run_pipeline(injected; threshold=5.0, match_radius=8.0, mask_fill_values=true,
                                      refine_astrometry=true, reference_stars=catalog, difference_stack=true,
                                      min_sharpness=0.3)
            rec = injection_recovery(candidates, truth; radius=1.0)
            # bright and clear of the field edges, but placed at random: one
            # can still land on a star or another synthetic object and be
            # masked out — a real loss, which is what this machinery measures
            @test count(rec.recovered) >= 3
            # and measured at the brightness they were injected with
            for r in filter(r -> r.recovered, collect(rec))
                rows = filter(c -> c.id == r.tracklet, collect(candidates))
                mags = candidate_magnitudes(injected, Table(rows); reference_stars=catalog)
                @test median(filter(isfinite, mags)) ≈ r.mag atol=0.15
            end
            # a run that finds nothing recovers nothing
            @test !any(injection_recovery(candidates[candidates.id .== -1], truth).recovered)
        end
    end

    @testset "MPCORB orbits and :mpcorb crossmatch" begin
        P = AsteroidPipeline
        @test P._unpack_mpc_epoch("K2669") ≈ datetime2julian(DateTime(2026, 6, 9))
        @test P._unpack_mpc_epoch("J9611") ≈ datetime2julian(DateTime(1996, 1, 1))
        @test P._unpack_mpc_epoch("K26A1") ≈ datetime2julian(DateTime(2026, 10, 1))
        @test P._unpack_mpc_epoch("X2669") === nothing

        # a circular, ecliptic, 1 AU orbit: starts on the x axis, a quarter
        # period later sits on the ecliptic y axis (tilted into equatorial)
        circ = (packed="T", name="T", H=10.0, G=0.15, epoch=2451545.0, M=0.0, peri=0.0, node=0.0,
                incl=0.0, e=0.0, a=1.0, n=0.9856076686)
        @test all(isapprox.(P._heliocentric_equatorial(circ, 2451545.0), (1.0, 0.0, 0.0); atol=1e-9))
        ε = deg2rad(P._OBLIQUITY_J2000)
        q = P._heliocentric_equatorial(circ, 2451545.0 + 90 / 0.9856076686)
        @test all(isapprox.(q, (0.0, cos(ε), sin(ε)); atol=1e-9))

        # a real MPCORB.DAT line (2019 UA185, MPC orbit of 2026-10-04) + header + junk
        line = "K19UI5A 19.38  0.15 K2669 342.15513  316.67443   41.65155    2.91030  0.1379969  " *
               "0.27806787   2.3247103  0 E2026-T35    39   7 2005-2026 1.03 M-v 3Ek MPCORBFIT  0000" *
               "          2019 UA185         20261004"
        mktempdir() do dir
            path = joinpath(dir, "MPCORB.DAT")
            write(path, "MINOR PLANET CENTER ORBIT DATABASE (MPCORB)\n" * "-"^160 * "\n" * line * "\nshort line\n")
            orbits = read_mpcorb(path)
            @test length(orbits) == 1
            o = only(orbits)
            @test o.name == "2019 UA185" && o.packed == "K19UI5A"
            @test o.a ≈ 2.3247103 && o.e ≈ 0.1379969 && o.incl ≈ 2.91030 && o.H ≈ 19.38
            @test_throws ArgumentError crossmatch_catalog([(id=1, ra=0.0, dec=0.0, epoch=2461317.9)], :mpcorb; radius=5.0)
            try
                # Astrometrica's measurement on a real IASC frame (XY72_p01, 2026-10-04
                # 10:17:00.6 UTC): older orbits put it degrees away, this one ~4"
                ra, dec = 15 * (1 + 35 / 60 + 44.183 / 3600), 7 + 31 / 60 + 33.29 / 3600
                m = crossmatch_catalog([(id=7, ra=ra, dec=dec, epoch=2461317.928479)], :mpcorb;
                                       radius=15.0, mpcorb=orbits, observatory="F51")
                @test only(m).name == "2019 UA185" && only(m).id == 7
                @test only(m).distance_arcsec < 10
                @test 19 < only(m).mv < 23
            catch e
                AsteroidPipeline._is_unavailable(e) || rethrow()
                @test_skip "network unavailable"
            end
        end
    end

    @testset "detector-fixed and bright-background filters; set_quality" begin
        P = AsteroidPipeline
        det(xs, ys) = Table(x=round.(Int, xs), y=round.(Int, ys), peak=fill(10.0, length(xs)), flux=fill(100.0, length(xs)),
                            flux_err=fill(5.0, length(xs)), xcen=Float64.(xs), ycen=Float64.(ys), sharpness=fill(0.6, length(xs)))
        # a hot pixel at (40, 40) in every frame's own grid; a mover elsewhere
        frames = [det([40.0, 100.0 + 10k], [40.0, 60.0]) for k in 1:4]
        kept = P._drop_detector_fixed(frames, 2.0)
        @test all(f -> !any(d -> d.xcen == 40.0, f), kept)
        @test all(length.(kept) .== 1)
        # only twice: could be chance, kept
        twice = [det([40.0], [40.0]), det([40.0], [40.0]), det([70.0], [70.0]), det([90.0], [90.0])]
        @test length.(P._drop_detector_fixed(twice, 2.0)) == [1, 1, 1, 1]

        Random.seed!(16)
        img = 100.0 .+ 5.0 .* randn(120, 120)
        for i in 1:120, j in 1:120                       # a bright star's halo around (90, 90)
            img[i, j] += 400.0 * exp(-hypot(i - 90, j - 90) / 15)
        end
        d = det([30.0, 85.0], [30.0, 85.0])               # clean sky; inside the halo
        @test [x for x in P._drop_on_bright_background(d, img, nothing, 1.0).xcen] == [30.0]
        @test length(P._drop_on_bright_background(d, img, nothing, Inf)) == 2

        @test P._check_frame_quality("ok.fits", 900, 170) < 30
        @test (@test_logs (:warn,) P._check_frame_quality("bad.fits", 160_000, 170)) > 30
    end

    @testset "the MPC orbit file only near its epoch" begin
        orbits = [(epoch=2461200.5, name="x")]                          # 2026-06-09
        w = String[]
        @test AsteroidPipeline._mpcorb_for(orbits, [2461330.5], w) === orbits   # 130 days: used
        @test isempty(w)
        @test AsteroidPipeline._mpcorb_for(orbits, [2458779.7], w) === nothing  # a 2019 image: not
        @test occursin("SkyBoT", only(w))
        @test AsteroidPipeline._mpcorb_for(nothing, [2458779.7], w) === nothing
    end

    @testset "analyze_set on a synthetic set" begin
        # Four frames of one field, 15 min apart: 60 stars, a header 3″ off,
        # and one asteroid at 30″/h. The reference stars are given, so only
        # the known-object services (SkyBoT) use the network — and an outage
        # of theirs is a warning.
        Random.seed!(29)
        nx, ny, scale = 500, 500, 0.26
        truth = WCSTransform(2; crpix=[250.0, 250.0], crval=[210.0, 45.0],
                             cdelt=[-scale / 3600, scale / 3600], ctype=["RA---TAN", "DEC--TAN"])
        star_xy = [(20 + 460rand(), 20 + 460rand()) for _ in 1:60]
        gmag = 15 .+ 4 .* rand(60)
        sky = [pix_to_world(truth, [x, y]) for (x, y) in star_xy]
        stars = Table(ra=first.(sky), dec=last.(sky), gmag=gmag)
        zp, sigma = 27.0, 1.6
        add!(img, x, y, flux) = for i in max(1, floor(Int, x) - 8):min(nx, ceil(Int, x) + 8),
                                    j in max(1, floor(Int, y) - 8):min(ny, ceil(Int, y) + 8)
            img[i, j] += flux * exp(-((i - x)^2 + (j - y)^2) / (2sigma^2)) / (2pi * sigma^2)
        end
        mover(k) = (150.0 + 29.0 * (k - 1), 200.0 + 6.0 * (k - 1))          # px, ~30″/h
        mover_mag = 18.0
        dir = mktempdir()
        for k in 1:4
            img = 200.0 .+ 10.0 .* randn(nx, ny)
            for ((x, y), g) in zip(star_xy, gmag)
                add!(img, x, y, 10^((zp - g) / 2.5))
            end
            add!(img, mover(k)..., 10^((zp - mover_mag) / 2.5))
            header = FITSHeader(
                ["CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                 "DATE-OBS", "EXPTIME"],
                Any["RA---TAN", "DEC--TAN", 250.0 + 12, 250.0 - 8, 210.0, 45.0, -scale / 3600, scale / 3600,
                    "2020-03-01T08:$(lpad(15 * (k - 1), 2, '0')):00", 45.0],
                fill("", 10))
            FITS(f -> write(f, img; header), joinpath(dir, "frame$k.fits"), "w")
        end

        r = analyze_set(dir; workdir=mktempdir(), observatory="F51", reference_stars=stars)
        @test r.reliable && r.frames_used == 1:4
        @test r.epochs[1] ≈ datetime2julian(DateTime(2020, 3, 1, 8)) + 22.5 / 86400 atol=1e-7
        @test length(unique(r.candidates.id)) == 1                # the asteroid, and nothing else
        @test r.report_kind === :candidates && only(r.unknown_ids) == r.candidates.id[1]
        for c in r.candidates
            ra, dec = pix_to_world(truth, collect(mover(c.frame)))
            @test hypot((c.ra - ra) * cosd(dec), c.dec - dec) * 3600 < 0.3   # through the refitted WCS
            @test c.mag ≈ mover_mag atol=0.3
        end
        @test count(l -> length(l) == 80 && endswith(l, "F51"), split(r.report, "\n")) == length(r.candidates)   # one 80-column line each
        @test isempty(r.dropped)
    end

    @testset "moving_light_curve" begin
        # 24 frames over 6 h; an asteroid at 10″/h whose brightness swings ±0.3 mag
        # with a 0.15 d period, given only by its positions in the first 4 frames.
        # In frame 10 it is missing: measured where predicted, not recentred.
        Random.seed!(31)
        nx, ny, scale = 500, 500, 0.26
        truth = WCSTransform(2; crpix=[250.0, 250.0], crval=[150.0, 10.0],
                             cdelt=[-scale / 3600, scale / 3600], ctype=["RA---TAN", "DEC--TAN"])
        mover(h) = (130.0 + 38.0 * h, 150.0 + 12.0 * h)                    # px; ~10″/h
        # no star within 20 px of the track: crossing one is real, but not what this measures
        clear(x, y) = minimum(hypot(x - mover(h)[1], y - mover(h)[2]) for h in 0:0.05:6) > 20
        star_xy = Tuple{Float64,Float64}[]
        while length(star_xy) < 60
            xy = (20 + 460rand(), 20 + 460rand())
            clear(xy...) && push!(star_xy, xy)
        end
        gmag = 15 .+ 4 .* rand(60)
        push!(star_xy, mover(5.75) .+ (2.0, 1.0)); push!(gmag, 17.5)    # one on the track, in frame 24
        sky = [pix_to_world(truth, [x, y]) for (x, y) in star_xy]
        stars = Table(ra=first.(sky), dec=last.(sky), gmag=gmag)
        zp, sigma = 27.0, 1.6
        add!(img, x, y, flux) = for i in max(1, floor(Int, x) - 8):min(nx, ceil(Int, x) + 8),
                                    j in max(1, floor(Int, y) - 8):min(ny, ceil(Int, y) + 8)
            img[i, j] += flux * exp(-((i - x)^2 + (j - y)^2) / (2sigma^2)) / (2pi * sigma^2)
        end
        hours = [0.25 * (k - 1) for k in 1:24]
        period = 0.15
        mover_mag(h) = 17.0 + 0.3 * sin(2pi * h / 24 / period)
        start = datetime2julian(DateTime(2020, 3, 1, 4))
        dir = mktempdir()
        paths = map(enumerate(hours)) do (k, h)
            img = 200.0 .+ 10.0 .* randn(nx, ny)
            for ((x, y), g) in zip(star_xy, gmag)
                add!(img, x, y, 10^((zp - g) / 2.5))
            end
            k == 10 || add!(img, mover(h)..., 10^((zp - mover_mag(h)) / 2.5))
            minutes = round(Int, 60h)
            header = FITSHeader(
                ["CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "DATE-OBS"],
                Any["RA---TAN", "DEC--TAN", 250.0 + 4, 250.0 - 3, 150.0, 10.0, -scale / 3600, scale / 3600,
                    "2020-03-01T$(lpad(4 + minutes ÷ 60, 2, '0')):$(lpad(minutes % 60, 2, '0')):00"],
                fill("", 9))
            path = joinpath(dir, "frame$k.fits")
            FITS(f -> write(f, img; header), path, "w")
            path
        end
        tracklet = [(ra=p[1], dec=p[2], epoch=start + h / 24)
                    for (h, p) in ((h, pix_to_world(truth, collect(mover(h)))) for h in hours[1:4])]

        lc = moving_light_curve(paths, tracklet; reference_stars=stars)
        @test length(lc) == 24 && lc.epoch ≈ start .+ hours ./ 24
        @test findall(lc.star_nearby) == [24]
        seen = [k != 10 for k in 1:24]
        @test all(lc.recentred[seen]) && !lc.recentred[10]
        seen[24] = false
        for r in lc[seen]
            x, y = mover(hours[r.frame])
            @test hypot(r.x - x, r.y - y) < 0.3                            # extrapolated 6 h, then recentred
            @test r.mag ≈ mover_mag(hours[r.frame]) atol = 0.05
        end
        @test abs(lc.flux[10]) < 3 * lc.flux_err[10]                      # nothing there, measured anyway
        @test hypot(lc.x[10] - mover(hours[10])[1], lc.y[10] - mover(hours[10])[2]) < 0.3
        p = recover_rotation_period(lc.epoch[seen], lc.mag[seen]; minimum_period=0.05, maximum_period=0.5)
        @test p.period ≈ period rtol = 0.05
        @test_throws ArgumentError moving_light_curve(paths, NamedTuple[])
        @test all(isfinite, lc.zeropoint_err)
        @test all(lc.mag_err[seen] .>= lc.zeropoint_err[seen])

        # comparison stars: each with its own colour offset, frames with their own
        # transparency, one star saturated (scattering) and one missing in a frame
        Random.seed!(5)
        transparency = [0.0, 0.05, -0.03, 0.12, 0.02]
        colour = 0.3 .* randn(8)
        v = [[25.0 + transparency[k] + colour[i] + 0.002randn() for i in 1:8] for k in 1:5]
        for k in 1:5
            v[k][7] += 0.2randn()
        end
        v[3][8] = NaN
        zp, err = AsteroidPipeline._comparison_zeropoints(v)
        @test zp .- zp[1] ≈ transparency .- transparency[1] atol = 0.003
        @test all(e -> 0 < e < 0.003, err)
        @test AsteroidPipeline._comparison_zeropoints([[25.0, NaN, 24.0, NaN] for _ in 1:3]) == (nothing, nothing)
    end

    @testset "reduced_magnitudes / object geometry" begin
        # H = 10 seen at three geometries, through the H-G model the MPC uses
        H, G = 10.0, 0.15
        geometry = [(r=2.0, delta=1.0, phase=0.0), (r=2.5, delta=1.6, phase=12.0), (r=3.0, delta=2.2, phase=19.0)]
        mags = [H + 5log10(g.r * g.delta) + AsteroidPipeline._hg_phase_dimming(deg2rad(g.phase), G) for g in geometry]
        epochs = [2460000.5, 2460010.5, 2460020.5]
        red = reduced_magnitudes(epochs, mags, geometry; G)
        @test red.mag ≈ fill(H, 3)
        @test (epochs[1] - red.epoch[1]) * 86400 ≈ 499.0 atol = 0.1       # light time for 1 AU
        @test reduced_magnitudes(epochs, mags, geometry; G=nothing).mag ≈ mags .- 5log10.([2.0, 4.0, 6.6])
        @test AsteroidPipeline._hg_phase_dimming(0.0, G) ≈ 0 atol = 1e-12
        @test_throws ArgumentError reduced_magnitudes(epochs[1:2], mags, geometry)

        response = """
            *******************************************************************************
             Date__(UT)__HR:MN:SC.fff, , ,               r,       rdot,             delta,     deldot,     S-T-O,
            *******************************************************************************
            \$\$SOE
             2022-Aug-09 00:00:00.000, , ,  5.268747694354, -0.2950553,  5.76723809407072,-21.2813335,    9.1862,
             2022-Sep-28 04:48:00.000,*, ,  5.260182936355, -0.2954448,  5.05448490840793,-25.7287738,   10.9186,
            \$\$EOE
            """
        rows = AsteroidPipeline._horizons_geometry(response)
        @test length(rows) == 2
        @test rows[2] == (r=5.260182936355, delta=5.05448490840793, phase=10.9186)
        @test_throws ErrorException AsteroidPipeline._horizons_geometry("No matches found.")
        @test length(object_geometry(624, Float64[])) == 0
    end

    @testset "image sets, MPC reports and their checks" begin
        P = AsteroidPipeline
        # packed designations, as Astrometrica wrote them in real reports
        @test pack_designation("2016 SB115") == "K16SB5B"
        @test pack_designation("2019 UA185") == "K19UI5A"
        @test pack_designation("(656932) 2016 FA12") == "~09bg"
        @test pack_designation("(550107) 2012 AZ30") == "t0107"
        @test pack_designation("433") == "00433"
        @test pack_designation("2016 SX110") == "K16SB0X"
        @test pack_designation("2008 EK13") == "K08E13K"
        @test pack_designation("4311 T-1") == "T1S4311"
        @test pack_designation("Hasslacher") === nothing

        cfg = IASC_PS1_CONFIG
        @test cfg.match_radius ≈ 2.0 / 0.2563
        @test pipeline_config(1.0).match_radius ≈ 2.0 && pipeline_config(1.0).min_speed ≈ 72.0

        header = mpc_report_header("F51"; observers="J. Bulger", measurers="A. Name", telescope="1.8-m",
                                   net="Gaia DR3", ack_time=DateTime(2026, 10, 18, 12, 57, 58))
        @test header == ["COD F51", "OBS J. Bulger", "MEA A. Name", "TEL 1.8-m",
                         "ACK MPCReport file updated 2026.10.18 12:57:58", "NET Gaia DR3"]
        @test_throws ArgumentError mpc_report_header("F5")
        blank = blank_mpc_report("F51")
        @test startswith(blank, "COD F51\nACK ") && endswith(blank, "----- end -----\n")
        rows = [(id=2, frame=k, ra=183.0 + k / 1000, dec=-10.7, epoch=2460934.9 + k / 100, mag=21.4) for k in 1:3]
        rep = mpc_report(rows, "F51"; trksub_prefix="P")
        obs = filter(l -> startswith(l, "     P2"), split(rep, '\n'))
        @test length(obs) == 3 && all(length.(obs) .== 80)
        @test mpc_report(rows[1:0], "F51") |> r -> !occursin("P2", r)

        # the reports of a real campaign day (2026-10-18), and what was wrong with them
        line_sb115 = "     K16SB5B  C2026 10 05.53420201 27 04.617-04 38 46.72         19.6 G      F51"
        line_fa12 = "~09bg         C2026 10 05.54588701 27 40.068-04 54 40.86         19.0 G      F51"
        p = P._parse_mpc80_line(line_fa12)
        @test p.number == "~09bg" && p.station == "F51"
        @test p.ra ≈ 15 * (1 + 27 / 60 + 40.068 / 3600) && p.dec ≈ -(4 + 54 / 60 + 40.86 / 3600)
        @test p.jd ≈ datetime2julian(DateTime(2026, 10, 5)) + 0.545887
        good = join([header; "----- end -----"], "\n")
        c = check_mpc_report(good; observatory="F51")
        @test c.ok && isempty(c.issues)
        c = check_mpc_report("ACK MPCReport file updated 2026.10.18 12:33:05\n----- end -----\n"; observatory="F51")
        @test !c.ok && [i.code for i in c.issues] == [:missing_cod]
        @test startswith(c.cleaned, "COD F51\nACK")
        c = check_mpc_report(replace(good, "COD F51" => "COD XXX"); observatory="F51")
        @test !c.ok && only(c.issues).code == :wrong_code && only(c.issues).detail == "XXX"
        with_obs = join([header; line_sb115; line_fa12; "----- end -----"], "\n")
        c = check_mpc_report(with_obs; observatory="F51", known=[(name="(656932) 2016 FA12", packed="~09bg")])
        @test c.ok
        @test [(i.code, i.line) for i in c.issues] == [(:single_observation, 7), (:known_object, 8)]
        @test occursin("K16SB5B", c.cleaned) && !occursin("~09bg", c.cleaned)
        # same exposures, another chip: the time can't tell, the footprint can
        c = check_mpc_report(with_obs; observatory="F51", epochs=[2461319.0342, 2461319.5459],
                             field=(ra, dec) -> ra > 21.9)
        @test !c.ok && count(i -> i.code == :other_set, c.issues) == 2
        @test check_mpc_report(""; observatory="F51").issues[1].code == :empty
        @test any(i -> i.code == :missing_end, check_mpc_report("COD F51\n"; observatory="F51").issues)
        @test any(i -> i.code == :malformed_line, check_mpc_report("COD F51\nhello\n"; observatory="F51").issues)

        mktempdir() do dir
            settings = joinpath(dir, "PS1.cfg")
            write(settings, "[Settings]\r\nCatalog=10\r\n[ObservingSite]\r\nMPCCode=F51\r\nObserver=J. Bulger, T. Lowe\r\n" *
                            "Measurer=T. Vorobjov\r\nTelescope=1.8-m f/4.4 Ritchey-Chretien + CCD\r\n")
            s = astrometrica_settings(settings)
            @test s.code == "F51" && s.observers == "J. Bulger, T. Lowe" && s.net == "Gaia DR1"
            @test s.telescope == "1.8-m f/4.4 Ritchey-Chretien + CCD"

            # a .zip of FITS frames with a folder inside, BLANK pixels, and a stray file
            src = mkpath(joinpath(dir, "src", "XY99_p00"))
            for k in 1:3
                FITS(joinpath(src, "frame$k.fits"), "w") do f
                    img = fill(Int16(100 + k), 20, 10); img[1, 1] = Int16(-32768)
                    write(f, img)
                    write_key(f[1], "BLANK", -32768)   # FITSIO drops it from a header passed to write
                end
            end
            write(joinpath(src, "notes.txt"), "not a frame")
            zip = joinpath(dir, "ps1-20990101_9_XY99_p00.zip")
            run(pipeline(`$(P.p7zip_jll.p7zip()) a -tzip $zip $src`; stdout=devnull))
            set = extract_set(zip, joinpath(dir, "out"))
            @test set.name == "XY99_p00" && basename.(set.paths) == ["frame1.fits", "frame2.fits", "frame3.fits"]
            @test isempty(filter(n -> startswith(n, ".extract-"), readdir(joinpath(dir, "out"))))
            @test extract_set(zip, joinpath(dir, "out")).paths == set.paths       # reused, not re-extracted
            cleaned = clean_blank_pixels(set.paths[1]; outdir=joinpath(dir, "clean"))
            @test cleaned != set.paths[1] && FITS(f -> read(f[1])[1, 1], cleaned) ≈ 101

            # an old orbit file is kept when the MPC can't be reached
            old = joinpath(dir, "MPCORB.DAT"); write(old, "old copy\n")
            touch(old); mtime_before = mtime(old)
            @test (@test_logs (:warn,) match_mode=:any download_mpcorb(old; max_age_hours=0,
                                                                      url="http://127.0.0.1:9/nothing.gz")) == old
            @test read(old, String) == "old copy\n" && !isfile(old * ".part")
            @test download_mpcorb(old; max_age_hours=24) == old
        end

        w = WCSTransform(2; ctype=["RA---TAN", "DEC--TAN"], crval=[30.0, 0.0], crpix=[50.0, 50.0],
                         cdelt=[-1 / 3600, 1 / 3600])
        on = frame_footprint([w], (100, 100))
        @test on(30.0, 0.0) && !on(30.1, 0.0)
        @test P._pixel_scale_arcsec(w, (100, 100)) ≈ 1.0 atol = 1e-6
    end

    @testset "observatory from the frames' headers" begin
        mktempdir() do dir
            # lines of the MPC's ObsCodes.html, as published (2026-10)
            path = joinpath(dir, "ObsCodes.html")
            write(path, "<pre>\nCode  Long.   cos      sin    Name\n" *
                        "250                           Hubble Space Telescope\n" *
                        "644 243.140220.836325+0.546877Palomar Mountain/NEAT\n" *
                        "675 243.137460.836357+0.546831Palomar Mountain\n" *
                        "F51 203.744090.936241+0.351543Pan-STARRS 1, Haleakala\n" *
                        "F52 203.744090.936239+0.351545Pan-STARRS 2, Haleakala\n" *
                        "I41 243.140220.836322+0.546875Palomar Mountain--ZTF\n" *
                        "T05 203.744440.936241+0.351543ATLAS-HKO, Haleakala\n</pre>\n")
            codes = read_obscodes(path)
            @test [c.code for c in codes] == ["644", "675", "F51", "F52", "I41", "T05"]   # Hubble: no site
            @test codes[3].lon ≈ 203.74409 && codes[3].rho_cos ≈ 0.936241 && codes[3].name == "Pan-STARRS 1, Haleakala"

            # Pan-STARRS's own keywords: west longitude in hours
            ps1 = FITSHeader(["TELESCOP", "INSTRUME", "ORIGIN", "FPA.LONGITUDE", "FPA.LATITUDE", "FPA.ELEVATION"],
                             Any["PS1", "gpc1", "PS1", 10.4170602777931, 20.7070994446013, 3048.0], fill("", 6))
            s = frame_site(ps1)
            @test s.telescope == "PS1" && isapprox(s.lon, -156.2559; atol=1e-4) && isapprox(s.lat, 20.7071; atol=1e-4)
            best = observatory_candidates(s, codes)
            @test best[1].code == "F51" && best[1].distance_km < 0.1    # F52 is 10 m away; the name decides
            @test Set(c.code for c in best) == Set(["F51", "F52", "T05"])
            # ZTF's: east longitude in degrees, names spelled out
            ztf = FITSHeader(["TELESCOP", "INSTRUME", "ORIGIN", "OBSERVAT", "OBSLON", "OBSLAT", "OBSALT"],
                             Any["Palomar 48-inch", "ZTF/MOSAIC", "Zwicky Transient Facility", "Palomar Observatory",
                                 -116.8598, 33.3573, 1668.0], fill("", 7))
            @test observatory_candidates(frame_site(ztf), codes)[1].code == "I41"
            # sexagesimal strings, and no site at all: names alone
            @test frame_site(FITSHeader(["SITELONG", "SITELAT"], Any["-116 51 35.3", "+33:21:26"], ["", ""])).lon ≈
                  -(116 + 51 / 60 + 35.3 / 3600)
            nosite = frame_site(FITSHeader(["TELESCOP"], Any["PS1"], [""]))
            @test nosite.lon === nothing && observatory_candidates(nosite, codes)[1].code == "F51"
            @test isempty(observatory_candidates(frame_site(FITSHeader(["X"], Any[1], [""])), codes))
        end
    end

    @testset "motion of a tracklet, and known objects near a candidate" begin
        P = AsteroidPipeline
        # 36"/h due west (PA 270°), on the equator, every 15 minutes
        line = [(id=1, ra=10.0 - k * 9 / 3600, dec=0.0, epoch=2461000.5 + k / 96) for k in 0:3]
        m = tracklet_motion(line)
        @test m.rate_arcsec_h ≈ 36 rtol = 1e-6
        @test m.pa_deg ≈ 270 atol = 1e-6
        @test m.rms_arcsec < 1e-6
        bent = [merge(p, (dec=p.dec + (k == 2) * 2 / 3600,)) for (k, p) in enumerate(line)]
        @test tracklet_motion(bent).rms_arcsec > 0.5
        @test tracklet_motion(line[1:1]).rate_arcsec_h == 0

        # the same object displaced along its track (an orbit's timing error)
        # agrees; across it, or moving otherwise, doesn't
        pred = (rate_arcsec_h=36.0, pa_deg=270.0)
        @test P._same_motion(m, pred, (-120.0, 0.0))          # 2' ahead, westward
        @test !P._same_motion(m, pred, (0.0, 60.0))           # 1' north: across
        @test !P._same_motion(m, (rate_arcsec_h=30.0, pa_deg=270.0), (-120.0, 0.0))
        @test !P._same_motion(m, (rate_arcsec_h=36.0, pa_deg=262.0), (-120.0, 0.0))

        # with no catalogue to ask, every candidate is new
        c = nearby_known_objects(line; observatory="F51", skybot=false)
        @test c.checks[1].verdict === :new && isempty(c.checks[1].objects)
        @test isempty(nearby_known_objects(line[1:0]; observatory="F51", skybot=false).checks)

        # live: 2019 UA185 on a real 2026-10-04 IASC frame (see the MPCORB testset),
        # measured where its MPC orbit puts it, and 2' further along its track
        orbit = (packed="K19UI5A", name="2019 UA185", H=19.38, G=0.15, epoch=P._unpack_mpc_epoch("K2669"),
                 M=342.15513, peri=316.67443, node=41.65155, incl=2.91030, e=0.1379969, a=2.3247103, n=0.27806787)
        try
            epochs = [2461317.928479 + k * 0.0122 for k in 0:3]
            obs = P._observer_heliocentric(epochs, "F51")
            at = [P._mpcorb_radec(orbit, e, o) for (e, o) in zip(epochs, obs)]
            measured = [(id=7, ra=q.ra, dec=q.dec, epoch=e) for (q, e) in zip(at, epochs)]
            c = nearby_known_objects(measured; observatory="F51", mpcorb=[orbit], skybot=false)
            @test c.checks[7].verdict === :known && only(c.checks[7].objects).agreement === :position
            mo = tracklet_motion(measured)
            ahead = [merge(p, (ra=p.ra + 120 * sind(mo.pa_deg) / 3600 / cosd(p.dec),
                               dec=p.dec + 120 * cosd(mo.pa_deg) / 3600)) for p in measured]
            c = nearby_known_objects(ahead; observatory="F51", mpcorb=[orbit], skybot=false)
            @test c.checks[7].verdict === :possibly_known
            @test only(c.checks[7].objects).separation_arcsec ≈ 120 atol = 5
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end
    end

    @testset "own designations; the NEO Confirmation Page" begin
        P = AsteroidPipeline
        rows = [(id=2, frame=k, ra=183.0 + k / 1000, dec=-10.7, epoch=2460934.9 + k / 100, mag=21.4) for k in 1:2]
        rep = mpc_report(rows, "F51"; designations=Dict(2 => "NHU0002"))
        @test count(l -> startswith(l, "     NHU0002  C"), split(rep, '\n')) == 2
        @test_throws ArgumentError mpc80_report(rows, "F51"; designations=Dict(2 => "NHU 2"))
        @test_throws ArgumentError mpc80_report(rows, "F51"; designations=Dict(2 => "NHU00002"))
        @test startswith(split(ades_psv(rows, "F51"; designations=Dict(2 => "ZTF10Go")), '\n')[2], "ZTF10Go|")

        # the MPC's ephemeris page, as it answered (2026-10-18), abridged
        html = """<pre>
         ZTF10Go
        <a href="x">Get the observations</a>
        Date       UT   *  R.A. (J2000) Decl.  Elong.  V        Motion     Object
        2026 10 18 2200   01 42 01.5 -08 46 22 161.4  17.9   45.22  078.4  160  -77   +60    0.52  082  -16   Map/Offsets   !
        2026 10 18 2230   01 43 29.8 -08 41 54 161.4  17.9   43.82  078.5  193  -78   +59    0.53  082  -10   Map/Offsets   !
         ABC0069
        2026 10 18 2200   05 22 59.7 +37 32 48 119.7  19.0   54.86  094.7  135  -08   +60    0.34  165  +03   Map/Offsets
        </pre>"""
        eph = P._parse_neocp_ephemerides(html, ["ZTF10Go", "ABC0069"])
        @test length(eph["ZTF10Go"]) == 2 && length(eph["ABC0069"]) == 1
        r = eph["ZTF10Go"][1]
        @test r.jd ≈ datetime2julian(DateTime(2026, 10, 18, 22, 0))
        @test r.ra ≈ 15 * (1 + 42 / 60 + 1.5 / 3600) && r.dec ≈ -(8 + 46 / 60 + 22 / 3600)
        @test r.rate_arcsec_h ≈ 45.22 * 60 && r.pa_deg ≈ 78.4 && r.uncertain
        @test !eph["ABC0069"][1].uncertain
        mid = P._interpolate_ephemeris(eph["ZTF10Go"], r.jd + 15 / 1440)
        @test mid.ra ≈ (eph["ZTF10Go"][1].ra + eph["ZTF10Go"][2].ra) / 2
        @test P._interpolate_ephemeris(eph["ZTF10Go"], r.jd - 1) === nothing

        # old images aren't checked: the NEOCP lists what is being observed now
        old = [(id=1, ra=10.0, dec=0.0, epoch=2460000.5 + k / 96) for k in 0:3]
        c = neocp_check(old; observatory="F51", now_jd=2460100.5)
        @test isempty(c.checks) && occursin("100 days old", only(c.warnings))

        # live: a tracklet built from an object now on the NEOCP is recognised
        try
            objs = neocp_objects()
            if !isempty(objs)
                o = objs[1]
                t0 = datetime2julian(P.Dates.now(P.Dates.UTC)) - 3 / 24
                e = P._neocp_ephemerides([o.designation], "F51", t0)[o.designation]
                track = [(id=1, P._interpolate_ephemeris(e, t0 + k / 96)...) for k in 0:3]
                c = neocp_check(track; observatory="F51")
                @test c.checks[1].verdict === :on_neocp
                @test first(c.checks[1].objects).designation == o.designation
            end
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end
    end

    @testset "crossmatch_catalog" begin
        @test_throws ArgumentError crossmatch_catalog([], :unknown_catalog; radius=5.0)

        try
            vega = [(id=1, ra=279.23473479, dec=38.78368896)]
            matches = crossmatch_catalog(vega, :simbad; radius=5.0)
            @test length(matches) >= 1
            @test matches[1].id == 1
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end

        try
            # a real VSX variable (Gaia DR3 2501737571391422464, type RS) in the
            # real_data_demo.jl field: a positive control
            star = [(id=1, ra=36.2344, dec=2.06997)]
            matches = crossmatch_catalog(star, :vsx; radius=5.0)
            @test any(==("RS"), matches.class)
            @test matches[1].id == 1
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end

        try
            # batched cross-matching: 2 real VSX variables and 1 non-match in one
            # request, each row attributed to the right candidate id
            batch = [(id=10, ra=36.2344, dec=2.06997),
                     (id=20, ra=83.1937, dec=5.41603),
                     (id=30, ra=0.0, dec=89.9)]
            matches = crossmatch_catalog(batch, :vsx; radius=5.0)
            @test length(matches) == 2
            @test Set(matches.id) == Set([10, 20])
            @test only(matches[matches.id.==10].class) == "RS"
            @test only(matches[matches.id.==20].class) == "EW"
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end

        try
            # near the pole, tiny radius: an empty TAP result must not crash
            empty_candidates = [(id=1, ra=0.0, dec=89.9)]
            matches = crossmatch_catalog(empty_candidates, :simbad; radius=1.0)
            @test length(matches) == 0
            @test isempty(matches.id)
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end

        try
            # a known object (127319 "2002 JB99") ~134" away: a positive control,
            # since an all-empty result can't tell "correctly empty" from broken
            known = [(id=1, ra=36.521300910000, dec=2.127107410000, epoch=2459288.617372700)]
            matches = crossmatch_catalog(known, :skybot; radius=200.0)
            @test any(==("2002 JB99"), matches.name)
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end

        try
            # near the pole, nothing inside: an empty result must not crash
            empty_candidates = [(id=1, ra=0.0, dec=89.9, epoch=2459288.6174)]
            matches = crossmatch_catalog(empty_candidates, :skybot; radius=1.0)
            @test length(matches) == 0
            @test isempty(matches.id)
        catch e
            AsteroidPipeline._is_unavailable(e) || rethrow()
            @test_skip "network unavailable"
        end
    end

    @testset "plate_solve" begin
        # pure request/response helpers, no network needed
        body = AsteroidPipeline._astrometry_login_body("mykey123")
        @test body == "request-json=" * HTTP.escapeuri("""{"apikey":"mykey123"}""")

        success = AsteroidPipeline._parse_astrometry_response(
            """{"status":"success","session":"abc123"}""", "login")
        @test success["session"] == "abc123"

        @test_throws ErrorException AsteroidPipeline._parse_astrometry_response(
            """{"status":"error","errormessage":"bad apikey"}""", "login")
        try
            AsteroidPipeline._parse_astrometry_response(
                """{"status":"error","errormessage":"bad apikey"}""", "login")
        catch e
            @test occursin("bad apikey", sprint(showerror, e))
        end

        # _astrometry_poll: succeeds once `check` returns something, times out if
        # never (the other astrometry.net calls need the network)
        calls = Ref(0)
        result = AsteroidPipeline._astrometry_poll(0.3, 0.05, "test condition") do
            calls[] += 1
            calls[] >= 3 ? :done : nothing
        end
        @test result == :done
        @test calls[] == 3
        @test_throws ErrorException AsteroidPipeline._astrometry_poll(
            () -> nothing, 0.1, 0.05, "test condition")

        # live round trip, only with an API key. It needs a real star field (one
        # synthetic source can't be solved), so it uses a real ZTF frame.
        api_key = get(ENV, "ASTROMETRY_API_KEY", nothing)
        real_frame = joinpath(@__DIR__, "..", "data", "real", "science",
                               "ztf_20191023195590_000451_zr_c01_o_q1_sciimg.fits")
        if api_key === nothing
            @test_skip "ASTROMETRY_API_KEY not set"
        elseif !isfile(real_frame)
            @test_skip "real ZTF test frame not present — run examples/fetch_data.sh"
        else
            try
                solved = plate_solve(real_frame; api_key=api_key)
                @test solved isa WCSTransform
            catch e
                # the service being down skips; a refused key still fails
                AsteroidPipeline._is_service_failure(e) || rethrow()
                @test_skip "astrometry.net unavailable"
            end
        end
    end

    @testset "julian_date_to_iso8601" begin
        @test AsteroidPipeline.julian_date_to_iso8601(2451545.0) == "2000-01-01T12:00:00.000Z"  # J2000.0
        @test AsteroidPipeline.julian_date_to_iso8601(2440587.5) == "1970-01-01T00:00:00.000Z"  # Unix epoch

        # 2451545.0 + 1 s must round up to the whole second, not truncate to :59.999
        one_second = 1.0 / 86400
        @test AsteroidPipeline.julian_date_to_iso8601(2451545.0 + one_second) == "2000-01-01T12:00:01.000Z"

        # rounding up to midnight on the last day of a month moves to the next month
        just_before = 2461345.5 - 1e-4 / 86400                # 2026-10-31 23:59:59.9999
        @test AsteroidPipeline.julian_date_to_iso8601(just_before) == "2026-11-01T00:00:00.000Z"
        @test AsteroidPipeline._mpc_date(just_before) == "2026 11 01.000000"
    end

    @testset "ades_psv" begin
        candidates = Table(id=[1, 1, 2], frame=[1, 2, 1],
                            x=[10.0, 12.0, 20.0], y=[10.0, 12.0, 20.0],
                            ra=[150.123456789, 150.123556789, 200.5],
                            dec=[20.987654321, 20.987754321, -10.25],
                            epoch=[2451545.0, 2451545.01, 2451545.0])

        @test_throws ArgumentError ades_psv(candidates, "XX")  # not 3 characters

        psv = ades_psv(candidates, "I41")
        lines = split(strip(psv), "\n")
        @test lines[1] == "trkSub|mode|stn|obsTime|ra|dec"
        @test length(lines) == 4  # header + 3 observations

        row1 = split(lines[2], "|")
        @test row1[1] == uppercase(string(1; base=36))  # trkSub from id=1
        @test row1[2] == "CCD"
        @test row1[3] == "I41"
        @test row1[4] == "2000-01-01T12:00:00.000Z"
        @test parse(Float64, row1[5]) ≈ 150.123456789 atol=1e-6
        @test parse(Float64, row1[6]) ≈ 20.987654321 atol=1e-6

        # both rows sharing id=1 must share the same trkSub — that's how
        # MPC correlates them into one tracklet on their end
        row2 = split(lines[3], "|")
        @test row2[1] == row1[1]
        row3 = split(lines[4], "|")
        @test row3[1] != row1[1]  # id=2 gets a distinct trkSub

        # optional columns only appear when actually supplied
        psv_with_cat = ades_psv(candidates, "I41"; astCat="Gaia2", band="G")
        header_with_cat = split(split(strip(psv_with_cat), "\n")[1], "|")
        @test "astCat" in header_with_cat
        @test "band" in header_with_cat
        @test "photCat" ∉ header_with_cat

        # trkSub length limit: a prefix long enough to push a real id
        # over 8 characters must error rather than silently truncate
        @test_throws ArgumentError ades_psv(candidates, "I41"; trksub_prefix="TOOLONGPREFIX")
    end

    @testset "mpc80_report" begin
        candidates = Table(id=[1, 1, 2], frame=[1, 2, 1],
                            x=[10.0, 12.0, 20.0], y=[10.0, 12.0, 20.0],
                            ra=[150.123456789, 150.123556789, 200.5],
                            dec=[20.987654321, 20.987754321, -10.25],
                            epoch=[2451545.0, 2451545.01, 2451545.0])

        @test_throws ArgumentError mpc80_report(candidates, "XX")  # not 3 characters

        report = mpc80_report(candidates, "I41")
        # chomp, not strip: line 1's 5 leading blank columns are content
        lines = split(chomp(report), "\n")
        @test length(lines) == 3
        @test all(l -> length(l) == 80, lines)  # the whole point of the format

        line1 = lines[1]
        @test line1[1:5] == "     "                                # no permanent number
        @test strip(line1[6:12]) == uppercase(string(1; base=36))  # temp designation from id=1
        @test line1[15] == 'C'                                     # note2 default (CCD)
        @test line1[78:80] == "I41"
        @test line1[66:71] == " "^6                               # no mag column -> blank

        # with magnitudes: one decimal in 66-70, band in 71, exactly as
        # IASC's Astrometrica reports write them ("21.2 G")
        with_mag = Table(candidates; mag=[21.2345, NaN, 19.0])
        mlines = split(chomp(mpc80_report(with_mag, "F51")), "\n")
        @test all(l -> length(l) == 80, mlines)
        @test mlines[1][66:71] == "21.2 G"
        @test mlines[2][66:71] == " "^6                            # NaN -> blank
        @test mlines[3][66:71] == "19.0 G"
        @test split(chomp(mpc80_report(with_mag, "F51"; band='V')), "\n")[1][71] == 'V'
        @test mlines[1][1:65] == split(chomp(mpc80_report(candidates, "F51")), "\n")[1][1:65]

        # date/RA/Dec parsed back and compared within the format's precision
        date_field = line1[16:32]
        @test date_field[1:4] == "2000"
        @test date_field[6:7] == "01"
        @test parse(Float64, date_field[9:end]) ≈ 1.5 atol=1e-6  # J2000.0 = Jan 1.5

        ra_field = line1[33:44]
        ra_deg = 15 * (parse(Float64, ra_field[1:2]) + parse(Float64, ra_field[4:5]) / 60 +
                        parse(Float64, ra_field[7:end]) / 3600)
        # atol matches the field's own precision (0.001s of time ≈
        # 4.2e-6°), not tighter than what the format can represent
        @test ra_deg ≈ 150.123456789 atol=1e-5

        dec_field = line1[45:56]
        dec_sign = dec_field[1] == '-' ? -1 : 1
        dec_deg = dec_sign * (parse(Float64, dec_field[2:3]) + parse(Float64, dec_field[5:6]) / 60 +
                               parse(Float64, dec_field[8:end]) / 3600)
        @test dec_deg ≈ 20.987654321 atol=1e-5

        # both rows sharing id=1 must share the same temp designation
        @test strip(lines[2][6:12]) == strip(line1[6:12])
        @test strip(lines[3][6:12]) != strip(line1[6:12])  # id=2 gets a distinct one

        # sign is mandatory even for a negative declination (row 3, dec=-10.25)
        @test lines[3][45] == '-'

        # temp designation length limit: mirrors ades_psv's trkSub cap,
        # but one character tighter (7, not 8) since this format's
        # designation slot has no eighth column to spare
        @test_throws ArgumentError mpc80_report(candidates, "I41"; trksub_prefix="TOOLONG")
    end

    @testset "digest2_score" begin
        # digest2 is an external C program: skipped when it isn't installed
        digest2_path = get(ENV, "DIGEST2_PATH", "digest2")
        resolved = isfile(digest2_path) ? digest2_path : Sys.which(digest2_path)

        candidates = Table(id=[1, 1, 1], frame=[1, 2, 3], x=[1.0, 2.0, 3.0], y=[1.0, 2.0, 3.0],
                            ra=[150.0, 150.05, 150.10], dec=[20.0, 20.02, 20.04],
                            epoch=[2451545.0, 2451545.04, 2451545.08])

        @test_throws ArgumentError digest2_score(candidates, "XX")  # not 3 characters
        @test_throws ArgumentError digest2_score(candidates, "I41"; digest2_path="not-a-real-digest2-binary")

        if resolved === nothing
            @test_skip "digest2 not installed"
        else
            scores = digest2_score(candidates, "I41"; digest2_path=resolved)
            @test length(scores) == 1
            @test scores.id[1] == 1
            @test 0.0 <= scores.neo_score[1] <= 100.0  # a real score, not this test asserting a specific value
        end
    end

    @testset "link_across_nights" begin
        T0 = 2451545.0
        RATE_RA, RATE_DEC = 33.83, 18.0  # arcsec/day, in the tangent plane around (150,20)

        # a real object: exactly linear motion over 3 nights, so the cross-night
        # extrapolation has zero residual
        function object_row(t)
            dt = t - T0
            dec = 20.0 + RATE_DEC * dt / 3600
            ra = 150.0 + RATE_RA * dt / (cosd(20.0) * 3600)
            return (ra=ra, dec=dec, epoch=t)
        end
        night1_obj = [object_row(T0 - 0.05), object_row(T0 + 0.05)]
        night2_obj = [object_row(T0 + 4.95), object_row(T0 + 5.05)]
        night3_obj = [object_row(T0 + 9.95), object_row(T0 + 10.05)]

        # noise tracklets: consistent motion within their night, at an unrelated
        # position and rate
        noise(t0, ra0, dec0) = [(ra=ra0, dec=dec0, epoch=t0), (ra=ra0 + 0.01, dec=dec0 - 0.01, epoch=t0 + 0.1)]
        night1_noise = noise(T0, 200.0, -10.0)
        night2_noise = noise(T0 + 5.0, 44.0, 60.0)
        night3_noise = noise(T0 + 10.0, 300.0, 0.0)

        to_table(obj, noise) = Table(
            id=vcat(fill(1, length(obj)), fill(2, length(noise))),
            frame=vcat(1:length(obj), 1:length(noise)),
            x=zeros(length(obj) + length(noise)), y=zeros(length(obj) + length(noise)),
            ra=vcat([r.ra for r in obj], [r.ra for r in noise]),
            dec=vcat([r.dec for r in obj], [r.dec for r in noise]),
            epoch=vcat([r.epoch for r in obj], [r.epoch for r in noise]))

        night_movers = [to_table(night1_obj, night1_noise), to_table(night2_obj, night2_noise),
                         to_table(night3_obj, night3_noise)]

        groups = link_across_nights(night_movers)
        @test length(groups) == 1  # only the real object spans >= min_nights=2
        real_group = groups[1]
        @test length(real_group) == 3  # one tracklet per night
        @test Set(g.night for g in real_group) == Set([1, 2, 3])
        @test all(g.id == 1 for g in real_group)  # never picks up a noise id=2 tracklet

        # min_nights: raising it above what the real object spans drops it too
        @test isempty(link_across_nights(night_movers; min_nights=4))

        # a real 3" offset from the exact extrapolation is still within
        # the default 5" tolerance, but not a tightened 1" one
        offset_row(t) = merge(object_row(t), (dec=object_row(t).dec + 3.0 / 3600,))
        night2_offset = [offset_row(T0 + 4.95), offset_row(T0 + 5.05)]
        offset_movers = [night_movers[1], to_table(night2_offset, night2_noise), night_movers[3]]

        offset_groups = link_across_nights(offset_movers)
        @test length(offset_groups) == 1
        @test length(offset_groups[1]) == 3  # still spans all 3 nights within the default tolerance

        # tightened past the 3" offset: night 2 drops out, nights 1 and 3 still link
        tight_groups = link_across_nights(offset_movers; match_radius_arcsec=1.0)
        @test length(tight_groups) == 1
        @test Set(g.night for g in tight_groups[1]) == Set([1, 3])
    end

end
