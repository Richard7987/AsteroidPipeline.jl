using AsteroidPipeline
using TypedTables
using HTTP
using WCS
using FITSIO
using Random
using Test
using Statistics
using Reproject

@testset "AsteroidPipeline.jl" begin

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

        noise_only = 100.0 .+ 5.0 .* randn(32, 32)
        @test length(detect_sources(noise_only; threshold=100.0)) == 0
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

        # velocity refit: frame 2's position has centroiding-like noise
        # ((+0.5,-0.3) off true (15,7)), so the 2-point seed velocity
        # (5.5,-3.3) drifts from the true (5,-3) enough that its
        # extrapolation misses frame 5 by 2.33 px — outside match_radius.
        # Only a refit (least squares over the seed's own matched points,
        # which pulls the fit back toward the 4 exact points) recovers it.
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

        # If the refit's own velocity would exceed max_speed, the refit is
        # discarded and the tracklet keeps its original (already
        # max_speed-checked) seed points, rather than being dropped
        # entirely — a seed under-estimating speed (frame 2 offset toward
        # frame 1: (-1.0,+0.6) off true (15,7)) passes max_speed=5.0 at
        # ~4.66, but refitting with the exact frame 3 point pulls the
        # velocity back toward the true ~5.83, which exceeds it.
        under1 = Table(x=[10.0], y=[10.0])
        under2 = Table(x=[14.0], y=[7.6])
        under3 = Table(x=[20.0], y=[4.0])
        under_tracklets = link_candidates([under1, under2, under3], [0.0, 1.0, 2.0];
                                           match_radius=3.0, min_frames=1, max_speed=5.0)
        @test length(under_tracklets) == 1
        @test length(under_tracklets[1]) == 3
        @test under_tracklets[1][2].x == 14.0 && under_tracklets[1][2].y == 7.6  # unrefit seed point kept
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

        # NaN marks out-of-footprint pixels near the edge (real, not a
        # bug); argmax must not be fooled by it — Julia's `isless(NaN, x)`
        # is false for every x, so a raw `argmax` over a NaN-containing
        # array returns the first NaN it meets, not the true peak.
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
        # empirical FWHM (in pixels, from the recovered radial profile)
        # should match the injected Gaussian's, sigma * 2sqrt(2ln2). Note
        # findfirst over a range returns a *position*, not the range's
        # value at that position (off by one here since the range starts
        # at 0) — e.g. findfirst(iseven, 0:9) is 1 (0 is at position 1),
        # not 0.
        half_max_pos = findfirst(r -> psf[c+r, c] < psf[c, c] / 2, 0:(c - 1))
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

        # pure Gaussian noise, no sources -> std(S_corr) ~ 1. gain is set
        # very high to isolate the FFT/normalization chain being tested
        # here from _variance_map's Poisson-term heuristic (which, given
        # only noise and no real source model, would otherwise inflate
        # the variance estimate with spurious "signal" from noise
        # excursions above the median — a known approximation, not what
        # this test is checking).
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

        # asymmetric case: different PSF widths and different noise levels
        # for N and R, as real science-vs-stacked-reference frames always
        # have — the symmetric checks above cannot exercise the general
        # (f_r^2 f_n, f_r f_n^2) terms in _k_n/_k_r, since psf_n == psf_r
        # cancels most of the asymmetry out algebraically.
        gaussian_psf(s) = (p = [exp(-((i - 5)^2 + (j - 5)^2) / (2s^2)) for i in 1:9, j in 1:9]; p ./ sum(p))
        psf_n2, psf_r2 = gaussian_psf(1.5), gaussian_psf(2.2)
        sigma_n2, sigma_r2 = 5.0, 2.0
        Random.seed!(15)
        n_noise2 = 100.0 .+ sigma_n2 .* randn(sz)
        r_noise2 = 100.0 .+ sigma_r2 .* randn(sz)
        _, s_corr_asym = zogy_subtract(n_noise2, r_noise2; psf_n=psf_n2, psf_r=psf_r2,
                                        sigma_n=sigma_n2, sigma_r=sigma_r2, gain_n=1e8, gain_r=1e8)
        @test isapprox(std(s_corr_asym), 1.0; atol=0.15)

        # mismatched sky background between N and R (e.g. different lunar
        # phase/airglow on different nights — routine in real data, and
        # unrelated to build_reference's *photometric* zeropoint matching,
        # which aligns star flux scales, not sky brightness). Regression
        # test for a real bug found on real ZTF data: without subtracting
        # each image's own background first, a raw FFT's DC term is a
        # *sum* over the whole image (~10^6 pixels here), so even a small
        # background mismatch blew up into a near-constant offset swamping
        # S_corr almost everywhere (observed: mean ~35-90 instead of ~0,
        # and >99.9% of a real frame reading "above 6 sigma").
        Random.seed!(16)
        n_bg_noise = 400.0 .+ sigma .* randn(sz)   # science: brighter sky
        r_bg_noise = 100.0 .+ sigma .* randn(sz)   # reference: darker sky
        _, s_corr_bg = zogy_subtract(n_bg_noise, r_bg_noise; psf_n=psf, psf_r=psf,
                                      sigma_n=sigma, sigma_r=sigma, gain_n=1e8, gain_r=1e8)
        @test isapprox(mean(s_corr_bg), 0.0; atol=0.2)
        @test isapprox(std(s_corr_bg), 1.0; atol=0.15)
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

            # science triplet: same static star, plus a moving source far
            # enough away that it never falls within estimate_psf's
            # isolation radius of the static star used for PSF fitting.
            # Each frame is also dithered a few pixels from the others
            # (real same-night ZTF dither is ~5-6 px — see src/reference.jl),
            # a regression case for reprojection's resulting NaN border:
            # left unsanitized before zogy_subtract, a single NaN poisons
            # its whole-array median/fft and turns all of S_corr to NaN.
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

            # quality_max_std=Inf: this test is about NaN-border
            # sanitization under dithering, not the quality gate — at this
            # image's small (100x100) size, the gate's plain-std metric
            # would otherwise be skewed by this test's own genuine
            # injected source (see "run_pipeline quality gate" below,
            # which is sized specifically to test the gate correctly).
            candidates = run_pipeline(scipaths; threshold=6.0, match_radius=2.0,
                                       reference=reference,
                                       psf_threshold=20.0, psf_min_separation=15.0,
                                       quality_max_std=Inf)

            @test length(unique(candidates.id)) >= 1
            @test !any(isnan, candidates.x) && !any(isnan, candidates.ra)
            by_frame = Dict(row.frame => row for row in candidates if row.id == candidates[1].id)
            @test haskey(by_frame, 1) && haskey(by_frame, 3)
            if haskey(by_frame, 1) && haskey(by_frame, 3)
                # expected positions on reference.wcs's grid, correcting
                # for the injected star's own frame's CRPIX dither
                # (dither(mjd, 4.0k, -3.0k) above) — the same physical sky
                # position lands at a different pixel index once the two
                # frames' pixel grids no longer share one CRPIX.
                @test by_frame[1].x ≈ x0 atol=1.5
                @test by_frame[1].y ≈ y0 atol=1.5
                @test by_frame[3].x ≈ x0 + 2dx - 4.0 * 2 atol=1.5
                @test by_frame[3].y ≈ y0 + 2dy - (-3.0 * 2) atol=1.5
            end
        end
    end

    @testset "run_pipeline quality gate" begin
        mktempdir() do dir
            # 600x600, not the 100x100 used elsewhere in this file: this
            # gate uses plain Statistics.std (see src/pipeline.jl for why
            # — a MAD-based spread turned out blind to the real anomaly it
            # exists to catch), and a plain std over a *small* image is
            # skewed by wherever a genuine bright source happens to be,
            # which would wrongly penalize the very detection this
            # pipeline exists to make. At 100x100 a real ~60 sigma source
            # pushed a clean frame's std to ~4; at 600x600 (the same
            # source, now a much smaller fraction of the image, closer to
            # real ~10^6 px ZTF data) it stays ~1.1, correctly separated
            # from the deliberately-bad frame's ~3.5.
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

            # Frames 1-2: clean, plus a genuine bright transient (250 flux,
            # present in science only) — must NOT trip the gate. Frame 3:
            # clean pixel noise plus a broad, low-contrast "glow" patch — a
            # stand-in for the real anomalous-frame case in git log
            # `3e31b95` (a spatially localized issue, not simply higher
            # uniform noise, which zogy_subtract's sigma_n already corrects
            # for and so would *not* trip this gate).
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

    @testset "crossmatch_catalog" begin
        @test_throws ArgumentError crossmatch_catalog([], :unknown_catalog; radius=5.0)

        try
            vega = [(id=1, ra=279.23473479, dec=38.78368896)]
            matches = crossmatch_catalog(vega, :simbad; radius=5.0)
            @test length(matches) >= 1
            @test matches[1].id == 1
        catch e
            e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
            @test_skip "network unavailable"
        end

        try
            # A known object (127319 "2002 JB99") ~134" from this position
            # at this epoch (verified independently against SkyBoT). This
            # is a positive control: it caught a real bug where the Julian
            # Date epoch, formatted with Julia's default Float64 `string`,
            # rendered in scientific notation (e.g. "2.46e6") and SkyBoT
            # silently rejected it as empty, making every match vanish
            # regardless of true content — an all-empty test wouldn't have
            # caught this, since it can't distinguish "correctly empty"
            # from "silently broken".
            known = [(id=1, ra=36.521300910000, dec=2.127107410000, epoch=2459288.617372700)]
            matches = crossmatch_catalog(known, :skybot; radius=200.0)
            @test any(==("2002 JB99"), matches.name)
        catch e
            e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
            @test_skip "network unavailable"
        end

        try
            # near the celestial pole, an epoch/radius/id chosen so no
            # minor planet falls inside: regression test for an empty
            # result crashing Table construction (found via real ZTF data).
            empty_candidates = [(id=1, ra=0.0, dec=89.9, epoch=2459288.6174)]
            matches = crossmatch_catalog(empty_candidates, :skybot; radius=1.0)
            @test length(matches) == 0
            @test isempty(matches.id)
        catch e
            e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
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

        # _astrometry_poll: succeeds once `check` stops returning nothing,
        # and times out (erroring) if it never does. (_astrometry_submission_job
        # and _astrometry_job_done both make real HTTP requests, so aren't
        # unit-testable without network — only _astrometry_poll's own
        # retry/timeout logic, generic over any `check` closure, is.)
        calls = Ref(0)
        result = AsteroidPipeline._astrometry_poll(0.3, 0.05, "test condition") do
            calls[] += 1
            calls[] >= 3 ? :done : nothing
        end
        @test result == :done
        @test calls[] == 3
        @test_throws ErrorException AsteroidPipeline._astrometry_poll(
            () -> nothing, 0.1, 0.05, "test condition")

        # live end-to-end round trip against the real service — only runs
        # if a real API key is available; validated once against a real
        # ZTF frame (see INVESTIGATION_LOG.md), not re-run on every CI pass
        api_key = get(ENV, "ASTROMETRY_API_KEY", nothing)
        if api_key === nothing
            @test_skip "ASTROMETRY_API_KEY not set"
        else
            mktempdir() do dir
                nx, ny = 60, 60
                wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                    cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])
                img = 100.0 .+ 5.0 .* randn(nx, ny)
                for i in 1:nx, j in 1:ny
                    r2 = (i - 30)^2 + (j - 30)^2
                    img[i, j] += 2000.0 * exp(-r2 / (2 * 2.0^2))
                end
                path = joinpath(dir, "test.fits")
                FITS(path, "w") do f
                    write(f, img)  # deliberately no WCS: this is what plate_solve is for
                end
                solved = plate_solve(path; api_key=api_key)
                @test solved isa WCSTransform
            end
        end
    end

end
