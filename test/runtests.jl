using AsteroidPipeline
using TypedTables
using HTTP
using WCS
using FITSIO
using Random
using Test
using Statistics
using Reproject
using Distributed

# Loads .env (repo root, gitignored) into ENV for local test runs — e.g.
# ASTROMETRY_API_KEY, which enables the live plate_solve round-trip test
# below instead of skipping it. Real environment variables always win: a
# line here is only applied if that key isn't already set. Not read by
# src/ itself — plate_solve takes api_key as an explicit argument; this
# is test-only convenience.
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

        # gain=nothing (default) uses background noise only; a finite gain
        # adds the bright source's own Poisson (shot) noise, which must
        # only ever widen flux_err, never narrow it, and only for a real
        # source (a pure-noise frame's flux_err is unaffected since there's
        # no positive signal to contribute shot noise).
        with_gain = detect_sources(image; threshold=5.0, gain=2.0)
        best_gain = with_gain[argmax(with_gain.peak)]
        @test best_gain.flux_err > best.flux_err

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

    @testset "find_variable_sources" begin
        # frame column 1: constant star (flux unchanged every frame);
        # column 2: variable star (flux swings well beyond flux_err);
        # column 3: mover (position advances each frame, so it can never
        # match within position_tolerance past frame 1) — a negative case
        # proving find_variable_sources doesn't also pick up movers.
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
        # 5 constant stars, well above the default S/N floor, whose flux
        # in frame 2 is uniformly 0.8x frame 1's (standing in for a real
        # zeropoint/transparency mismatch between exposures — measured as
        # large as this on real ZTF frames, see INVESTIGATION_LOG.md) plus
        # one genuinely variable source whose flux does *not* follow that
        # scaling; frame 3 returns to frame 1's scale.
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

        # without normalization, frame 2's uniform 0.8x mismatch alone can
        # push a merely-rescaled constant star's chi2 over threshold too —
        # normalization is what keeps that from happening above, so
        # disabling it must not find *fewer* candidates here.
        raw_groups = find_variable_sources([frame1, frame2, frame3], timestamps; normalize=false)
        @test length(raw_groups) >= length(groups)
    end

    @testset "photometric_scale min_stars guard" begin
        # only 2 matched stars per frame — below the default min_stars=5 —
        # so the scale must stay uncorrected (1.0) rather than trust a
        # factor derived from 2 stars (measured on real data: only 2 of
        # 120 detections cleared the default S/N floor in every frame).
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

        # Regression test for a real bug found on real Pan-STARRS1 (IASC
        # practice campaign) headers: CNPIX1/CNPIX2 (a legacy IRAF/DSS
        # plate-astrometry keyword pair, unrelated to this header's real
        # CTYPE/CRVAL/CRPIX/CDELT WCS) makes wcslib build a separate,
        # degenerate implicit WCS and raise "Linear transformation matrix
        # is singular" for the *whole* header, even though the real WCS
        # parses fine alone — confirmed by bisecting a real PS1 header
        # down to this exact keyword pair; see load_wcs's docstring.
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

    @testset "build_reference with workers" begin
        # Same synthetic setup as the sequential test above, run through
        # the `workers` path instead — proves the pmap/WCS-header-string
        # round trip actually works (not just that it doesn't crash) by
        # requiring it to reproduce the sequential result exactly.
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

    @testset "fit_moffat_psf / estimate_psf fallback" begin
        Random.seed!(13)
        nx, ny = 120, 120
        img = 100.0 .+ 3.0 .* randn(nx, ny)
        true_alpha, true_beta = 3.0, 2.5
        # centers all within ~25 px of each other, closer than the default
        # min_separation=40 — every star fails estimate_psf's isolation
        # filter, forcing the analytic fallback. relaxation_attempts=0
        # below keeps this true regardless of the relaxation added for
        # sparser fields (see the "min_separation relaxation" testset) —
        # this test is specifically about the fallback itself, not about
        # how hard estimate_psf tries before reaching it, and these
        # centers are also close enough that a relaxed min_separation
        # would make some pass isolation while still leaving overlapping
        # fit_moffat_psf stamps (stamp_size=25 => half=12, closer than
        # 2*12=24 apart) that would corrupt the analytic fit anyway.
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
        # 25 px apart: fails the default min_separation=40, but 25 >= 20
        # (the first relaxed attempt, 40/2), so both should pass isolation
        # once relaxed — this should recover a real empirical PSF, not
        # fall back to the analytic Moffat fit.
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

            # same real transient, but with min_frames *omitted* entirely
            # (its new default, nothing, auto-subtracts the 1 gated frame
            # from length(scipaths)=3, giving an effective min_frames=2) —
            # must reach the same real result as the explicit min_frames=1
            # above, without the caller having to know a frame was gated.
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

    @testset "search_field" begin
        mktempdir() do dir
            # Larger canvas than the "run_pipeline" testset above needs
            # (80x60): a bright-enough variable star to clear the S/N floor
            # below (amplitude up to 8000) contaminates
            # BackgroundMeshes' background/RMS estimate at that canvas
            # size — found empirically as ~175-450 spurious detections
            # per frame at 80x60, dropping to a clean 2 once the source is
            # a small enough fraction of the frame — the same
            # source-fraction-of-image effect already documented in
            # INVESTIGATION_LOG.md's quality-gate test (600x600 for the
            # same reason). Real ZTF frames (~10^6 px) aren't at risk of
            # this; a small synthetic frame with a very bright source is.
            nx, ny = 300, 225
            wcs = WCSTransform(2; crpix=[nx / 2, ny / 2], crval=[150.0, 20.0],
                                cdelt=[-1 / 3600, 1 / 3600], ctype=["RA---TAN", "DEC--TAN"])

            # true asteroid track (as in the "run_pipeline" testset above)
            x0, y0, dx, dy = 20.0, 45.0, 5.0, -3.0
            # stationary variable star, well separated from the track,
            # with a flux swing far beyond photon noise and flux_err/flux
            # comfortably under find_variable_sources's default 10% S/N
            # floor in every frame.
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
            # a real VSX variable (Gaia DR3 2501737571391422464, type RS)
            # inside the real_data_demo.jl field, verified independently
            # via a direct VizieR TAP query during the migration off the
            # CDS X-Match service (see INVESTIGATION_LOG.md) — a positive
            # control, not just an absence-of-error check.
            star = [(id=1, ra=36.2344, dec=2.06997)]
            matches = crossmatch_catalog(star, :vsx; radius=5.0)
            @test any(==("RS"), matches.class)
            @test matches[1].id == 1
        catch e
            e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
            @test_skip "network unavailable"
        end

        try
            # regression test for batched (OR-chained, one request for all
            # candidates) cross-matching: 3 real, independently-verified
            # positions in one call — 2 real VSX variables (same RS star
            # as above, plus V0651 Ori, type EW) and 1 real non-match (near
            # the celestial pole) — confirms each returned row is
            # attributed to the *correct* candidate id from a single
            # shared request, not just that matches exist somewhere.
            batch = [(id=10, ra=36.2344, dec=2.06997),
                     (id=20, ra=83.1937, dec=5.41603),
                     (id=30, ra=0.0, dec=89.9)]
            matches = crossmatch_catalog(batch, :vsx; radius=5.0)
            @test length(matches) == 2
            @test Set(matches.id) == Set([10, 20])
            @test only(matches[matches.id.==10].class) == "RS"
            @test only(matches[matches.id.==20].class) == "EW"
        catch e
            e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
            @test_skip "network unavailable"
        end

        try
            # near the celestial pole, tiny radius: regression test for
            # empty CDS TAP results (mirrors the SkyBoT empty-result case
            # below) — table construction must not crash on zero rows.
            empty_candidates = [(id=1, ra=0.0, dec=89.9)]
            matches = crossmatch_catalog(empty_candidates, :simbad; radius=1.0)
            @test length(matches) == 0
            @test isempty(matches.id)
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
        # if a real API key is available. Needs a real star field: tried
        # once with a single-Gaussian-source-in-noise synthetic image and
        # the job legitimately failed to solve — astrometry.net matches
        # geometric patterns across several stars, which one fake point
        # source can never provide, no matter how good the key is. Uses
        # a real downloaded ZTF frame instead (its own real WCS is simply
        # never read by plate_solve, which only uploads pixels) — this is
        # the same frame the docstring's one-time manual validation used.
        api_key = get(ENV, "ASTROMETRY_API_KEY", nothing)
        real_frame = joinpath(@__DIR__, "..", "data", "real", "science",
                               "ztf_20191023195590_000451_zr_c01_o_q1_sciimg.fits")
        if api_key === nothing
            @test_skip "ASTROMETRY_API_KEY not set"
        elseif !isfile(real_frame)
            @test_skip "real ZTF test frame not present — run examples/fetch_data.sh"
        else
            solved = plate_solve(real_frame; api_key=api_key)
            @test solved isa WCSTransform
        end
    end

    @testset "julian_date_to_iso8601" begin
        # Two independent, exactly known reference points — not just
        # trusted from the algorithm's own derivation.
        @test AsteroidPipeline.julian_date_to_iso8601(2451545.0) == "2000-01-01T12:00:00.000Z"  # J2000.0
        @test AsteroidPipeline.julian_date_to_iso8601(2440587.5) == "1970-01-01T00:00:00.000Z"  # Unix epoch

        # A fractional-second value that must round, not truncate, to
        # avoid landing at HH:MM:59.999 instead of rolling to the next
        # minute — 2451545.0 + 1 second (in days) rounds up to whole ms.
        one_second = 1.0 / 86400
        @test AsteroidPipeline.julian_date_to_iso8601(2451545.0 + one_second) == "2000-01-01T12:00:01.000Z"
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

end
