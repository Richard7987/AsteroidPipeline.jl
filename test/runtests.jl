using AsteroidPipeline
using TypedTables
using HTTP
using WCS
using Random
using Test

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
    end

end
