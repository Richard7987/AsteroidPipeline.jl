#=
Rotation periods of several asteroids from public ZTF images, against their
published ones: the same steps as moving_light_curve_demo.jl, for more objects.

    julia --project=. examples/ztf_periods.jl

Declared before running, with the objects and periods chosen from the
Asteroid Lightcurve Database (LCDB, Warner et al. 2009, rev. 2023-October, as
JPL's Small-Body Database gives it), without looking at their ZTF data:

- objects: four Jupiter trojans, two main-belt asteroids and two near-Earth
  asteroids, with periods of 4-11 h (ZTF visits an object a few times a night);
- images: every ZTF r exposure of 2022 in which the object is predicted at
  V 13.5-19 (brighter saturates; fainter is noise);
- recovered: the Fourier period, doubled (two maxima per rotation), within 1%
  of the published one. Every object is reported, recovered or not.

The first run fetches the frames from IRSA into data/real/ztf_periods.
=#
using AsteroidPipeline, HTTP, Printf, Random, Statistics

const OBJECTS = [(588, "Achilles", 7.306), (884, "Priamus", 6.8605), (911, "Agamemnon", 6.592),
                 (3451, "Mentor", 7.702), (243, "Ida", 4.634), (216, "Kleopatra", 5.385),
                 (1627, "Ivar", 4.795), (1036, "Ganymed", 10.297)]
const IRSA = "https://irsa.ipac.caltech.edu"
const ROOT = joinpath(@__DIR__, "..", "data", "real", "ztf_periods")

function read_ipac(text)
    lines = split(text, '\n')
    header = first(l for l in lines if startswith(l, "|"))
    bars = findall('|', header)
    names = [strip(header[bars[i]+1:bars[i+1]-1]) for i in 1:length(bars)-1]
    data = [l for l in lines if !isempty(strip(l)) && !startswith(l, '\\') && !startswith(l, '|')]
    return [Dict(names[i] => String(strip(l[bars[i]:min(bars[i+1], lastindex(l))], [' ', '|']))
                 for i in eachindex(names)) for l in data]
end

function cutout(r, dir)
    path = joinpath(dir, r["expid"] * ".fits")
    isfile(path) && return path
    ff = r["filefracday"]
    name = @sprintf("ztf_%s_%06d_%s_c%02d_%s_q%s_sciimg.fits", ff, parse(Int, r["field"]), r["filtercode"],
                    parse(Int, r["ccdid"]), r["imgtypecode"], r["qid"])
    url = "$IRSA/ibe/data/ztf/products/sci/$(ff[1:4])/$(ff[5:8])/$(ff[9:end])/$name" *
          "?center=$(r["ra_obj"]),$(r["dec_obj"])&size=600arcsec&gzip=false"
    body = try
        HTTP.get(url; retry=true, retries=8, read_idle_timeout=120).body
    catch e
        @warn "cutout not fetched" r["expid"] exception = e
        return nothing
    end
    startswith(String(body[1:min(6, end)]), "SIMPLE") || return nothing
    write(path * ".part", body)
    mv(path * ".part", path)
    return path
end

function exposures(number, dir)
    most = joinpath(dir, "most.tbl")
    if !isfile(most)
        url = "$IRSA/cgi-bin/MOST/nph-most?catalog=ztf&input_type=name_input&obj_name=$number&obj_type=Asteroid" *
              "&ephem_type=0&output_mode=Brief&obs_begin=2022-01-01&obs_end=2022-12-31"
        write(most, HTTP.get(url; read_idle_timeout=600).body)
    end
    best = Dict{String,Dict{String,String}}()
    for r in read_ipac(read(most, String))
        r["filtercode"] == "zr" && 13.5 <= parse(Float64, r["vmag"]) <= 19 || continue
        prev = get(best, r["expid"], nothing)
        if prev === nothing || parse(Float64, r["dist_ctr"]) < parse(Float64, prev["dist_ctr"])
            best[r["expid"]] = r
        end
    end
    return sort(collect(values(best)); by=r -> parse(Float64, r["obsjd"]))
end

function fold_rms(time, mag, P)
    ω = 2π / P
    X = hcat(ones(length(time)), sin.(ω .* time), cos.(ω .* time), sin.(2ω .* time), cos.(2ω .* time))
    return sqrt(mean((mag .- X * (X \ mag)) .^ 2))
end

night(r) = floor(Int, parse(Float64, r["obsjd"]) - 0.25)

function measure(number, name, published_h)
    dir = mkpath(joinpath(ROOT, string(number)))
    obs = exposures(number, dir)
    paths = map(r -> cutout(r, dir), obs)
    have = paths .!== nothing
    obs, paths = obs[have], String.(paths[have])
    rows = NamedTuple[]
    for n in unique(night.(obs))
        k = findall(r -> night(r) == n, obs)
        positions = [(ra=parse(Float64, obs[i]["ra_obj"]), dec=parse(Float64, obs[i]["dec_obj"]),
                      epoch=parse(Float64, obs[i]["obsjd"])) for i in k]
        lc = try
            moving_light_curve(paths[k], positions)
        catch e
            @warn "night not measured" number n exception = e
            continue
        end
        foreach(m -> push!(rows, (; m.epoch, m.mag, m.star_nearby)), lc)
    end
    t = [r for r in rows if isfinite(r.mag) && !r.star_nearby]
    if length(t) < 30
        @printf("RES (%d) %s: %d exposures, %d points: too few\n", number, name, length(obs), length(t))
        return
    end
    # ZTF observes from Palomar (I41)
    reduced = reduced_magnitudes([r.epoch for r in t], [r.mag for r in t],
                                 object_geometry(number, [r.epoch for r in t]; observatory="I41"))
    open(joinpath(dir, "reduced.csv"), "w") do io
        println(io, "epoch,mag")
        foreach(k -> @printf(io, "%.6f,%.4f\n", reduced.epoch[k], reduced.mag[k]), eachindex(reduced.epoch))
    end
    P = published_h / 24
    f = fourier_period(reduced.epoch, reduced.mag; minimum_period=2 / 24, maximum_period=20 / 24, order=2)
    Random.seed!(1)
    others = [fold_rms(reduced.epoch, reduced.mag, p) for p in (2 .+ 18 .* rand(2000)) ./ 24]
    off = abs(2f.period - P) / P
    @printf("RES (%d) %s: %d exposures, %d points over %d nights; Fourier %.4f h, doubled %.4f ± %.4f h, published %.4f h (%.2f%% off: %s); amplitude %.2f mag; %.1f%% of random periods fit better than the published\n",
            number, name, length(obs), length(t), length(unique(night.(obs))), 24f.period, 48f.period,
            48f.period_error, published_h, 100off, off <= 0.01 ? "recovered" : "not recovered", f.amplitude,
            100 * mean(others .< fold_rms(reduced.epoch, reduced.mag, P)))
end

for (number, name, period) in OBJECTS
    try
        measure(number, name, period)
    catch e
        @printf("RES (%d) %s: failed: %s\n", number, name, sprint(showerror, e)[1:min(end, 200)])
    end
    flush(stdout)
end
