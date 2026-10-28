#=
Validates moving_light_curve on real ZTF frames of (624) Hektor, a Jupiter
trojan with a published rotation period of 6.9205 h, in r from July to
December 2022. The first run fetches the frames from IRSA (no
authentication) into data/real/hektor: about 110 cutouts of 600″.

Declared before running: the photometry within a night should scatter as
much as its reported errors. ZTF visits an object a few times a night, so
the period is expected to fit better than random ones, not to come out of
the periodogram by itself.
=#
using AsteroidPipeline, HTTP, Printf, Random, Statistics

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "hektor")
const IRSA = "https://irsa.ipac.caltech.edu"
const PERIOD = 6.9205 / 24

function read_ipac(text)
    lines = split(text, '\n')
    header = first(l for l in lines if startswith(l, "|"))
    bars = findall('|', header)
    names = [strip(header[bars[i]+1:bars[i+1]-1]) for i in 1:length(bars)-1]
    data = [l for l in lines if !isempty(strip(l)) && !startswith(l, '\\') && !startswith(l, '|')]
    return [Dict(names[i] => String(strip(l[bars[i]:min(bars[i+1], lastindex(l))], [' ', '|']))
                 for i in eachindex(names)) for l in data]
end

mkpath(DATA_DIR)
most = joinpath(DATA_DIR, "most.tbl")
if !isfile(most)
    println("Querying IRSA's moving object search for Hektor...")
    url = "$IRSA/cgi-bin/MOST/nph-most?catalog=ztf&input_type=name_input&obj_name=624&obj_type=Asteroid" *
          "&ephem_type=0&output_mode=Brief&obs_begin=2022-07-01&obs_end=2022-12-31"
    write(most, HTTP.get(url; read_idle_timeout=600).body)
end

# one quadrant per exposure: the one with Hektor nearest its centre, so the cutout isn't clipped
exposures = Dict{String,Dict{String,String}}()
for r in read_ipac(read(most, String))
    r["filtercode"] == "zr" || continue
    prev = get(exposures, r["expid"], nothing)
    if prev === nothing || parse(Float64, r["dist_ctr"]) < parse(Float64, prev["dist_ctr"])
        exposures[r["expid"]] = r
    end
end
obs = sort(collect(values(exposures)); by=r -> parse(Float64, r["obsjd"]))

function cutout(r)
    path = joinpath(DATA_DIR, r["expid"] * ".fits")
    isfile(path) && return path
    ff = r["filefracday"]
    name = @sprintf("ztf_%s_%06d_%s_c%02d_%s_q%s_sciimg.fits", ff, parse(Int, r["field"]), r["filtercode"],
                    parse(Int, r["ccdid"]), r["imgtypecode"], r["qid"])
    url = "$IRSA/ibe/data/ztf/products/sci/$(ff[1:4])/$(ff[5:8])/$(ff[9:end])/$name" *
          "?center=$(r["ra_obj"]),$(r["dec_obj"])&size=600arcsec&gzip=false"
    body = try
        HTTP.get(url; retry=true, retries=3, read_idle_timeout=120).body
    catch e
        @warn "cutout not fetched" r["expid"] exception = e
        return nothing
    end
    startswith(String(body[1:min(6, end)]), "SIMPLE") || return nothing
    write(path * ".part", body)
    mv(path * ".part", path)
    return path
end

println("Fetching $(length(obs)) cutouts...")
paths = asyncmap(cutout, obs; ntasks=4)
have = paths .!== nothing
@printf("%d of %d exposures available\n", count(have), length(obs))
obs, paths = obs[have], String.(paths[have])

night(r) = floor(Int, parse(Float64, r["obsjd"]) - 0.25)
rows = NamedTuple[]
for n in unique(night.(obs))
    k = findall(r -> night(r) == n, obs)
    positions = [(ra=parse(Float64, obs[i]["ra_obj"]), dec=parse(Float64, obs[i]["dec_obj"]),
                  epoch=parse(Float64, obs[i]["obsjd"])) for i in k]
    lc = moving_light_curve(paths[k], positions)
    for m in lc
        push!(rows, (; m.epoch, m.mag, m.mag_err, m.star_nearby, night=n))
    end
end
@printf("%d frames over %d nights; %d near a star, %d not measured\n", length(rows), length(unique(night.(obs))),
        count(r -> r.star_nearby, rows), count(r -> !isfinite(r.mag), rows))
t = [r for r in rows if isfinite(r.mag) && !r.star_nearby]
# ZTF observes from Palomar (I41)
reduced = reduced_magnitudes([r.epoch for r in t], [r.mag for r in t],
                             object_geometry(624, [r.epoch for r in t]; observatory="I41"))

# precision without a rotation model: scatter about a straight line through each night's exposures
within = Float64[]
for n in unique(r.night for r in t)
    s = [r for r in t if r.night == n]
    length(s) >= 3 || continue
    X = hcat(ones(length(s)), [r.epoch - s[1].epoch for r in s])
    y = [r.mag for r in s]
    append!(within, (y .- X * (X \ y)) .* sqrt(length(s) / (length(s) - 2)))
end
@printf("within a night: %d points, RMS about a line %.3f mag; reported errors %.3f–%.3f (10th–90th percentile)\n",
        length(within), sqrt(mean(within .^ 2)), quantile([r.mag_err for r in t], [0.1, 0.9])...)

function fold_rms(k, P)
    time, mag = reduced.epoch[k], reduced.mag[k]
    ω = 2π / P
    X = hcat(ones(length(k)), sin.(ω .* time), cos.(ω .* time), sin.(2ω .* time), cos.(2ω .* time))
    return sqrt(mean((mag .- X * (X \ mag)) .^ 2))
end

alias = 2 / (2 / PERIOD - 1)
Random.seed!(1)
random_periods = (2 .+ 18 .* rand(2000)) ./ 24
for (label, s) in [("Jul–Dec", eachindex(t)), ("Jul–Aug", findall(r -> r.epoch < 2459823.5, t))]
    others = [fold_rms(s, P) for P in random_periods]
    @printf("%s, %d points: published %.3f mag, alias %.4f h %.3f, random periods median %.3f; %.1f%% fit better than the published\n",
            label, length(s), fold_rms(s, PERIOD), 24alias, fold_rms(s, alias), median(others),
            100 * mean(others .< fold_rms(s, PERIOD)))
end

p = recover_rotation_period(reduced.epoch, reduced.mag; minimum_period=2 / 24, maximum_period=20 / 24)
@printf("periodogram's highest peak: %.3f h (half the period: %.3f h; its daily alias: %.3f h)\n",
        24p.period, 12PERIOD, 24 / (2 / PERIOD - 1))
