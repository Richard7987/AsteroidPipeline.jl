#=
Chooses and downloads the sets of other telescopes that telescope_test.jl
searches, by rules fixed before any of their images was looked at:

    julia --project=. examples/telescope_data.jl ztf|lco|decam

ZTF (Palomar 48-inch, 1.01"/pixel, 30 s exposures). IRSA's metadata of every
science exposure of CCD 1, quadrant 1 on the nights of October 2023 (JD
2460218.5 to 2460249.5, one query per night). A sequence is a field in one
filter on one night with 4 or more exposures spanning 10 min to 3 h. Its
position is that quadrant's centre; kept if its ecliptic latitude is within
15 deg, then split by galactic latitude: |b| >= 20 deg and |b| < 10 deg (a
crowded field). Three sequences drawn from each with Xoshiro(11), and for
each a quadrant (CCD 1-16, quadrant 1-4) drawn with the same generator. Of
that quadrant's exposures in the sequence, 5 by number (first, last, three
evenly between) or all if fewer. Full quadrant images (sciimg) from IRSA.

Las Cumbres Observatory (0.4 m telescopes, QHY600 cameras, 0.74"/pixel; 1 m
telescopes, Sinistro cameras, 0.39"/pixel). The public archive's frames
reduced by its pipeline (BANZAI, reduction level 91) on the nights of
2024-12-01 to 2024-12-14, none used before. A sequence is one target on one
camera on one night in one filter, exposures of 15 s or more; its exposures
within 2 h of its first, 5 or more spanning 30 min or more, of which 5 by
number (first, last, three evenly between). Sequences drawn in a random
order (Xoshiro(12)) and kept while the first frame's header says the
telescope tracked the stars (RATRACK and DECTRACK 0) and its position is
within 15 deg of the ecliptic, until 3 of 0.4 m and 3 of 1 m telescopes are
kept.

DECam (Blanco 4 m, 0.263"/pixel, one of its 2k x 4k CCDs per set). NOIRLab's
Astro Data Archive, frames calibrated by the community pipeline (InstCal,
images) of the nights of March 2023. A sequence is one filter on one night,
exposures of 30 s or more whose centres are within 30" of the first's; 4 or
more spanning 10 min to 3 h, of which 5 by number as above. Kept if within
15 deg of the ecliptic. Three drawn with Xoshiro(13), and for each one CCD
(extension 1-59) drawn with the same generator; the archive returns that
extension with the primary header.

Frames are kept in data/real/telescopes/<instrument>/<set>/, the selection in
data/real/telescopes/<instrument>/selection.txt.
=#
using CSV, Dates, HTTP, JSON, Printf, Random

const ROOT = joinpath(@__DIR__, "..", "data", "real", "telescopes")
const IRSA = "https://irsa.ipac.caltech.edu/ibe"

galactic_latitude(ra, dec) = asind(sind(dec) * sind(27.12825) + cosd(dec) * cosd(27.12825) * cosd(ra - 192.85948))
ecliptic_latitude(ra, dec) = asind(sind(dec) * cosd(23.4393) - cosd(dec) * sind(23.4393) * sind(ra))

function irsa_query(where; columns="ra,dec,field,ccdid,qid,filtercode,obsjd,filefracday,maglimit")
    url = IRSA * "/search/ztf/products/sci?" *
          HTTP.URIs.escapeuri(Dict("WHERE" => where, "COLUMNS" => columns, "ct" => "csv"))
    body = String(HTTP.get(url; retry=true, retries=6, readtimeout=600).body)
    startswith(body, "<?xml") && error("IRSA: " * first(body, 300))
    CSV.File(IOBuffer(body))
end

spread(v, n) = length(v) <= n ? v : v[unique(round.(Int, range(1, length(v); length=n)))]

function ztf()
    dir = mkpath(joinpath(ROOT, "ztf"))
    rows = NamedTuple[]
    for jd in 2460218.5:1:2460249.5
        for r in irsa_query("obsjd BETWEEN $jd AND $(jd + 1) AND ccdid=1 AND qid=1")
            push!(rows, (night=jd, field=r.field, filter=String(r.filtercode), obsjd=r.obsjd, ra=r.ra, dec=r.dec))
        end
    end
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for r in rows
        push!(get!(groups, (r.night, r.field, r.filter), NamedTuple[]), r)
    end
    sequences = NamedTuple[]
    for (key, g) in sort(collect(groups); by=first)
        t = sort([r.obsjd for r in g])
        span = (t[end] - t[1]) * 24
        (length(t) >= 4 && 1 / 6 <= span <= 3) || continue
        ra, dec = g[1].ra, g[1].dec
        abs(ecliptic_latitude(ra, dec)) <= 15 || continue
        b = abs(galactic_latitude(ra, dec))
        stratum = b >= 20 ? :high : b < 10 ? :crowded : continue
        push!(sequences, (night=key[1], field=key[2], filter=key[3], n=length(t), span=span, b=b, stratum=stratum))
    end
    rng = Xoshiro(11)
    chosen = NamedTuple[]
    for stratum in (:high, :crowded)
        pool = filter(s -> s.stratum == stratum, sequences)
        println(stratum, ": ", length(pool), " sequences")
        for s in pool[randperm(rng, length(pool))[1:3]]
            push!(chosen, merge(s, (ccd=rand(rng, 1:16), qid=rand(rng, 1:4))))
        end
    end
    open(joinpath(dir, "selection.txt"), "w") do io
        for s in chosen
            q = irsa_query("obsjd BETWEEN $(s.night) AND $(s.night + 1) AND field=$(s.field) AND ccdid=$(s.ccd) " *
                           "AND qid=$(s.qid) AND filtercode='$(s.filter)'")
            exposures = spread(sort(collect(q); by=r -> r.obsjd), 5)
            name = @sprintf("f%04d_c%02d_q%d_%s_%d", s.field, s.ccd, s.qid, s.filter, Int(s.night - 0.5))
            @printf(io, "%s night %.1f n %d span %.2f h |b| %.1f %s maglimit %s\n", name, s.night, s.n, s.span, s.b,
                    s.stratum, join([@sprintf("%.2f", r.maglimit) for r in exposures], ","))
            setdir = mkpath(joinpath(dir, name))
            for r in exposures
                ffd = string(r.filefracday)
                file = @sprintf("ztf_%s_%06d_%s_c%02d_o_q%d_sciimg.fits", ffd, s.field, s.filter, s.ccd, s.qid)
                out = joinpath(setdir, file)
                isfile(out) && continue
                url = "$IRSA/data/ztf/products/sci/$(ffd[1:4])/$(ffd[5:8])/$(ffd[9:end])/$file"
                write(out * ".part", HTTP.get(url; retry=true, retries=6, readtimeout=600).body)
                mv(out * ".part", out)
            end
            println(name, ": ", length(exposures), " exposures")
        end
    end
end

const LCO = "https://archive-api.lco.global/frames/"
lco_json(url) = JSON.parse(String(HTTP.get(url; retry=true, retries=8, read_idle_timeout=300).body))

function lco()
    dir = mkpath(joinpath(ROOT, "lco"))
    frames = Any[]
    for day in Date(2024, 12, 1):Day(1):Date(2024, 12, 14)
        url = LCO * "?" * HTTP.URIs.escapeuri(Dict("public" => "true", "reduction_level" => "91",
              "configuration_type" => "EXPOSE", "start" => string(day), "end" => string(day + Day(1)), "limit" => "100"))
        while url !== nothing
            page = lco_json(url)
            append!(frames, filter(r -> r["exposure_time"] >= 15 && r["observation_day"] == string(day), page["results"]))
            url = page["next"]
        end
    end
    t(r) = DateTime(first(r["observation_date"], 19))
    groups = Dict{Tuple,Vector{Any}}()
    for r in frames
        push!(get!(groups, (r["target_name"], r["instrument_id"], r["observation_day"], r["primary_optical_element"],
                            r["telescope_id"]), Any[]), r)
    end
    sequences = Any[]
    for (key, g) in sort(collect(groups); by=first)
        g = sort(g; by=t)
        g = filter(r -> t(r) - t(g[1]) <= Hour(2), g)
        (length(g) >= 5 && t(g[end]) - t(g[1]) >= Minute(30)) || continue
        push!(sequences, (key=key, frames=spread(g, 5)))
    end
    rng = Xoshiro(12)
    kept = Dict("0m4" => Any[], "1m0" => Any[])
    for s in sequences[randperm(rng, length(sequences))]
        class = first(s.key[5], 3)
        haskey(kept, class) && length(kept[class]) < 3 || continue
        h = lco_json(LCO * "$(s.frames[1]["id"])/headers/")["data"]
        tracked = Float64(get(h, "RATRACK", NaN)) == 0 && Float64(get(h, "DECTRACK", NaN)) == 0
        ra, dec = Float64(h["CRVAL1"]), Float64(h["CRVAL2"])
        (tracked && abs(ecliptic_latitude(ra, dec)) <= 15) || continue
        push!(kept[class], s)
    end
    open(joinpath(dir, "selection.txt"), "w") do io
        println("$(length(sequences)) sequences")
        for class in ("0m4", "1m0"), s in kept[class]
            target, camera, day, band = s.key
            name = replace("$(target)_$(camera)_$(day)_$(band)", r"[^A-Za-z0-9_.-]" => "_")
            println(io, name, " ", class, " ", join([r["basename"] for r in s.frames], ","))
            setdir = mkpath(joinpath(dir, name))
            for r in s.frames
                out = joinpath(setdir, r["basename"] * ".fits.fz")
                isfile(out) && continue
                write(out * ".part", HTTP.get(r["url"]; retry=true, retries=8, read_idle_timeout=300).body)
                mv(out * ".part", out)
            end
            println(name, ": ", length(s.frames), " exposures")
        end
    end
end
const NOIRLAB = "https://astroarchive.noirlab.edu/api"

function decam()
    dir = mkpath(joinpath(ROOT, "decam"))
    body = JSON.json(Dict("outfields" => ["md5sum", "caldat", "dateobs_min", "ra_center", "dec_center", "ifilter", "exposure"],
                          "search" => [["instrument", "decam"], ["proc_type", "instcal"], ["prod_type", "image"],
                                       ["obs_type", "object"], ["caldat", "2023-03-01", "2023-03-31"]]))
    found = JSON.parse(String(HTTP.post("$NOIRLAB/adv_search/find/?limit=100000", ["Content-Type" => "application/json"],
                                        body; retry=true, retries=6, readtimeout=900).body))[2:end]
    t(r) = DateTime(first(r["dateobs_min"], 23))
    exposures = sort(filter(r -> r["exposure"] >= 30 && r["ra_center"] !== nothing, found); by=t)
    sep(a, b) = 3600 * acosd(clamp(sind(a["dec_center"]) * sind(b["dec_center"]) + cosd(a["dec_center"]) *
                                   cosd(b["dec_center"]) * cosd(a["ra_center"] - b["ra_center"]), -1, 1))
    used = falses(length(exposures))
    sequences = Any[]
    for i in eachindex(exposures)
        used[i] && continue
        a = exposures[i]
        idx = [j for j in i:length(exposures) if !used[j] && exposures[j]["caldat"] == a["caldat"] &&
               exposures[j]["ifilter"] == a["ifilter"] && sep(a, exposures[j]) <= 30 && t(exposures[j]) - t(a) <= Hour(3)]
        used[idx] .= true
        g = exposures[idx]
        (length(g) >= 4 && t(g[end]) - t(g[1]) >= Minute(10)) || continue
        abs(ecliptic_latitude(a["ra_center"], a["dec_center"])) <= 15 || continue
        push!(sequences, g)
    end
    rng = Xoshiro(13)
    open(joinpath(dir, "selection.txt"), "w") do io
        println(length(sequences), " sequences")
        for g in sequences[randperm(rng, length(sequences))[1:3]]
            ext = rand(rng, 1:59)
            g = spread(g, 5)
            band = first(split(g[1]["ifilter"]))
            name = "$(g[1]["caldat"])_$(band)_$(round(Int, g[1]["ra_center"] * 10))_ccd$(ext)"
            println(io, name, " ", join([r["md5sum"] for r in g], ","))
            setdir = mkpath(joinpath(dir, name))
            for r in g
                out = joinpath(setdir, r["md5sum"] * ".fits")
                isfile(out) && continue
                write(out * ".part", HTTP.get("$NOIRLAB/retrieve/$(r["md5sum"])/?hdus=0,$ext"; retry=true, retries=6,
                                              readtimeout=600).body)
                mv(out * ".part", out)
            end
            println(name, ": ", length(g), " exposures")
        end
    end
end
instrument = get(ARGS, 1, "")
instrument == "ztf" ? ztf() : instrument == "lco" ? lco() : instrument == "decam" ? decam() :
    error("usage: telescope_data.jl ztf|lco|decam")
