#=
Tests moving_light_curve and fourier_period on asteroids chosen by a rule, with
the method used for (624) Hektor (moving_light_curve_demo.jl) and nothing tuned.

Declared before measuring any of them:

Sample. From the Asteroid Lightcurve Database (Warner et al. 2009): numbered,
period quality U 3 or 3-, 2-20 h, largest amplitude at least 0.3 mag, not from a
survey, not binary, not tumbling. Hektor is left out. In three groups (near-Earth,
Jupiter trojan, main belt), shuffled with a fixed seed, the first 10 with at least
60 ZTF r exposures at V 14-19 within 180 days, which is also the window measured.
Only the number of visits is looked at, never a measurement.

Outcome. The fit's period, 2-20 h with 2 harmonics, within 1% of the published one
or of half of it counts as recovered. Twice the period and daily aliases are
reported as what they are, and don't count. Every selected object is reported.

Control. Each object's magnitudes shuffled among its epochs and searched the same way.

The first run fetches the database and the frames from IRSA (no authentication)
into data/real/rotation_sample: about 1.3 MB per exposure.
=#
using AsteroidPipeline, HTTP, Printf, Random, Statistics
import Dates
import p7zip_jll

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "rotation_sample")
const LCDB = "https://minplanobs.org/MPInfo/datazips/LCLIST_PUB_CURRENT.zip"
const IRSA = "https://irsa.ipac.caltech.edu"
const SEED = 20261105
const GROUPS = ["near-Earth" => ["9101"], "trojan" => ["9202"], "main belt" => ["9104", "9105", "9106"]]
const PER_GROUP = 10
const MIN_EXPOSURES = 60
const WINDOW_DAYS = 180
const V_RANGE = (14.0, 19.0)
const TOLERANCE = 0.01
const SIDEREAL_DAY_H = 23.9345

read_table(path) = [Dict(zip(split(first(eachline(path)), ','), split(l, ','))) for l in Iterators.drop(eachline(path), 1)]

# column positions from the summary's header line; names are Latin-1, so by byte
function read_lcdb(path)
    rows = NamedTuple[]
    for line in Iterators.drop(eachline(path), 5)
        b = codeunits(line)
        length(b) < 209 && continue
        col(i, j) = String(strip(String(b[i:j])))
        number, period, amax = tryparse(Int, col(1, 7)), tryparse(Float64, col(146, 158)), tryparse(Float64, col(183, 186))
        (number === nothing || period === nothing || amax === nothing) && continue
        push!(rows, (; number, name=col(11, 40), family=col(63, 70), period, amax,
                     U=col(188, 189), notes=col(191, 195), binary=col(197, 199), survey=col(205, 209)))
    end
    return rows
end

eligible(r) = r.number > 0 && r.number != 624 && r.U in ("3", "3-") && 2 <= r.period <= 20 && r.amax >= 0.3 &&
              isempty(r.survey) && isempty(r.binary) && !occursin('T', r.notes)

function read_ipac(text)
    lines = split(text, '\n')
    header = first(l for l in lines if startswith(l, "|"))
    bars = findall('|', header)
    names = [strip(header[bars[i]+1:bars[i+1]-1]) for i in 1:length(bars)-1]
    data = [l for l in lines if !isempty(strip(l)) && !startswith(l, '\\') && !startswith(l, '|')]
    return [Dict(names[i] => String(strip(l[bars[i]:min(bars[i+1], lastindex(l))], [' ', '|']))
                 for i in eachindex(names)) for l in data]
end

# The months JPL Horizons predicts the object within half a magnitude of V_RANGE, as
# date spans of at most a year: IRSA answers a query of the whole 2018-2024
# survey for a near-Earth asteroid with nothing, after ~7 minutes.
function bright_spans(number)
    epochs = 2458209.5:10:2460675.5                     # 2018-04-01 to 2024-12-31
    q = Dict("COMMAND" => "'$number;'", "EPHEM_TYPE" => "OBSERVER", "CENTER" => "'500@399'", "QUANTITIES" => "'9'",
             "TIME_TYPE" => "UT", "CSV_FORMAT" => "YES", "OBJ_DATA" => "NO")
    bright = Float64[]
    for (chunk, result) in zip(Iterators.partition(epochs, 50), AsteroidPipeline._horizons(q, [@sprintf("%.1f", t) for t in epochs]))
        rows = split(strip(match(r"\$\$SOE(.*)\$\$EOE"s, result)[1]), '\n')
        for (t, row) in zip(chunk, rows)
            v = tryparse(Float64, strip(split(row, ',')[4]))
            v !== nothing && V_RANGE[1] - 0.5 <= v <= V_RANGE[2] + 0.5 && push!(bright, t)
        end
    end
    spans = Tuple{Float64,Float64}[]
    for t in bright
        if !isempty(spans) && t - spans[end][2] <= 10 && spans[end][2] - spans[end][1] < 365
            spans[end] = (spans[end][1], t)
        else
            push!(spans, (t, t))
        end
    end
    return [(a - 10, b + 10) for (a, b) in spans]
end

function most_table(number)
    path = joinpath(DATA_DIR, "most", "$number.tbl")
    if !isfile(path)
        parts = String[]
        for (a, b) in bright_spans(number)
            day(jd) = Dates.format(Dates.julian2datetime(jd), "yyyy-mm-dd")
            url = "$IRSA/cgi-bin/MOST/nph-most?catalog=ztf&input_type=name_input&obj_name=$number&obj_type=Asteroid" *
                  "&ephem_type=0&output_mode=Brief&obs_begin=$(day(a))&obs_end=$(day(b))"
            push!(parts, String(HTTP.get(url; read_idle_timeout=900, retry=false).body))
        end
        # one header, then every span's rows
        header = isempty(parts) ? "" : join(filter(l -> startswith(l, '\\') || startswith(l, '|'), split(first(parts), '\n')), '\n')
        rows = [l for p in parts for l in split(p, '\n') if !isempty(strip(l)) && !startswith(l, '\\') && !startswith(l, '|')]
        write(path * ".part", isempty(rows) ? "" : header * "\n" * join(rows, '\n') * "\n")
        mv(path * ".part", path)
    end
    return read(path, String)
end

# r exposures in the visibility range, one quadrant each (the object nearest its centre, so the cutout isn't clipped)
function exposures(number)
    best = Dict{String,Dict{String,String}}()
    text = most_table(number)
    occursin(r"^\|"m, text) || return Dict{String,String}[]
    for r in read_ipac(text)
        r["filtercode"] == "zr" || continue
        v = tryparse(Float64, r["vmag"])
        (v === nothing || !(V_RANGE[1] <= v <= V_RANGE[2])) && continue
        prev = get(best, r["expid"], nothing)
        if prev === nothing || parse(Float64, r["dist_ctr"]) < parse(Float64, prev["dist_ctr"])
            best[r["expid"]] = r
        end
    end
    return sort(collect(values(best)); by=r -> parse(Float64, r["obsjd"]))
end

function best_window(obs)
    jd = [parse(Float64, r["obsjd"]) for r in obs]
    isempty(jd) && return obs
    counts = [searchsortedlast(jd, t + WINDOW_DAYS) - i + 1 for (i, t) in enumerate(jd)]
    i = argmax(counts)
    return obs[i:i+counts[i]-1]
end

function select_sample()
    path = joinpath(DATA_DIR, "selection.csv")
    isfile(path) && return read_table(path)
    zip_path = joinpath(DATA_DIR, "lcdb.zip")
    isfile(zip_path) || write(zip_path, HTTP.get(LCDB; read_idle_timeout=300).body)
    summary = joinpath(DATA_DIR, "lc_summary_pub.txt")
    isfile(summary) || run(pipeline(`$(p7zip_jll.p7zip()) e -so $zip_path lc_summary_pub.txt`; stdout=summary))
    rows = filter(eligible, read_lcdb(summary))
    log = open(joinpath(DATA_DIR, "selection_log.csv"), "w")
    println(log, "group,number,name,exposures,taken")
    selected = NamedTuple[]
    for (group, families) in GROUPS
        pool = shuffle(MersenneTwister(SEED), sort(filter(r -> r.family in families, rows); by=r -> r.number))
        @printf("%s: %d eligible\n", group, length(pool))
        taken = 0
        for r in pool
            taken >= PER_GROUP && break
            n = try
                length(best_window(exposures(r.number)))
            catch e
                @warn "query failed" r.number exception = e
                -1
            end
            ok = n >= MIN_EXPOSURES
            println(log, join((group, r.number, r.name, n, ok), ','))
            flush(log)
            @printf("  %6d %-20s %4d exposures%s\n", r.number, r.name, n, ok ? "  taken" : "")
            ok || continue
            taken += 1
            push!(selected, (; group, r.number, r.name, r.period, r.U, r.amax, exposures=n))
        end
    end
    close(log)
    open(path * ".part", "w") do io
        println(io, "group,number,name,period,U,amax,exposures")
        foreach(s -> println(io, join((s.group, s.number, s.name, s.period, s.U, s.amax, s.exposures), ',')), selected)
    end
    mv(path * ".part", path)
    return read_table(path)
end

function cutout(dir, r)
    path = joinpath(dir, r["expid"] * ".fits")
    isfile(path) && return path
    ff = r["filefracday"]
    name = @sprintf("ztf_%s_%06d_%s_c%02d_%s_q%s_sciimg.fits", ff, parse(Int, r["field"]), r["filtercode"],
                    parse(Int, r["ccdid"]), r["imgtypecode"], r["qid"])
    url = "$IRSA/ibe/data/ztf/products/sci/$(ff[1:4])/$(ff[5:8])/$(ff[9:end])/$name" *
          "?center=$(r["ra_obj"]),$(r["dec_obj"])&size=600arcsec&gzip=false"
    body = nothing
    for attempt in 1:10
        try
            body = HTTP.get(url; retry=false, read_idle_timeout=120).body
            break
        catch e
            # IRSA throttles with 503; back off instead of giving the frame up
            if e isa HTTP.StatusError && e.status == 503 && attempt < 10
                sleep(min(2^attempt, 300))
            else
                @warn "cutout not fetched" r["expid"] exception = e
                return nothing
            end
        end
    end
    startswith(String(body[1:min(6, end)]), "SIMPLE") || return nothing
    write(path * ".part", body)
    mv(path * ".part", path)
    return path
end

function outcome(period, published)
    for (ratio, label) in ((1, "P"), (1 / 2, "P/2"), (2, "2P"))
        abs(period / (ratio * published) - 1) < TOLERANCE && return label
    end
    for ratio in (1, 1 / 2), k in (-2, -1, 1, 2)
        f = 1 / (ratio * published) + k / SIDEREAL_DAY_H
        f > 0 && abs(period * f - 1) < TOLERANCE && return "daily alias"
    end
    return "other"
end

function measure(s)
    number, published = parse(Int, s["number"]), parse(Float64, s["period"])
    dir = joinpath(DATA_DIR, s["number"])
    mkpath(dir)
    obs = best_window(exposures(number))
    # one at a time: IRSA answers parallel requests with HTTP 503
    paths = map(r -> cutout(dir, r), obs)
    have = paths .!== nothing
    obs, paths = obs[have], String.(paths[have])

    night(r) = floor(Int, parse(Float64, r["obsjd"]) - 0.25)
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
        foreach(m -> push!(rows, (; m.epoch, m.mag, m.mag_err, m.star_nearby)), lc)
    end
    t = [r for r in rows if isfinite(r.mag) && !r.star_nearby]
    # ZTF observes from Palomar (I41)
    reduced = reduced_magnitudes([r.epoch for r in t], [r.mag for r in t],
                                 object_geometry(number, [r.epoch for r in t]; observatory="I41"))
    open(joinpath(dir, "reduced.csv"), "w") do io
        println(io, "epoch,mag,mag_err")
        foreach(k -> @printf(io, "%.6f,%.4f,%.4f\n", reduced.epoch[k], reduced.mag[k], t[k].mag_err), eachindex(t))
    end
    search(mag) = fourier_period(reduced.epoch, mag; minimum_period=2 / 24, maximum_period=20 / 24, order=2)
    f = search(reduced.mag)
    control = search(shuffle(MersenneTwister(SEED + number), reduced.mag))
    return (; available=length(obs), used=length(t),
            span=maximum(reduced.epoch) - minimum(reduced.epoch), period=24f.period, period_error=24f.period_error,
            amplitude=f.amplitude, rms=f.rms, median_err=median(r.mag_err for r in t),
            result=outcome(24f.period, published), control=outcome(24control.period, published))
end

mkpath(joinpath(DATA_DIR, "most"))
sample = select_sample()
"--select" in ARGS && exit()

results = joinpath(DATA_DIR, "results.csv")
isfile(results) || write(results, "group,number,name,published,U,amax,available,used,span,period,period_error,amplitude,rms,median_err,result,control\n")
done = Set(r["number"] for r in read_table(results))
for s in sample
    s["number"] in done && continue
    head = (s["group"], s["number"], s["name"], s["period"], s["U"], s["amax"])
    line = try
        m = measure(s)
        @sprintf("%d,%d,%.2f,%.6f,%.6f,%.3f,%.4f,%.4f,%s,%s", m.available, m.used, m.span, m.period, m.period_error,
                 m.amplitude, m.rms, m.median_err, m.result, m.control)
    catch e
        @warn "not measured" s["number"] exception = e
        "0,0,0,NaN,NaN,NaN,NaN,NaN,not measured,not measured"
    end
    open(io -> println(io, join(head, ','), ',', line), results, "a")
    println(join(head, ' '), " -> ", line)
end
