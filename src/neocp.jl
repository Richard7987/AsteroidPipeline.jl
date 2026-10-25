#=
Is a candidate already on the MPC's NEO Confirmation Page? The NEOCP lists
objects reported in the last days or weeks that have no orbit yet — so no
catalogue (SkyBoT, MPCORB) knows them, and a candidate that is one of them
would look new. Recent images only: the NEOCP lists what is being observed
now, and its ephemerides, from arcs of hours, are only good near the
observations they come from.
=#

const NEOCP_LIST_URL = "https://www.minorplanetcenter.net/Extended_Files/neocp.json"
const NEOCP_EPHEMERIS_URL = "https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi"

"""
    neocp_objects() -> Vector{NamedTuple}

The objects on the MPC's [NEO Confirmation Page](https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html)
now (a live request): `designation` (temporary), `score` (the MPC's 0-100
likelihood of being a NEO), `ra`, `dec` (degrees, at the last observation),
`mag` (V), `n_obs`, `arc_days` and `not_seen_days`.
"""
function neocp_objects()
    body = String(HTTP.get(NEOCP_LIST_URL; read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body)
    return [(designation=String(o["Temp_Desig"]), score=Float64(o["Score"]), ra=15 * Float64(o["R.A."]),
             dec=Float64(o["Decl."]), mag=Float64(something(o["V"], NaN)), n_obs=Int(o["NObs"]),
             arc_days=Float64(o["Arc"]), not_seen_days=Float64(o["Not_Seen_dys"]))
            for o in JSON.parse(body)]
end

const _NEOCP_ROW = r"^(\d{4}) (\d\d) (\d\d) (\d\d)(\d\d)\s+(\d\d) (\d\d) (\d\d\.\d)\s+([+-])(\d\d) (\d\d) (\d\d)\s+\S+\s+\S+\s+([\d.]+)\s+([\d.]+)"

# Ephemerides of `designations` from `observatory`, every 30 minutes for ten
# hours from the hour before `first_jd` (UTC): designation => rows of `jd`,
# `ra`, `dec`, `rate_arcsec_h`, `pa_deg`, `uncertain` (the MPC's "!" marks).
function _neocp_ephemerides(designations, observatory::AbstractString, first_jd::Real;
                            now_jd::Real=datetime2julian(Dates.now(Dates.UTC)))
    start = floor(Int, (first_jd - now_jd) * 24) - 1          # hours from now
    form = ["mb=-30", "mf=30", "dl=-90", "du=%2B90", "nl=0", "nu=100", "sort=d", "W=j",
            ("obj=" * HTTP.escapeuri(d) for d in designations)...,
            "Parallax=1", "obscode=$observatory", "long=", "lat=", "alt=", "int=1", "start=$start",
            "raty=a", "mot=m", "dmot=p", "out=f", "sun=x", "oalt=-90"]
    html = String(HTTP.post(NEOCP_EPHEMERIS_URL, ["Content-Type" => "application/x-www-form-urlencoded"],
                            join(form, "&"); read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body)
    return _parse_neocp_ephemerides(html, designations)
end

# The ephemeris page the MPC answers with -> designation => rows.
function _parse_neocp_ephemerides(html::AbstractString, designations)
    out = Dict{String,Vector{NamedTuple}}()
    current = nothing
    for raw in split(replace(html, r"<[^>]*>" => ""), '\n')
        line = rstrip(raw)
        if (m = match(r"^\s*([A-Za-z0-9]{7})\s*$", line)) !== nothing && m[1] in designations
            current = String(m[1])
            out[current] = NamedTuple[]
            continue
        end
        current === nothing && continue
        m = match(_NEOCP_ROW, line)
        m === nothing && continue
        jd = datetime2julian(DateTime(parse.(Int, (m[1], m[2], m[3], m[4], m[5]))...))
        ra = 15 * (parse(Float64, m[6]) + parse(Float64, m[7]) / 60 + parse(Float64, m[8]) / 3600)
        dec = (m[9] == "-" ? -1 : 1) * (parse(Float64, m[10]) + parse(Float64, m[11]) / 60 + parse(Float64, m[12]) / 3600)
        push!(out[current], (jd=jd, ra=ra, dec=dec, rate_arcsec_h=60 * parse(Float64, m[13]),
                             pa_deg=parse(Float64, m[14]), uncertain=occursin('!', line)))
    end
    return out
end

# Linear interpolation of an ephemeris at `jd`; `nothing` outside it.
function _interpolate_ephemeris(rows, jd)
    i = findlast(r -> r.jd <= jd, rows)
    (i === nothing || i == length(rows)) && return nothing
    a, b = rows[i], rows[i+1]
    f = (jd - a.jd) / (b.jd - a.jd)
    dra = mod(b.ra - a.ra + 180, 360) - 180
    return (ra=mod(a.ra + f * dra, 360), dec=a.dec + f * (b.dec - a.dec), epoch=jd)
end

"""
    neocp_check(candidates; observatory, radius_arcmin=30, max_age_days=30,
                now_jd=<now>) -> (checks, warnings)

For every tracklet in `candidates`, the objects of the MPC's NEO
Confirmation Page predicted within `radius_arcmin` of it at its own epochs,
from `observatory` (the MPC's own ephemerides, a live request). `checks`
maps each tracklet id to `(verdict, objects)`: `objects` are `designation`,
`score`, `separation_arcsec` (median over the frames), the object's predicted
`rate_arcsec_h` and `pa_deg`, `uncertain` (the MPC marks the ephemeris as
uncertain) and `agreement`: `:motion` if the object moves like the tracklet
(rate within 30%, direction within 15°), else `:none`. `verdict` is
`:on_neocp` if one agrees, else `:none`.

A candidate on the NEOCP is someone else's recent, unconfirmed find: report
it under that designation ([`mpc_report`](@ref)'s `designations`) — the
observations help confirm it. The radius is wide because NEOCP ephemerides,
from arcs of hours, are uncertain by arcminutes.

Images older than `max_age_days` are not checked (a warning says so): the
NEOCP lists only what is being observed now.
"""
function neocp_check(candidates; observatory::AbstractString, radius_arcmin::Real=30, max_age_days::Real=30,
                     now_jd::Real=datetime2julian(Dates.now(Dates.UTC)))
    rows = collect(candidates)
    warnings = String[]
    ObjT = NamedTuple{(:designation, :score, :separation_arcsec, :rate_arcsec_h, :pa_deg, :uncertain, :agreement),
                      Tuple{String,Float64,Float64,Float64,Float64,Bool,Symbol}}
    checks = Dict{Int,Any}()
    isempty(rows) && return (checks=checks, warnings=warnings)
    epochs = [r.epoch for r in rows]
    age = now_jd - maximum(epochs)
    if age > max_age_days
        push!(warnings, "NEO Confirmation Page not checked: the images are $(round(Int, age)) days old, " *
                        "and the NEOCP lists only objects being observed now")
        return (checks=checks, warnings=warnings)
    end
    tracklets = Dict(id => sort(filter(r -> r.id == id, rows); by=r -> r.epoch) for id in unique(r.id for r in rows))
    try
        objects = neocp_objects()
        # objects listed within 45° of a candidate (they move fast, but not that fast in a few days)
        near = filter(o -> any(t -> _angular_distance_arcsec(o.ra, o.dec, t[1].ra, t[1].dec) < 45 * 3600,
                               values(tracklets)), objects)
        ephemerides = isempty(near) ? Dict{String,Vector{NamedTuple}}() :
                      _neocp_ephemerides([o.designation for o in near], observatory, minimum(epochs); now_jd)
        scores = Dict(o.designation => o.score for o in near)
        for (id, pts) in tracklets
            motion = tracklet_motion(pts)
            found = ObjT[]
            for (desig, eph) in ephemerides
                pred = [_interpolate_ephemeris(eph, p.epoch) for p in pts]
                any(isnothing, pred) && continue
                sep = median(_angular_distance_arcsec(p.ra, p.dec, q.ra, q.dec) for (p, q) in zip(pts, pred))
                sep <= radius_arcmin * 60 || continue
                m = tracklet_motion(pred)
                same = abs(m.rate_arcsec_h - motion.rate_arcsec_h) <= 0.3 * max(motion.rate_arcsec_h, 1.0) &&
                       abs(mod(m.pa_deg - motion.pa_deg + 180, 360) - 180) <= 15
                unc = any(r -> r.uncertain && pts[1].epoch - 0.05 <= r.jd <= pts[end].epoch + 0.05, eph)
                push!(found, (designation=desig, score=get(scores, desig, NaN), separation_arcsec=sep,
                              rate_arcsec_h=m.rate_arcsec_h, pa_deg=m.pa_deg, uncertain=unc,
                              agreement=same ? :motion : :none))
            end
            sort!(found; by=o -> (o.agreement !== :motion, o.separation_arcsec))
            checks[id] = (verdict=any(o -> o.agreement === :motion, found) ? :on_neocp : :none, objects=found)
        end
    catch e
        (_is_unavailable(e) || e isa ArgumentError) || rethrow()
        push!(warnings, "NEO Confirmation Page unavailable: $(sprint(showerror, e))")
    end
    return (checks=checks, warnings=warnings)
end
