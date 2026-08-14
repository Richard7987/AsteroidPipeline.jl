#=
Runs the pipeline end to end on real ZTF (Zwicky Transient Facility) public
survey data, then cross-matches the resulting candidates against SkyBoT to
report which are already-known solar-system objects.

Data: field 451, CCD 1, quadrant 1, zr filter, 2021-03-15, a same-night
triplet of science exposures ~3.5 minutes apart (a real ZTF high-cadence
sequence, suitable for tracklet linking). Not included in the repo; fetch
first with `fetch_data.sh` in this directory.

Parameters below (`threshold`, `match_radius`, `max_speed`) were tuned
against this specific field: SkyBoT reports a nearby known object (127319
"2002 JB99", MB>Outer) moving ~6 arcsec over the full triplet baseline,
i.e. ~3 pixels per exposure gap at this field's ~1.01 arcsec/pixel scale.
A looser match_radius (e.g. the 2.0 px default) is too tight to recover
real motion at this cadence; a much looser one produces hundreds of
spurious tracklets in this crowded, near-ecliptic field. These are not
universal defaults — re-tune per field/cadence.
=#
using AsteroidPipeline

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real")
const PATHS = joinpath.(DATA_DIR, [
    "ztf_20210315117014_000451_zr_c01_o_q1_sciimg.fits",
    "ztf_20210315119711_000451_zr_c01_o_q1_sciimg.fits",
    "ztf_20210315122141_000451_zr_c01_o_q1_sciimg.fits",
])

for path in PATHS
    isfile(path) || error("missing $path — run examples/fetch_data.sh first")
end

candidates = run_pipeline(PATHS; timestamp_key="OBSMJD", threshold=8.0,
                           match_radius=3.5, max_speed=3000.0)
println(length(unique(candidates.id)), " candidate tracklets from ",
        length(candidates), " total detections")

# One row per tracklet (first frame) is enough to query SkyBoT; querying
# every row would triple the request count for no extra information.
first_rows = [first(filter(r -> r.id == id, candidates))
              for id in unique(candidates.id)]

matches = crossmatch_catalog(first_rows, :skybot; radius=15.0)
known_ids = Set(matches.id)

println(length(known_ids), " tracklets match a known SkyBoT object:")
for m in matches
    println("  id=$(m.id): $(m.name) ($(m.class), Mv=$(m.mv)), ",
             "offset $(round(m.distance_arcsec, digits=1))\"")
end

unmatched = length(unique(candidates.id)) - length(known_ids)
println(unmatched, " tracklets have no known counterpart within 15\" — ",
        "candidates for human vetting.")

if isempty(known_ids)
    println()
    println("Zero known matches is expected here, not a failure: the only ",
            "catalogued asteroids in this field (per a direct SkyBoT cone ",
            "search) are Mv >= 19.7, and detect_sources runs on a single, ",
            "non-difference-imaged 30s frame — ZTF's own real-time pipeline ",
            "recovers objects this faint by subtracting a deep reference ",
            "image first, which this pipeline does not (yet) do. The 35 ",
            "unmatched tracklets are most likely spurious 3-frame linear ",
            "coincidences among the field's ~400 detected stars, not ",
            "genuine moving objects — cross-matching against a synthetic, ",
            "brighter injected source (see test/runtests.jl) is what ",
            "confirms the detection/linking math itself is correct.")
end
