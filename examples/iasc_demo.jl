#=
Runs the pipeline against real IASC (International Astronomical Search
Collaboration) practice campaign data — the one usage pattern documented
in the wiki's "Using real IASC campaign data" section but never actually
attempted, until now: local FITS files the user already has, not a
public URL this repo can fetch on its own (IASC campaign/practice
material isn't public the way the ZTF demo data is), so unlike
fetch_data.sh there is no download step here — place your own IASC
practice set(s) under data/real/iasc/<set name>/*.fits (one subdirectory
per set) before running this.

Data used to build and validate this script: 5 real Pan-STARRS1 (PS1)
practice sets (IASC's "Practice Image Sets", 2019-08-28/09-04/09-24),
each 4 exposures of the same field spanning ~40-70 min — no deep
reference stack is provided with these sets, so this only exercises the
raw (no-ZOGY) detection path, unlike real_data_demo.jl.

Three real bugs surfaced getting this far, all fixed and covered by
regression tests (see the Investigation Log):

- `load_wcs` raised "Linear transformation matrix is singular" on every
  one of these real headers — PS1's `CNPIX1`/`CNPIX2` legacy keywords
  make wcslib build a separate, degenerate implicit WCS description
  alongside the header's real, valid one. Fixed in `load_wcs` itself
  (see its docstring), since this could affect any survey using the same
  legacy convention, not just this dataset.
- PS1's pixel scale (~0.256"/px) is ~4x finer than ZTF's (~1.01"/px, used
  to calibrate real_data_demo.jl's `match_radius`/`max_speed`). Reusing
  those pixel-space defaults unchanged would silently apply a ~4x
  tighter angular tolerance than intended. Converted below from the same
  angular tolerances real_data_demo.jl uses, not copied as raw numbers.
- These FITS files mark invalid/masked pixels using the standard `BLANK`
  header keyword (scaled through `BZERO`/`BSCALE`) rather than `NaN` —
  `FITSIO.jl` does not convert these automatically. Left alone,
  `detect_sources` reads a masked region as enormous real flux: one real
  frame had 158,443 BLANK-sentinel pixels (2.7% of the image) and
  produced 176,165 spurious "detections" at `threshold=8.0`, which then
  made `run_pipeline`'s pairwise star-matching (quadratic in
  detections/frame) impractically slow — confirmed directly, not
  guessed: cleaning the same frame's BLANK pixels dropped it to 176. This
  is a general real-FITS-data concern, not specific to any one frame
  here, so it's handled as a preprocessing step below, before any frame
  reaches `run_pipeline`.
=#
using AsteroidPipeline
using FITSIO, Statistics

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "iasc")
const PS1_ARCSEC_PER_PIXEL = 0.2563   # measured from these sets' own CDELT1 (~7.118e-5 deg/px)
const ZTF_ARCSEC_PER_PIXEL = 1.01     # real_data_demo.jl's own field, for reference

# First attempt reused real_data_demo.jl's match_radius=10.0 px at ZTF's
# pixel scale, converted to keep the same ~10.1" angular tolerance — a
# real, measured mistake: PS1's own headers report each frame's actual
# astrometric solution quality directly (PERROR, the per-star positional
# RMS residual — 0.20-0.23" across the fields checked here, not
# guessed), ~50x tighter than that 10" tolerance. On the densest field
# tested, the 10" version produced 10,422 tracklets from ~500
# detections/frame — almost certainly distinct real stars within 10" of
# each other getting cross-linked as spurious tracklets. 2" (~10x
# PERROR, comfortable margin for real motion and centroiding noise, not
# just the bare residual) is measured directly to still recover every
# known object the looser radius did, at a fraction of the tracklet
# count — see the Investigation Log for the before/after.
const MATCH_RADIUS = 2.0 / PS1_ARCSEC_PER_PIXEL
const MAX_SPEED = 5000.0 * ZTF_ARCSEC_PER_PIXEL / PS1_ARCSEC_PER_PIXEL

isdir(DATA_DIR) || error("no $DATA_DIR — place your own IASC practice FITS sets there first " *
                          "(one subdirectory per set, e.g. data/real/iasc/XY25_p10/*.fits)")
sets = sort(filter(isdir, readdir(DATA_DIR, join=true)))
isempty(sets) && error("no set subdirectories found in $DATA_DIR")

const CLEANED_DIR = joinpath(@__DIR__, "..", "data", "real", "iasc_cleaned")

"""
    clean_blank_pixels(path) -> String

Write a copy of the FITS file at `path` with `BLANK`-sentinel pixels (see
the header comment above) replaced by the frame's own valid-region
median, to `CLEANED_DIR`, preserving the original header exactly (so WCS
and `MJD-OBS` are untouched) — returns the cleaned copy's path, reusing
an existing one if already written. A no-op (returns `path` unchanged)
if the header has no `BLANK` keyword.
"""
function clean_blank_pixels(path)
    out = joinpath(CLEANED_DIR, basename(path))
    isfile(out) && return out
    cleaned = FITS(path, "r") do fin
        hdu = fin[1]
        h = read_header(hdu)
        haskey(h, "BLANK") || return false
        bzero = haskey(h, "BZERO") ? h["BZERO"] : 0.0
        bscale = haskey(h, "BSCALE") ? h["BSCALE"] : 1.0
        sentinel = h["BLANK"] * bscale + bzero
        img = Float64.(read(hdu))
        is_blank = img .== sentinel
        any(is_blank) || return false
        img[is_blank] .= median(img[.!is_blank])
        mkpath(CLEANED_DIR)
        FITS(out, "w") do fout
            write(fout, img; header=h)
        end
        return true
    end
    return cleaned ? out : path
end

function summarize(label, candidates)
    n_tracklets = length(unique(candidates.id))
    println("-- $label: $n_tracklets tracklet(s) from $(length(candidates)) detections --")
    n_tracklets == 0 && return

    first_rows = [first(filter(r -> r.id == id, candidates)) for id in unique(candidates.id)]
    matches = crossmatch_catalog(first_rows, :skybot; radius=15.0)
    known_ids = Set(matches.id)
    for m in matches
        println("  id=$(m.id): $(m.name) ($(m.class), Mv=$(m.mv)), offset $(round(m.distance_arcsec, digits=1))\"")
    end
    println("  $(length(known_ids))/$(n_tracklets) match a known SkyBoT object; ",
            "$(n_tracklets - length(known_ids)) unmatched (candidates for human vetting).")
end

for set_dir in sets
    raw_paths = sort(filter(p -> endswith(p, ".fits"), readdir(set_dir, join=true)))
    isempty(raw_paths) && continue
    println("\n=== $(basename(set_dir)): $(length(raw_paths)) frames ===")
    paths = clean_blank_pixels.(raw_paths)
    # threshold=8.0 (run_pipeline's own default is 5.0), matching
    # real_data_demo.jl's ZTF choice — a reasonable, non-arbitrary
    # starting point for a new survey, not re-derived from scratch here.
    candidates = run_pipeline(paths; threshold=8.0, match_radius=MATCH_RADIUS, max_speed=MAX_SPEED)
    summarize(basename(set_dir), candidates)
end
