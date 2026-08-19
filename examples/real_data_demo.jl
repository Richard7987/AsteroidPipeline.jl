#=
Runs the pipeline against real ZTF (Zwicky Transient Facility) public
survey data twice — once on raw science frames, once with ZOGY difference
imaging against a deep reference stack (see src/zogy.jl, src/reference.jl)
— and cross-matches both against SkyBoT, to measure whether differencing
actually helps rather than assume it does.

Data: field 451, CCD 1, quadrant 1, zr filter. Science: 2019-10-23, 5
exposures spread over the night's 6.25 h span (seeing 1.8-2.0 px, maglim
20.1-20.5). Reference: 30 frames from other nights, best seeing first.
Not included in the repo; fetch first with `fetch_data.sh` in this
directory. See the plan this was built from (or git log) for how this
field/night/reference set was chosen, including a dead end: an earlier
choice of night turned out to have no catalogued object above that
night's own detection limit, which would have made ZOGY look like it
failed when the real problem was upstream.

Ground truth: SkyBoT reports 135992 "2002 UY45" (MB>Middle, Mv 18.4) in
this field, 1.9 mag above the night's own 5-sigma limit and thus
detectable even without differencing — this is what "recovered by both"
below should include. It moves ~222" (~219 px) over the full baseline.

Building the 30-frame reference takes a while (each frame is individually
reprojected onto the science grid — roughly half an hour on a laptop,
sequentially); this only has to happen once, not once per science frame.
Spread across worker processes below via build_reference's `workers`
keyword — a real, measured 3.21x on this exact dataset (951s -> 296s,
8 worker processes; see the Investigation Log's Design refinements page).
=#
using AsteroidPipeline
using Distributed

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real")
const SCIENCE_PATHS = joinpath.(DATA_DIR, "science", [
    "ztf_20191023195590_000451_zr_c01_o_q1_sciimg.fits",
    "ztf_20191023267002_000451_zr_c01_o_q1_sciimg.fits",
    "ztf_20191023347164_000451_zr_c01_o_q1_sciimg.fits",
    "ztf_20191023425208_000451_zr_c01_o_q1_sciimg.fits",
    "ztf_20191023455822_000451_zr_c01_o_q1_sciimg.fits",
])

for path in SCIENCE_PATHS
    isfile(path) || error("missing $path — run examples/fetch_data.sh first")
end
reference_paths = readdir(joinpath(DATA_DIR, "reference"), join=true)
isempty(reference_paths) && error("no reference frames found — run examples/fetch_data.sh first")

function summarize(label, candidates)
    n_tracklets = length(unique(candidates.id))
    println("\n-- $label: $n_tracklets tracklet(s) from $(length(candidates)) detections --")
    n_tracklets == 0 && return Set{Int}()

    first_rows = [first(filter(r -> r.id == id, candidates)) for id in unique(candidates.id)]
    matches = crossmatch_catalog(first_rows, :skybot; radius=15.0)
    known_ids = Set(matches.id)
    for m in matches
        println("  id=$(m.id): $(m.name) ($(m.class), Mv=$(m.mv)), offset $(round(m.distance_arcsec, digits=1))\"")
    end
    println("  $(length(known_ids))/$(n_tracklets) match a known SkyBoT object; ",
            "$(n_tracklets - length(known_ids)) unmatched (candidates for human vetting).")
    return known_ids
end

println("Baseline: detection directly on science frames (no differencing)...")
baseline = run_pipeline(SCIENCE_PATHS; timestamp_key="OBSMJD", threshold=8.0,
                         match_radius=10.0, max_speed=5000.0)
baseline_known = summarize("baseline", baseline)

println("\nBuilding a deep reference from $(length(reference_paths)) frames (this is the slow part)...")
sci1 = load_frame(SCIENCE_PATHS[1])
refs = [load_frame(p) for p in reference_paths]
reference_workers = addprocs(Sys.CPU_THREADS; exeflags="--project=$(Base.active_project())")
@everywhere reference_workers using AsteroidPipeline
ref_image, ref_sigma, ref_mask = try
    build_reference(refs, sci1.wcs, size(sci1.image); workers=reference_workers)
finally
    rmprocs(reference_workers)
end
ref_psf = estimate_psf(ref_image; threshold=15.0)
reference = (image=ref_image, sigma=ref_sigma, mask=ref_mask, psf=ref_psf, wcs=sci1.wcs)
println("reference built: sigma=$(round(ref_sigma, digits=3)), ",
        "coverage=$(round(100 * count(ref_mask) / length(ref_mask), digits=1))%")

println("\nZOGY: detection on the difference image against that reference...")
# min_frames' default now automatically accounts for whatever frames the
# per-frame quality gate (quality_max_std, default 1.5) excludes — real
# on this dataset, see INVESTIGATION_LOG.md — no manual override needed.
zogy_result = run_pipeline(SCIENCE_PATHS; timestamp_key="OBSMJD", threshold=6.0,
                            match_radius=10.0, max_speed=5000.0, reference=reference)
zogy_known = summarize("ZOGY", zogy_result)

println("\n== Summary ==")
println("baseline: $(length(unique(baseline.id))) tracklets, $(length(baseline_known)) known")
println("ZOGY:     $(length(unique(zogy_result.id))) tracklets, $(length(zogy_known)) known")
if isempty(baseline_known) && !isempty(zogy_known)
    println("ZOGY recovered a real object the undifferenced baseline missed.")
elseif length(zogy_known) <= length(baseline_known)
    println("ZOGY did not recover more known objects than the baseline here — ",
            "worth investigating before trusting this path on fainter fields.")
end
