# Analyze image sets

An *image set* is three or more exposures of one field, taken over tens of
minutes to a few hours, in time order or not. [`analyze_set`](@ref) takes it
as any of:

| Given as | Example |
|:--|:--|
| a `.zip` of FITS files | `analyze_set("set.zip"; …)` — unpacked next to the `.zip` |
| a folder of FITS files | `analyze_set("images/"; …)` |
| a list of FITS files | `analyze_set(["a.fits", "b.fits", "c.fits"]; …)` |

## One set

```julia
using AsteroidPipeline

orbits = read_mpcorb(download_mpcorb())        # the MPC's orbits, refreshed daily
r = analyze_set("set.zip"; observatory="F51", workdir="work", mpcorb=orbits)
```

Reuse `orbits` across sets: reading the file (1.5 million orbits) takes ~10 s.

For an IASC campaign, pass also `itf=download_itf()`: candidates Pan-STARRS
already reported are marked, as IASC would reject them as duplicates. The
file is ~130 MB; reading one night of it takes ~11 s.

### What comes back

| Field | Meaning |
|:--|:--|
| `name`, `paths`, `epochs` | the set, its frames, their mid-exposure times (Julian Date, UTC) |
| `wcs` | each frame's astrometric solution, refitted to Gaia DR3 |
| `ratios`, `frames_used`, `dropped`, `reliable` | frame quality: 10σ sources per reference star; frames above 30 are left out; a set with fewer than 3 good frames is not searched |
| `candidates` | every measured point: tracklet `id`, `frame`, `ra`, `dec`, `epoch`, `mag` |
| `known` | tracklet id ⇒ known object it was identified as |
| `unknown_ids` | the rest: candidates for a person to check |
| `nearby` | for each candidate, the known objects within 10′ and whether one is probably it — see [Is a candidate new?](../concepts/known-objects.md) |
| `neocp` | for each candidate, the objects on the MPC's NEO Confirmation Page near it, for images up to 30 days old |
| `itf` | with `itf=download_itf()`, the candidates already reported to the MPC unlinked, and by whom — see [Already reported](../validation/itf.md) |
| `digest2` | with `digest2=<the executable>`, each candidate's digest2 scores, as the MPC computes them; 65 or more in `neo_score` is its threshold for the NEO Confirmation Page — see [Near-Earth or not](../validation/digest2.md) |
| `expected` | every known object predicted on the frames, detected or not |
| `report`, `report_kind` | the MPC report: `:blank`, `:candidates` or `:unreliable` |
| `warnings` | services that couldn't be reached, frames that couldn't be calibrated |

Progress can be followed with a callback:

```julia
analyze_set("set.zip"; observatory="F51", workdir="work",
            progress=(stage, fraction) -> println(stage, " ", round(Int, 100fraction), "%"))
```

## Many sets

`examples/analyze_sets.jl` analyzes every `.zip` (or every folder of FITS
images) in a folder and prints one line per set:

```sh
julia --project=. examples/analyze_sets.jl ~/Downloads/sets F51
```

```
set                      result                 candidates  known  frames
XY45_p01                 no new objects                  0      2       4
XY45_p10                 no new objects                  0      2       4
XY45_p11                 no new objects                  0      1       4
```

and writes each set's report next to it as `<set>_MPCReport.txt`.

## Without the network, and the same again

Every answer of Gaia DR3, SkyBoT and JPL's Small-Body Identification service
is kept on disk ([`gaia_cache_dir`](@ref), [`answers_dir`](@ref)). A set
analyzed once with the network can be analyzed again without it:

```julia
orbits = read_mpcorb(default_mpcorb_path())  # the copy download_mpcorb kept, not a new download
r = analyze_set("XY45_p01.zip"; workdir="work", observatory="F51", mpcorb=orbits, offline=true)
```

No service is asked, and the result is the one the answers kept give, even if
a service has since changed its answers or is down. A query never answered
before leaves a warning (known objects) or an error (Gaia DR3, which the
astrometry needs); the NEO Confirmation Page, which lists only what is being
observed now, is not checked.

## Lower level

[`analyze_set`](@ref) is built on [`run_pipeline`](@ref), which takes FITS
paths and keyword arguments for every stage and returns the calibrated
tracklets only. Use it to change a stage's settings, or to run without the
network:

```julia
candidates = run_pipeline(paths; pipeline_config(0.25)...)   # 0.25"/pixel
```

[`pipeline_config`](@ref) gives the settings validated on real data, scaled
to a pixel scale; [How it works](../concepts/pipeline.md) explains each.
