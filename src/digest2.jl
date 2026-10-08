"""
    digest2_score(candidates, station::AbstractString;
                  digest2_path::AbstractString="digest2",
                  config_dir::Union{Nothing,AbstractString}=nothing,
                  trksub_prefix::AbstractString="") -> Table

Score each tracklet in `candidates` (an [`astrometric_calibrate`](@ref)
table, same shape [`mpc80_report`](@ref) takes) for its likelihood of
belonging to several real orbit classes, via the Minor Planet Center's
own `digest2` — a real, external, already-validated classifier (short-arc
statistical ranging against a real solar-system population model), not
anything trained or guessed at here.

`digest2` is a separate C program (source:
`https://github.com/Smithsonian/digest2`), not a Julia dependency —
build it yourself (`make` in its `digest2/` directory; see
[MPC digest2 Scoring](@ref) for exact steps) and either put the
resulting `digest2` executable on `PATH` or pass its path via
`digest2_path`. `config_dir` defaults to the directory of the *real*
executable — symlinks resolved, so a link on `PATH` (say `~/.local/bin`)
to a build tree still finds that tree's files: `digest2`'s own documented
"same directory as the executable" layout — and must contain
`digest2.model.csv`, `digest2.obscodes`, and an `MPC.config` — all three
ship with the `digest2` repository itself and are used unmodified here
(this pipeline supplies none of its own: `digest2`'s upstream `MPC.config`
already sets a real, tuned per-observatory error, including ZTF's own
`I41`, and rewriting it would just be guessing at numbers `digest2`'s own
maintainers already calibrated).

Reuses [`mpc80_report`](@ref) verbatim to build `digest2`'s input (`digest2`
requires MPC 80-column input grouped by consecutive same-designation
lines, sorted by designation then time — `candidates` is sorted by
`(id, epoch)` locally first, since [`astrometric_calibrate`](@ref) only
guarantees grouping by `id`, not epoch order across the whole table, and
a caller's own row order should not silently change `digest2`'s result).
Piped via `stdin`/`stdout` through a `Cmd` built from an argument vector
(never a shell string), so `digest2_path`/`config_dir` cannot be
mis-parsed as shell syntax.

Returns one row per tracklet (`id`, matching `candidates`' own `id`) with
`digest2`'s `RMS` (its own linear-motion fit residual, arcsec — a
free byproduct neither [`link_candidates`](@ref) nor
[`astrometric_calibrate`](@ref) compute) and four raw 0-100
pseudo-probability scores: `int_score` (any orbit class of general MPC
interest), `neo_score` (near-Earth object), `n22_score`/`n18_score` (NEO
with absolute magnitude ≤ 22/≤ 18 respectively — i.e. how large the
object would have to be). These are `digest2`'s own default output
columns, present regardless of `MPC.config`'s requested class list (real
behavior, confirmed by actually building and running `digest2`, not
assumed from its docs) — this function does not parse the trailing
"Other Possibilities" free-text column (main-belt/Trojan/comet class
scores for whatever doesn't make the main four).

**Real finding, not a design assumption — `digest2`'s scores are an
orbit-class classifier, not a real/bogus one, and this matters in
practice**: run for real against `real_data_demo.jl`'s ZTF field 451
baseline tracklets (see [MPC digest2 Scoring](@ref) for the full
numbers), the two real, SkyBoT-confirmed known objects scored `neo_score`
5 and 7 (correctly low — they're Main Belt, not NEOs), while 131 of the
other 133 tracklets scored `neo_score=100`. Those 131 were not missed
NEOs: inspecting one directly, its 5 points jitter by under 1" across
the full 6.25 h baseline (`~3"/day` implied rate — this is a real,
stationary star, re-detected each frame and "linked" into a bogus
tracklet only because [`link_candidates`](@ref)'s `match_radius` here
(10", well above ZTF's real sub-arcsec centroiding precision) is loose
enough to also accept a star's own centroid jitter as if it were
motion). `digest2` has no way to know that near-zero, jittery apparent
motion is a stationary star rather than a real object at a range where
that motion is dynamically self-consistent — it scored the *hypothesis*
correctly, on input that was never a real tracklet to begin with.

**The actionable lesson**: don't sort by `neo_score` on raw
[`link_candidates`](@ref) output as a "most interesting first" filter —
tighten `match_radius` to the survey's own real astrometric precision
first (`iasc-campaign-validation.md` already established this exact
principle from PS1's `PERROR`), or otherwise vet tracklets for genuine,
consistent motion, before scoring; `digest2_score` is a real, valuable
tool for classifying the dynamical class of an *already-plausible*
tracklet, not a substitute for that upstream quality control.

Throws `ArgumentError` if `digest2_path` resolves to nothing runnable
(`Sys.which`-style — an explicit argument, not an environment variable
this function reads itself, matching [`plate_solve`](@ref)'s own
`api_key` convention) or if `station` isn't a 3-character MPC
observatory code.
"""
function digest2_score(candidates, station::AbstractString;
                        digest2_path::AbstractString="digest2",
                        config_dir::Union{Nothing,AbstractString}=nothing,
                        trksub_prefix::AbstractString="")
    length(station) == 3 || throw(ArgumentError("station must be a 3-character MPC observatory code"))

    resolved = isfile(digest2_path) ? digest2_path : Sys.which(digest2_path)
    resolved === nothing && throw(ArgumentError(
        "digest2 executable not found at or on PATH as \"$digest2_path\" — build it " *
        "from https://github.com/Smithsonian/digest2 and pass its path via digest2_path"))
    dir = config_dir === nothing ? dirname(realpath(resolved)) : config_dir

    sorted = sort(candidates; by=row -> (row.id, row.epoch))
    input = mpc80_report(sorted, station; trksub_prefix=trksub_prefix)

    out = IOBuffer()
    run(pipeline(`$resolved -p $dir -`, stdin=IOBuffer(input), stdout=out))
    lines = split(chomp(String(take!(out))), "\n")

    ids, rmss, ints, neos, n22s, n18s = Int[], Float64[], Float64[], Float64[], Float64[], Float64[]
    for line in lines[2:end]  # skip the "Desig. RMS Int NEO N22 N18 ..." header
        tokens = split(line)
        designation = tokens[1][(length(trksub_prefix)+1):end]
        push!(ids, parse(Int, designation; base=36))
        push!(rmss, parse(Float64, tokens[2]))
        push!(ints, parse(Float64, tokens[3]))
        push!(neos, parse(Float64, tokens[4]))
        push!(n22s, parse(Float64, tokens[5]))
        push!(n18s, parse(Float64, tokens[6]))
    end

    return Table(id=ids, rms=rmss, int_score=ints, neo_score=neos, n22_score=n22s, n18_score=n18s)
end
