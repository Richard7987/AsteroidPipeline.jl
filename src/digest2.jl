"""
    digest2_score(candidates, station::AbstractString;
                  digest2_path::AbstractString="digest2",
                  config_dir::Union{Nothing,AbstractString}=nothing,
                  trksub_prefix::AbstractString="") -> Table

Score each tracklet in `candidates` (an [`astrometric_calibrate`](@ref)
table, as [`mpc80_report`](@ref) takes) for its likelihood of belonging to
each orbit class, with the Minor Planet Center's `digest2` (short-arc
statistical ranging against a solar-system population model).

`digest2` is a separate C program
(`https://github.com/Smithsonian/mpc-public`, directory `digest2`): build it
(see [Follow-up tools](@ref "Orbit class: digest2")) and put the executable on
`PATH` or pass it as `digest2_path`. `config_dir` defaults to the directory
of the executable, symlinks resolved, and must hold `digest2.model.csv`, as it ships with `digest2`
(`digest2.obscodes` too, which `digest2` downloads from the MPC if missing).

It runs as the MPC runs it: with the astrometric uncertainty of each station
from `MPC.config` in `config_dir` (1″ for a station not listed, or with no
`MPC.config`), and seeded the same way for every tracklet, so a tracklet
always gets the same scores.

The input is [`mpc80_report`](@ref)'s output, sorted by `(id, epoch)` as
`digest2` requires. The command is built from an argument vector, never a
shell string.

Returns one row per tracklet (`id`) with `digest2`'s `RMS` (its
linear-motion fit residual, arcsec; `NaN` if too large to print) and four
0-100 pseudo-probabilities: `int_score` (any orbit class of MPC interest),
`neo_score` (near-Earth object), and `n22_score`/`n18_score` (NEO with H ≤
22 / ≤ 18). The MPC posts a tracklet with `neo_score` of 65 or more to its
NEO Confirmation Page.

These scores classify orbits; they don't tell real objects from bogus
ones. A stationary star linked into a tracklet by a loose `match_radius`
scores `neo_score=100`, so vet tracklets for consistent motion first.

Throws `ArgumentError` if `digest2_path` resolves to nothing runnable or
`station` isn't a 3-character MPC observatory code.
"""
function digest2_score(candidates, station::AbstractString;
                        digest2_path::AbstractString="digest2",
                        config_dir::Union{Nothing,AbstractString}=nothing,
                        trksub_prefix::AbstractString="")
    length(station) == 3 || throw(ArgumentError("station must be a 3-character MPC observatory code"))

    resolved = isfile(digest2_path) ? digest2_path : Sys.which(digest2_path)
    resolved === nothing && throw(ArgumentError(
        "digest2 executable not found at or on PATH as \"$digest2_path\" — build it " *
        "from https://github.com/Smithsonian/mpc-public and pass its path via digest2_path"))
    dir = config_dir === nothing ? dirname(realpath(resolved)) : config_dir

    sorted = sort(candidates; by=row -> (row.id, row.epoch))
    input = mpc80_report(sorted, station; trksub_prefix=trksub_prefix)

    mpc_config = joinpath(dir, "MPC.config")
    obserr = isfile(mpc_config) ? filter(l -> startswith(strip(l), "obserr"), readlines(mpc_config)) : String[]
    out = IOBuffer()
    mktempdir() do tmp
        config = joinpath(tmp, "digest2.config")
        write(config, join(["noheadings", "rms", "noid", "repeatable", "Int", "NEO", "N22", "N18", obserr...], "\n") * "\n")
        run(pipeline(`$resolved -p $dir -c $config -`, stdin=IOBuffer(input), stdout=out))
    end

    ids, rmss, ints, neos, n22s, n18s = Int[], Float64[], Float64[], Float64[], Float64[], Float64[]
    for line in split(String(take!(out)), "\n"; keepempty=false)
        tokens = split(line)
        length(tokens) == 6 || continue
        designation = tokens[1][(length(trksub_prefix)+1):end]
        push!(ids, parse(Int, designation; base=36))
        push!(rmss, something(tryparse(Float64, tokens[2]), NaN))
        push!(ints, parse(Float64, tokens[3]))
        push!(neos, parse(Float64, tokens[4]))
        push!(n22s, parse(Float64, tokens[5]))
        push!(n18s, parse(Float64, tokens[6]))
    end

    return Table(id=ids, rms=rmss, int_score=ints, neo_score=neos, n22_score=n22s, n18_score=n18s)
end
