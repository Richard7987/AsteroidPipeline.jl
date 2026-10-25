"""
    digest2_score(candidates, station::AbstractString;
                  digest2_path::AbstractString="digest2",
                  config_dir::Union{Nothing,AbstractString}=nothing,
                  trksub_prefix::AbstractString="") -> Table

Score each tracklet in `candidates` (an [`astrometric_calibrate`](@ref)
table, as [`mpc80_report`](@ref) takes) for its likelihood of belonging to
each orbit class, with the Minor Planet Center's `digest2` (short-arc
statistical ranging against a solar-system population model).

`digest2` is a separate C program (`https://github.com/Smithsonian/digest2`):
build it (`make` in its `digest2/` directory; see
[Follow-up tools](@ref "Orbit class: digest2")) and put the executable on
`PATH` or pass it as `digest2_path`. `config_dir` defaults to the directory
of the executable, symlinks resolved, and must hold `digest2.model.csv`,
`digest2.obscodes` and `MPC.config`, all used as they ship with `digest2`.

The input is [`mpc80_report`](@ref)'s output, sorted by `(id, epoch)` as
`digest2` requires. The command is built from an argument vector, never a
shell string.

Returns one row per tracklet (`id`) with `digest2`'s `RMS` (its
linear-motion fit residual, arcsec) and four 0-100 pseudo-probabilities:
`int_score` (any orbit class of MPC interest), `neo_score` (near-Earth
object), and `n22_score`/`n18_score` (NEO with H ≤ 22 / ≤ 18). The
"Other Possibilities" column is not parsed.

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
