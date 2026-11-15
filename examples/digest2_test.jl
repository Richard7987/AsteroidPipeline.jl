#=
How the MPC's digest2 scores the tracklets analyze_set finds on campaign
sets, known asteroids and unidentified alike:

    julia --project=. examples/digest2_test.jl <digest2 executable> <folder or .zip> ...

Declared before running:

Sets. The 19 Pan-STARRS1 sets of IASC's October 2026 campaign searched by
the NEO Hunters team (one folder or .zip each).

Search. analyze_set with its defaults, observatory F51, the MPC's orbits;
digest2 v1.0.0 from github.com/Smithsonian/mpc-public with its MPC.config,
through digest2_score.

Measured. For every identified tracklet, its NEO score and whether the
object is a near-Earth object by its MPC orbit (perihelion q < 1.3 au);
how many of those that are not score 65 or more, the score at which the
MPC posts a tracklet to its NEO Confirmation Page. For every unidentified
tracklet, its scores.
=#
using AsteroidPipeline, Printf, Statistics

length(ARGS) >= 2 || error("usage: julia --project=. examples/digest2_test.jl <digest2> <folder or .zip> ...")
digest2 = ARGS[1]
orbits = read_mpcorb(download_mpcorb())
bare(n) = strip(replace(n, r"^\(\d+\)\s*" => ""))
by_name = Dict{String,Any}()
for o in orbits
    by_name[bare(o.name)] = o
    by_name[o.packed] = o
end

known_rows = NamedTuple[]
unknown_rows = NamedTuple[]
for source in ARGS[2:end]
    work = mktempdir()
    r = analyze_set(abspath(source); workdir=work, extract_to=work, observatory="F51", mpcorb=orbits, digest2)
    @printf("%s: %s, %d known, %d unidentified\n", r.name, r.reliable ? "reliable" : "unreliable",
            length(r.known), length(r.unknown_ids))
    rows = filter(c -> haskey(r.known, c.id), collect(r.candidates))
    for s in (isempty(rows) ? [] : digest2_score(rows, "F51"; digest2_path=digest2))
        name = r.known[s.id]
        o = get(by_name, bare(name), nothing)
        q = o === nothing ? NaN : o.a * (1 - o.e)
        push!(known_rows, (set=r.name, name=name, q=q, neo=s.neo_score, int=s.int_score, rms=s.rms))
        @printf("  known %-22s q %5.2f au  NEO %3.0f  Int %3.0f  RMS %.2f\"\n", name, q, s.neo_score, s.int_score, s.rms)
    end
    for (id, s) in sort(collect(r.digest2); by=first)
        push!(unknown_rows, (set=r.name, id=id, neo=s.neo_score, int=s.int_score, n22=s.n22_score, rms=s.rms))
        @printf("  unidentified %d  NEO %3.0f  Int %3.0f  N22 %3.0f  RMS %.2f\"\n", id, s.neo_score, s.int_score,
                s.n22_score, s.rms)
    end
end

nonneo = filter(k -> isfinite(k.q) && k.q >= 1.3, known_rows)
neo = filter(k -> isfinite(k.q) && k.q < 1.3, known_rows)
@printf("\n%d identified tracklets: %d not NEOs, %d NEOs, %d without an MPC orbit found\n",
        length(known_rows), length(nonneo), length(neo), count(k -> !isfinite(k.q), known_rows))
isempty(nonneo) || @printf("not NEOs: NEO score median %.0f, max %.0f; %d at 65 or more\n",
                           median(k.neo for k in nonneo), maximum(k.neo for k in nonneo), count(k -> k.neo >= 65, nonneo))
for k in neo
    @printf("NEO %s: NEO score %.0f\n", k.name, k.neo)
end
@printf("%d unidentified tracklets, %d at NEO 65 or more\n", length(unknown_rows), count(u -> u.neo >= 65, unknown_rows))
