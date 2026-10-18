#=
Analyze every image set in a folder at once — .zip files of FITS frames, or
sub-folders of FITS frames — and print one line per set: what was found, and
which sets have candidates for a person to check. For each set, its report
(the blank one, or the candidates' observations for checking) is written
next to it as <set>_MPCReport.txt.

    julia --project=. examples/analyze_sets.jl ~/Downloads/set4 F51

The second argument is the observatory's MPC code (F51 = Pan-STARRS1). The
MPC's orbit file is downloaded once a day into ~/.cache/AsteroidPipeline (a
~70 MB download, unpacked to ~320 MB); without it, known objects come from
SkyBoT alone. The interactive interface in ui/ does the same per set, and
adds the images and the human check.
=#
using AsteroidPipeline
using Printf

length(ARGS) == 2 || error("usage: julia --project=. examples/analyze_sets.jl <folder> <observatory code>")
folder, observatory = abspath(ARGS[1]), uppercase(ARGS[2])
isdir(folder) || error("no folder $folder")

is_fits(f) = any(e -> endswith(lowercase(f), e), (".fits", ".fit", ".fts"))
# The .zip files if there are any — their sets get unpacked into sub-folders
# here, which must not count twice — else the sub-folders of FITS frames.
sources = [joinpath(folder, f) for f in readdir(folder) if endswith(lowercase(f), ".zip")]
if isempty(sources)
    sources = [joinpath(folder, f) for f in readdir(folder)
               if isdir(joinpath(folder, f)) && count(is_fits, readdir(joinpath(folder, f))) >= 3]
end
isempty(sources) && error("no .zip files or folders of FITS frames in $folder")

mpcorb = try
    read_mpcorb(download_mpcorb())
catch e
    @warn "no MPC orbit file; known objects from SkyBoT only" exception = e
    nothing
end

workdir = joinpath(homedir(), ".cache", "AsteroidPipeline", "sets")
summary = []
for (i, source) in enumerate(sources)
    println("[$i/$(length(sources))] $(basename(source))")
    r = try
        analyze_set(source; workdir, observatory, mpcorb)
    catch e
        @error "analysis failed" source exception = e
        push!(summary, (basename(source), "FAILED: $(sprint(showerror, e))", 0, 0, 0))
        continue
    end
    out = joinpath(folder, "$(r.name)_MPCReport.txt")
    write(out, r.report)
    verdict = r.report_kind == :unreliable ? "damaged image set" :
              r.report_kind == :blank ? "no new objects" : "CANDIDATES TO CHECK"
    push!(summary, (r.name, verdict, length(r.unknown_ids), length(r.known), length(r.frames_used)))
end

println()
@printf("%-24s %-22s %10s %6s %7s\n", "set", "result", "candidates", "known", "frames")
for (name, verdict, n_new, n_known, n_used) in summary
    @printf("%-24s %-22s %10d %6d %7d\n", name, verdict, n_new, n_known, n_used)
end
