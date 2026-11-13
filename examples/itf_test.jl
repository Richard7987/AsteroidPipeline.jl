#=
How many of the tracklets analyze_set leaves unidentified had already been
reported to the MPC, unlinked, and sit in its Isolated Tracklet File:

    julia --project=. examples/itf_test.jl <folder or .zip> ...

Declared before running:

Sets. The 19 Pan-STARRS1 sets of IASC's October 2026 campaign searched by
the NEO Hunters team (one folder or .zip each).

Search. analyze_set with its defaults, observatory F51, the MPC's orbits
and the ITF as downloaded on the day of the run (download_itf).

Measured. The unidentified tracklets, and those itf_check finds in the ITF
(2" from the tracklet's straight-line motion, within 30 min of it): for
each, the ITF tracklet, its station, the distance and how far in time each
ITF observation is from the nearest of the set's exposures (the same
exposure, reported by Pan-STARRS itself, is seconds away).
=#
using AsteroidPipeline, Printf, Statistics

isempty(ARGS) && error("usage: julia --project=. examples/itf_test.jl <folder or .zip> ...")
orbits = read_mpcorb(download_mpcorb())
itf_path = download_itf()
@printf("ITF of %s\n", Libc.strftime("%Y-%m-%d", mtime(itf_path)))

unidentified = 0
in_itf = 0
for source in ARGS
    work = mktempdir()
    r = analyze_set(abspath(source); workdir=work, extract_to=work, observatory="F51", mpcorb=orbits, itf=itf_path)
    global unidentified += length(r.unknown_ids)
    global in_itf += length(r.itf)
    @printf("%s: %s, %d known, %d unidentified, %d in the ITF\n", r.name, r.reliable ? "reliable" : "unreliable",
            length(r.known), length(r.unknown_ids), length(r.itf))
    for (id, matches) in sort(collect(r.itf); by=first)
        for m in matches
            dt = minimum(abs(m.jd - e) for e in r.epochs) * 86400
            @printf("  tracklet %d: %s %s, %.2f\" from its motion, %.0f s from an exposure\n",
                    id, m.trksub, m.station, m.distance_arcsec, dt)
        end
    end
end
@printf("\n%d unidentified tracklets on %d sets, %d of them in the ITF\n", unidentified, length(ARGS), in_itf)
