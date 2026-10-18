# Identify known objects

Most moving objects in an image set are asteroids someone has already
catalogued. Three tools tell them apart from new ones; [`analyze_set`](@ref)
uses all three.

## Is this tracklet a known object?

[`crossmatch_catalog`](@ref) compares each measured point with the positions
catalogued objects are predicted at, at that point's own time, from the
observatory:

```julia
skybot = crossmatch_catalog(candidates, :skybot; radius=5.0, observatory="F51")
orbits = read_mpcorb(download_mpcorb())
mpc = crossmatch_catalog(candidates, :mpcorb; radius=15.0, mpcorb=orbits, observatory="F51")
```

| Catalogue | Source | Up to date? | Use |
|:--|:--|:--|:--|
| `:skybot` | [SkyBoT](https://vo.imcce.fr/webservices/skybot/) (IMCCE), a web service | can lag days behind the MPC | quick, no download |
| `:mpcorb` | the MPC's own orbit file, propagated here (two-body, light-time corrected) | as of the file's day | the current orbits; needs the ~70 MB file |
| `:vsx`, `:simbad` | variable stars, other sources | — | for stationary variables |

A tracklet counts as an object when it matches in at least three frames — a
star near an asteroid's path matches in one.

!!! tip "Use both"
    On a real set, SkyBoT placed one asteroid (2019 UA185) far from where it
    was, using an orbit the MPC had replaced that day; the MPC's file put it
    within 4″. [`download_mpcorb`](@ref) keeps a copy no older than a day.

## Is a candidate *really* new?

A known object with an imprecise orbit can be predicted further away than
these radii, and then shows up as a candidate. [`nearby_known_objects`](@ref)
looks at everything known within 10′ and compares motions:

```julia
check = nearby_known_objects(candidate_rows; observatory="F51", mpcorb=orbits)
check.checks[2].verdict      # :new, :possibly_known or :known
check.checks[2].objects      # name, separation, predicted rate and direction, agreement
```

How the comparison works, and why it has to be strict, is in
[Is a candidate new?](../concepts/known-objects.md).

## Is someone else observing it?

Objects reported in the last days that have no orbit yet are on the MPC's
[NEO Confirmation Page](https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html)
(NEOCP) — and in no catalogue. [`neocp_check`](@ref) asks the MPC for their
positions at the candidate's own times, from the observatory:

```julia
c = neocp_check(candidate_rows; observatory="F51")
c.checks[2].verdict        # :on_neocp or :none
c.checks[2].objects        # NEOCP designation, separation, predicted motion, agreement
```

A candidate on the NEOCP is someone's recent, unconfirmed find: if it is
real, report it under the NEOCP's designation — the observations help
confirm it. Only recent images are checked (by default, up to 30 days old):
the NEOCP lists what is being observed now, and its ephemerides come from
arcs of hours, so they are uncertain by arcminutes (30′ are searched).

## Which known objects should be in my images?

[`known_objects_in_field`](@ref) lists every catalogued object predicted on
the frames — including the ones too faint to detect — with its position in
each frame's pixels: what to look for when checking a set by eye.

```julia
field = known_objects_in_field(r.wcs, r.epochs, r.image_size; mpcorb=orbits, observatory="F51")
[(o.name, round(o.mag; digits=1)) for o in field.objects]
```
