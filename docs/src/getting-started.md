# Getting started

This tutorial takes one image set from files on disk to an MPC report. It
needs about fifteen minutes, most of it Julia installing packages the first
time.

## What you need

- **Julia 1.11 or newer** — [julialang.org/install](https://julialang.org/install).
- **An image set**: three or more FITS images of the same field, taken over
  tens of minutes to a few hours, each with a WCS (astrometric solution) in
  its header. A `.zip` of them, or a folder, both work. Pan-STARRS1, ZTF,
  and most survey and amateur pipelines provide this.
- **The MPC code of the observatory** that took them — `F51` for
  Pan-STARRS1, `I41` for ZTF; the [MPC's list](https://minorplanetcenter.net/iau/lists/ObsCodes.html)
  has them all. If you don't know it, [Use a different telescope](guides/telescopes.md)
  shows how to find it from the images.
- **An internet connection**: positions are calibrated against Gaia DR3,
  and known asteroids come from the MPC and SkyBoT.

## 1. Install

```julia
using Pkg
Pkg.add(url="https://github.com/Richard7987/AsteroidPipeline.jl")
```

## 2. Analyze the set

```julia
using AsteroidPipeline

result = analyze_set("ps1-20261005_9_XY45_p10.zip";
                     observatory="F51",
                     workdir="work",
                     mpcorb=read_mpcorb(download_mpcorb()))
```

What happens, in order (each step is explained in [How it works](concepts/pipeline.md)):

| Step | Takes |
|:--|:--|
| unpack the `.zip` next to it | seconds |
| download the MPC's orbit file (once a day, ~70 MB) and read it | ~1 min the first time |
| refit each image's astrometry to Gaia DR3 stars | ~10 s per image |
| search for moving objects | ~30 s |
| identify the known ones, and check the rest against everything known nearby | ~30 s |

The `mpcorb` argument is optional; without it, known objects come from
SkyBoT alone, which can lag behind the MPC's newest orbits.

## 3. Read the result

```julia
julia> result.reliable                 # false: the images are too damaged to search
true

julia> result.known                    # tracklet id => the known object it is
Dict{Int64, String} with 2 entries:
  2 => "2015 BK618"
  1 => "(656932) 2016 FA12"

julia> result.unknown_ids              # candidates: nothing known moves like them
Int64[]

julia> [(o.name, round(o.mag; digits=1), o.detected) for o in result.expected]
2-element Vector{Tuple{String, Float64, Bool}}:
 ("(656932) 2016 FA12", 19.8, true)
 ("2015 BK618", 20.1, true)
```

`result.candidates` holds every measured position (right ascension,
declination, time, magnitude) of every tracklet; `result.nearby` says, for
each candidate, which known objects lie within 10′ and whether one of them
is probably it.

## 4. Write the report

With no candidates, `result.report` is the blank report — what is submitted
for an image set with nothing new in it:

```
COD F51
ACK MPCReport file updated 2026.10.18 13:16:23
NET Gaia DR3
----- end -----
```

With candidates, it holds their observations. **Look at every candidate in
the images before submitting it** — [Telling real objects from artifacts](concepts/artifacts.md)
shows what false ones look like. To add the report's header lines (who
observed, who measured, with what):

```julia
result = analyze_set("images.zip"; observatory="F51", workdir="work",
                     header=(observers="J. Bulger, T. Lowe", measurers="A. Name",
                             telescope="1.8-m f/4.4 Ritchey-Chretien + CCD"))
```

## Next

- Many image sets at once: [Analyze image sets](guides/image-sets.md).
- Is a candidate really new? [Identify known objects](guides/known-objects.md).
- Reports in ADES, and checking a report before submitting it:
  [Write and check MPC reports](guides/reports.md).
