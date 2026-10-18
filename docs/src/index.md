```@raw html
---
next:
  text: Getting started
  link: /getting-started
---
```

# AsteroidPipeline.jl

**Find moving objects in astronomical images.** AsteroidPipeline.jl takes a
short sequence of images of one patch of sky — a few exposures taken minutes
apart — and finds the asteroids (and comets) moving through it: it calibrates
the images against the Gaia catalogue, detects what changes from one image to
the next, links detections that move along straight lines, identifies the
ones that are already known, and writes the observations of the rest in the
format the Minor Planet Center (MPC) accepts.

![The pipeline's stages, from an image set to an MPC report](assets/pipeline-overview.png)

It works with the FITS images of any telescope that come with an astrometric
solution (a WCS) in their headers, or can be given one. It has been tested on
real data from Pan-STARRS1 and the Zwicky Transient Facility.

## In a few lines

```julia
using AsteroidPipeline

result = analyze_set("images.zip"; observatory="F51", workdir="work")

result.unknown_ids        # tracklets no catalogue knows: candidates
result.expected           # known asteroids predicted in the field
print(result.report)      # the MPC report (blank if nothing new was found)
```

## Where to go next

| If you want to… | Read |
|:--|:--|
| try it on your own images | [Getting started](getting-started.md) |
| do a specific task (reports, known objects, another telescope…) | the [Guides](guides/image-sets.md) |
| understand how it decides what is real | [How it works](concepts/pipeline.md) |
| know how well it works, and where it doesn't | [Validation](validation/summary.md) |
| look up a function | the [Reference](reference/sets.md) |

## Status

Version 0.2. Every stage has been run on real telescope data and checked
against independent results — catalogued asteroids, another program's
measurements, and synthetic asteroids injected into real images. See the
[validation summary](validation/summary.md). Limits are listed in
[Limitations](concepts/limitations.md).

AsteroidPipeline.jl is free software under the MIT license.
