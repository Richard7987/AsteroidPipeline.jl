# Telling real objects from artifacts

A difference image keeps everything that changed. Most of it isn't an
asteroid. Every filter below exists because a false candidate of that kind
reached a person on real data, and each was set by measuring what separates
it from real objects — not by guessing.

## Hot pixels and detector defects

A defective pixel is bright in every frame, at the same place **on the
detector**. Between exposures the telescope moves slightly (dithers), so on
the sky the defect moves — and a defect that "moves" a few pixels between
frames looks like a slow asteroid.

| | Movement on the detector across a set |
|:--|:--|
| real objects measured | 32 px or more |
| the hot pixel linked as a slow mover | 0 px |

`detector_artifact_radius` drops a detection that recurs within that many
pixels of the same detector position in two or more other frames.

## Light of bright stars

Around a bright star, the difference images keep residuals: its halo, an
internal reflection of it (a *ghost*), and its diffraction spikes, which
rotate with the telescope and so don't cancel in the median. They all sit on
light a faint real object doesn't have. Measured as the median of an 8–20 px
ring around the detection, in the raw frame, above the sky:

![Background excess around real objects and artifacts](../assets/background-excess.png)

`max_background_excess` (1σ) drops detections above the cut. It also costs
real objects that happen to pass close to bright stars: 2 of 180 injected
asteroids in a completeness test.

## A damaged image set

Some sets aren't worth searching: defective detector cells fill a frame with
speckles and gradients. Such a frame has far more sources than stars: on a
real set, 137–1126 per Gaia star, against 0.6–5.4 on every other set (see
[the pipeline](pipeline.md)). The frame is left out; a set without three good
frames is flagged `reliable = false` and not searched.

## What a real candidate looks like

When checking candidates by eye, a real moving object:

- is in every frame (or all but one), on a **straight line**, at **even
  spacing** (the exposures' timing);
- looks like the **stars around it** — same size and shape, or slightly
  trailed along its motion;
- keeps roughly the **same brightness**;
- is **not** on a star's halo or spike, at the edge of a gap, or on a
  column of bad pixels.

A candidate that fails one of these is almost always an artifact; one that
passes them all, and that no known object explains
([Is a candidate new?](known-objects.md)), is worth reporting.
