# Use a different telescope

Nothing in the pipeline is specific to one telescope; three things depend on
which one took the images.

## 1. The observatory code

Reports start with it, and known objects are predicted as seen from it. The
images usually say where they were taken, and [`frame_site`](@ref) reads it;
[`observatory_candidates`](@ref) finds the MPC codes that fit:

```julia
using FITSIO
site = FITS(f -> frame_site(read_header(f[1])), "frame1.fits")
codes = read_obscodes(download_obscodes())
observatory_candidates(site, codes)[1:3]
```

| Images | Header says | Best match |
|:--|:--|:--|
| Pan-STARRS1 | `TELESCOP=PS1`, site in `FPA.LONGITUDE` (west, hours)… | **F51** Pan-STARRS 1 (F52, Pan-STARRS 2, is 10 m away: the name decides) |
| ZTF | `TELESCOP=Palomar 48-inch`, `ORIGIN=Zwicky Transient Facility`, `OBSLON`… | **I41** Palomar Mountain--ZTF |

Codes within a few kilometres of the site are ranked by how well their names
match the telescope's; without a site in the header, by name alone. Check the
result: two telescopes on one mountain can only be told apart by name.

## 2. The pixel scale

Search tolerances are angles (a 2″ linking radius, rates of 3–210″/h); the
pipeline works in pixels. [`analyze_set`](@ref) measures the scale from the
images and calls [`pipeline_config`](@ref) with it. With
[`run_pipeline`](@ref), do the same:

```julia
run_pipeline(paths; pipeline_config(1.01)...)    # ZTF: 1.01"/pixel
```

## 3. Header conventions

| Need | Read from | If missing |
|:--|:--|:--|
| time of each exposure | the first of `DATE-OBS` (with `TIME-OBS`, `UT` or `UTC` if it holds only the date), `MJD-OBS`, `OBSMJD`, `MJD_OBS`, `JD`, `OBSJD`, `JD-OBS` — the start of the exposure; plus half of `EXPTIME` ([`frame_epoch`](@ref)) | pass `timestamp_key` to [`run_pipeline`](@ref) |
| astrometric solution | the header's WCS, refitted to Gaia DR3, all frames together; the header may be off by any amount up to ~75″ | [`plate_solve`](@ref) (see [Follow-up tools](follow-up.md)) |
| known objects | SkyBoT; and the MPC's orbit file for images within a year of its orbits' epoch — older images (an archive) use SkyBoT alone | — |
| invalid pixels | `BLANK`, constant fill values | handled by [`clean_blank_pixels`](@ref) and `mask_fill_values` |

Check the time keywords first: a header whose `MJD-OBS` is in TAI rather than
UTC (Pan-STARRS1's is) puts every position 37 s out — ~0.4″ for a typical
main-belt asteroid.
