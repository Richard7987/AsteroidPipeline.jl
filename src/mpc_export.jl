"""
    julian_date_to_iso8601(jd::Real) -> String

Convert a Julian Date (UTC) to an ISO 8601 UTC timestamp
(`"YYYY-MM-DDTHH:MM:SS.sssZ"`), the format ADES requires for `obsTime`.

Uses the standard Gregorian-calendar algorithm (Meeus, *Astronomical
Algorithms*, ch. 7) — verified here against two independent, exactly
known reference points, not just trusted from memory: `2451545.0` is
the J2000.0 epoch (2000-01-01T12:00:00Z) and `2440587.5` is the Unix
epoch (1970-01-01T00:00:00Z); both are exercised in the test suite.
"""
function julian_date_to_iso8601(jd::Real)
    Z = floor(Int, jd + 0.5)
    F = (jd + 0.5) - Z

    if Z < 2299161
        A = Z
    else
        alpha = floor(Int, (Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - floor(Int, alpha / 4)
    end

    B = A + 1524
    C = floor(Int, (B - 122.1) / 365.25)
    D = floor(Int, 365.25 * C)
    E = floor(Int, (B - D) / 30.6001)

    day_frac = B - D - floor(Int, 30.6001 * E) + F
    day = floor(Int, day_frac)
    month = E < 14 ? E - 1 : E - 13
    year = month > 2 ? C - 4716 : C - 4715

    # Round to the nearest millisecond, not truncate — otherwise a
    # fractional day landing at (e.g.) 23:59:59.9997 truncates to
    # 23:59:59.999 instead of correctly rolling to the next second.
    ms_of_day = round(Int, (day_frac - day) * 86_400_000)
    ms_of_day == 86_400_000 && (ms_of_day = 0; day += 1)  # (only possible from rounding at the boundary)
    hour, rem1 = divrem(ms_of_day, 3_600_000)
    minute, rem2 = divrem(rem1, 60_000)
    second, ms = divrem(rem2, 1000)

    return @sprintf("%04d-%02d-%02dT%02d:%02d:%02d.%03dZ", year, month, day, hour, minute, second, ms)
end

"""
    ades_psv(candidates, station::AbstractString; mode::AbstractString="CCD",
             trksub_prefix::AbstractString="", astCat=nothing,
             photCat=nothing, band=nothing) -> String

Format `candidates` (an [`astrometric_calibrate`](@ref) table — columns
`id`, `frame`, `x`, `y`, `ra`, `dec`, `epoch`) as an ADES PSV
(pipe-separated values) observation table — the format the Minor Planet
Center currently requires for astrometric submissions, superseding the
legacy fixed-width 80-column format.

Only the legacy 80-column format's replacement (ADES) is implemented
here, not the 80-column format itself: 80-column records pack a
provisional designation into a specific fixed encoding that needs a
real MPC-assigned designation to round-trip correctly, which nothing in
this pipeline has (candidates are locally-numbered tracklets, not
MPC-designated objects) — guessing at that packing without real
reference examples to check against risked producing output that reads
as well-formed but is subtly wrong, exactly the kind of mistake that
matters for a real submission. ADES has no such requirement: new,
undesignated objects are identified by `trkSub`, an observer-chosen
tracking label (here, each real tracklet's own `id`, base-36 encoded to
stay compact and prefixed with `trksub_prefix` if given), which is
exactly what this pipeline already produces.

One row per detection point (i.e. one row per `candidates` row, not one
per tracklet) — this is the granularity ADES observation records use;
the Minor Planet Center correlates same-`trkSub` rows into a tracklet on
its own end. `station` is the observer's MPC-assigned station/observatory
code (3 characters, e.g. ZTF's is `"I41"`) and must be supplied — there
is no way to derive it from pixel data. `astCat`/`photCat`/`band` are
the astrometric reference catalog, photometric reference catalog, and
photometric band used, respectively; all `nothing` (omitted from the
output) by default, since this pipeline does not itself calibrate a
photometric zeropoint or track which catalog `load_wcs`'s astrometric
solution was fit against — real gaps, not filled in with a guessed
value. A submitted ADES file without magnitudes is valid; MPC accepts
astrometry-only submissions.

Returns the PSV content as a `String`; write it to a `.psv` file
yourself (e.g. `write("submission.psv", ades_psv(candidates, "I41"))`).
"""
function ades_psv(candidates, station::AbstractString; mode::AbstractString="CCD",
                   trksub_prefix::AbstractString="",
                   astCat::Union{Nothing,AbstractString}=nothing,
                   photCat::Union{Nothing,AbstractString}=nothing,
                   band::Union{Nothing,AbstractString}=nothing)
    length(station) == 3 || throw(ArgumentError("station must be a 3-character MPC observatory code"))

    columns = ["trkSub", "mode", "stn", "obsTime", "ra", "dec"]
    astCat !== nothing && push!(columns, "astCat")
    band !== nothing && push!(columns, "band")
    photCat !== nothing && push!(columns, "photCat")

    lines = [join(columns, "|")]
    for row in candidates
        trksub = trksub_prefix * uppercase(string(row.id; base=36))
        length(trksub) <= 8 || throw(ArgumentError(
            "trkSub \"$trksub\" exceeds ADES's 8-character limit; use a shorter trksub_prefix"))

        fields = [trksub, mode, station, julian_date_to_iso8601(row.epoch),
                  @sprintf("%.7f", row.ra), @sprintf("%.7f", row.dec)]
        astCat !== nothing && push!(fields, astCat)
        band !== nothing && push!(fields, band)
        photCat !== nothing && push!(fields, photCat)
        push!(lines, join(fields, "|"))
    end

    return join(lines, "\n") * "\n"
end
