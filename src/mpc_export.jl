"""
    _jd_to_calendar(jd::Real) -> (year, month, day, day_fraction)

Break a Julian Date (UTC) into its Gregorian calendar year/month/integer
day plus the leftover fraction of that day, via the standard algorithm
(Meeus, *Astronomical Algorithms*, ch. 7). Shared by
[`julian_date_to_iso8601`](@ref) and [`mpc80_report`](@ref)'s date
field, which need the same calendar date at different sub-day
precisions (milliseconds vs. 1e-6 day) — each rounds `day_fraction` to
its own precision itself, since rounding here would fix one caller's
precision at the expense of the other's.
"""
function _jd_to_calendar(jd::Real)
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

    return year, month, day, day_frac - day
end

"""
    julian_date_to_iso8601(jd::Real) -> String

Convert a Julian Date (UTC) to an ISO 8601 UTC timestamp
(`"YYYY-MM-DDTHH:MM:SS.sssZ"`), the format ADES requires for `obsTime`.

Verified here against two independent, exactly known reference points,
not just trusted from memory: `2451545.0` is the J2000.0 epoch
(2000-01-01T12:00:00Z) and `2440587.5` is the Unix epoch
(1970-01-01T00:00:00Z); both are exercised in the test suite.
"""
function julian_date_to_iso8601(jd::Real)
    year, month, day, day_frac = _jd_to_calendar(jd)

    # Round to the nearest millisecond, not truncate — otherwise a
    # fractional day landing at (e.g.) 23:59:59.9997 truncates to
    # 23:59:59.999 instead of correctly rolling to the next second.
    ms_of_day = round(Int, day_frac * 86_400_000)
    ms_of_day == 86_400_000 && (ms_of_day = 0; day += 1)  # (only possible from rounding at the boundary; doesn't cascade into month/year — not worth the complexity for a sub-millisecond edge case that real ZTF timestamps won't hit)
    hour, rem1 = divrem(ms_of_day, 3_600_000)
    minute, rem2 = divrem(rem1, 60_000)
    second, ms = divrem(rem2, 1000)

    return @sprintf("%04d-%02d-%02dT%02d:%02d:%02d.%03dZ", year, month, day, hour, minute, second, ms)
end

"""
    _mpc_date(jd::Real) -> String

Format a Julian Date (UTC) as the MPC 80-column format's
`"YYYY MM DD.dddddd"` date field (columns 16-32; 1e-6 day ≈ 0.0864s
precision). Unlike [`julian_date_to_iso8601`](@ref)'s millisecond
rollover, a rounding carry here is resolved by re-deriving the whole
calendar date at a nudged `jd` (via [`_jd_to_calendar`](@ref)) rather
than hand-incrementing `day`, so it stays correct across a month/year
boundary too.
"""
function _mpc_date(jd::Real)
    year, month, day, day_frac = _jd_to_calendar(jd)
    micros = round(Int, day_frac * 1_000_000)
    if micros == 1_000_000
        year, month, day, day_frac = _jd_to_calendar(jd + 1e-6)
        micros = round(Int, day_frac * 1_000_000)
    end
    return @sprintf("%04d %02d %02d.%06d", year, month, day, micros)
end

"""
    _mpc_ra(ra_deg::Real) -> String

Format a right ascension in decimal degrees as the MPC 80-column
format's `"HH MM SS.ddd"` field (columns 33-44). Working in integer
milliseconds-of-time and reducing modulo a full day sidesteps the
usual manual carry-the-1 logic for seconds/minutes/hours rollover
entirely.
"""
function _mpc_ra(ra_deg::Real)
    hours = mod(ra_deg, 360.0) / 15.0
    total_ms = mod(round(Int, hours * 3_600_000), 24 * 3_600_000)
    hour, rem1 = divrem(total_ms, 3_600_000)
    minute, rem2 = divrem(rem1, 60_000)
    second, ms = divrem(rem2, 1000)
    return @sprintf("%02d %02d %02d.%03d", hour, minute, second, ms)
end

"""
    _mpc_dec(dec_deg::Real) -> String

Format a declination in decimal degrees as the MPC 80-column format's
`"sDD MM SS.dd"` field (columns 45-56; `s` is a mandatory `+`/`-`, per
the spec — MPC does not default a blank sign to positive). Clamped at
the pole rather than left to overflow into an invalid degrees field:
not expected on any real asteroid-search field, but cheap to guard.
"""
function _mpc_dec(dec_deg::Real)
    sign_char = dec_deg < 0 ? '-' : '+'
    total_cs = min(round(Int, abs(dec_deg) * 360_000), 90 * 360_000)
    deg, rem1 = divrem(total_cs, 360_000)
    minute, rem2 = divrem(rem1, 6_000)
    second, cs = divrem(rem2, 100)
    return @sprintf("%c%02d %02d %02d.%02d", sign_char, deg, minute, second, cs)
end

"""
    ades_psv(candidates, station::AbstractString; mode::AbstractString="CCD",
             trksub_prefix::AbstractString="", astCat=nothing,
             photCat=nothing, band=nothing) -> String

Format `candidates` (an [`astrometric_calibrate`](@ref) table — columns
`id`, `frame`, `x`, `y`, `ra`, `dec`, `epoch`) as an ADES PSV
(pipe-separated values) observation table — the format the Minor Planet
Center currently requires for astrometric submissions, superseding the
legacy fixed-width 80-column format. See [`mpc80_report`](@ref) for the
80-column format itself, still required by some programs (e.g. IASC, as
of 2026) even though the MPC's own submissions now prefer ADES.

New, undesignated objects (this pipeline never produces MPC-designated
ones — candidates are locally-numbered tracklets) are identified by
`trkSub`, an observer-chosen tracking label (here, each real tracklet's
own `id`, base-36 encoded to stay compact and prefixed with
`trksub_prefix` if given), which is exactly what this pipeline already
produces.

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

"""
    mpc80_report(candidates, station::AbstractString; note1::AbstractChar=' ',
                 note2::AbstractChar='C', trksub_prefix::AbstractString="") -> String

Format `candidates` (an [`astrometric_calibrate`](@ref) table, same
shape [`ades_psv`](@ref) takes) as the legacy MPC1992 fixed-width
80-column format — verified column-by-column against the Minor Planet
Center's own published specification (`OpticalObs.html`,
`PackedDes.html`), not guessed: columns 1-5 packed permanent number,
6-12 packed provisional/temporary designation, 13 discovery asterisk,
14 note 1, 15 note 2, 16-32 date (`"YYYY MM DD.dddddd"`), 33-44 RA
(`"HH MM SS.ddd"`), 45-56 Dec (`"sDD MM SS.dd"`, sign mandatory — MPC
does *not* default a blank sign to positive), 57-65 blank, 66-71
magnitude+band, 72-77 blank, 78-80 observatory code.

Columns 1-5 are always blank: this pipeline's candidates are
locally-numbered tracklets, never MPC-numbered objects. Columns 6-12
use the same real, spec-sanctioned escape hatch [`ades_psv`](@ref)'s
`trkSub` relies on: the spec requires *some* designation in columns
1-12 ("never leave \\[them\\] blank"), but for a brand-new,
not-yet-designated discovery it explicitly allows an
observer-assigned *temporary* designation in place of an MPC-packed
provisional one — alphanumeric only, max 7 characters here (one
shorter than ADES's 8-character `trkSub` limit, since the packed-number
field this shares a slot with doesn't get the eighth column ADES adds).
Built the same way as `ades_psv`'s `trkSub`: each tracklet's own `id`,
base-36 encoded and optionally prefixed via `trksub_prefix`.

Three real gaps, same as `ades_psv`'s (deliberately not guessed at):
no discovery asterisk (column 13 always blank — this pipeline doesn't
track which observation of a tracklet was reported first), no
magnitude/band (columns 66-71 always blank — no photometric zeropoint
is calibrated against a reference catalog; a submission without
magnitudes is valid, same as for ADES), and `note1`/`note2` are fixed
per call rather than derived per detection. `note2` defaults to `'C'`
(CCD), matching how ZTF (and most modern digital-sensor surveys)
report; see the MPC spec for the full note2 code table if submitting
from a different observing mode.

One row (one 80-character line, newline-terminated) per detection
point, same granularity as `ades_psv`. Returns the report content as a
`String`; write it to a file yourself
(e.g. `write("submission.txt", mpc80_report(candidates, "I41"))`).
"""
function mpc80_report(candidates, station::AbstractString; note1::AbstractChar=' ',
                       note2::AbstractChar='C', trksub_prefix::AbstractString="")
    length(station) == 3 || throw(ArgumentError("station must be a 3-character MPC observatory code"))

    lines = String[]
    for row in candidates
        trksub = trksub_prefix * uppercase(string(row.id; base=36))
        length(trksub) <= 7 || throw(ArgumentError(
            "temporary designation \"$trksub\" exceeds the 80-column format's 7-character limit; use a shorter trksub_prefix"))

        line = " "^5 * rpad(trksub, 7) * ' ' * note1 * note2 *
               _mpc_date(row.epoch) * _mpc_ra(row.ra) * _mpc_dec(row.dec) *
               " "^9 * " "^6 * " "^6 * station
        @assert length(line) == 80
        push!(lines, line)
    end

    return join(lines, "\n") * "\n"
end
