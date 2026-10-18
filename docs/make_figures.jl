#=
Generates the figures embedded in the documentation, from the
real numbers measured during this project's real-data investigations (not
synthetic/illustrative data) — see the corresponding prose sections for
where each number comes from. Run automatically by docs/make.jl before
`makedocs`, so CI regenerates these on every docs build rather than
committing static images that could drift from the numbers in the text.

Styled for the wiki's dark-only theme (docs/src/.vitepress/config.mts,
appearance: 'force-dark') — transparent background, white text — rather
than a generic light-background default.
=#
using CairoMakie

const ASSETS_DIR = joinpath(@__DIR__, "src", "assets")
mkpath(ASSETS_DIR)

set_theme!(Theme(
    backgroundcolor=:transparent,
    textcolor=:white,
    Axis=(
        backgroundcolor=:transparent,
        xlabelcolor=:white, ylabelcolor=:white, titlecolor=:white,
        xticklabelcolor=:white, yticklabelcolor=:white,
        xtickcolor=:white, ytickcolor=:white,
        leftspinecolor=:white, rightspinecolor=:white,
        bottomspinecolor=:white, topspinecolor=:white,
        xgridcolor=(:white, 0.15), ygridcolor=(:white, 0.15),
    ),
    Legend=(backgroundcolor=:transparent, labelcolor=:white, framecolor=(:white, 0.3)),
))

# --- find_variable_sources: systematic_error_fraction sweep -----------
# Real measurements from sweeping systematic_error_fraction against three
# independent real ZTF fields, each with its own real matched stationary
# stars (false-positive rate at chi2_threshold=10.0) and its own real,
# independently-confirmed variable star (reduced chi2): field 451 (152
# stars; variable ASASSN-V J183620.31, from a different real field), and
# two more found and validated later, V1012 Mon (197 stars) and ASASSN-V
# J072906.85-090518.2 (182 stars). See "Cross-validating
# find_variable_sources's systematic error floor..." in this log.
let
    floors = [0.0, 1.0, 2.0, 3.0, 5.0]  # percent, common to all three fields
    field451_fp = [13.2, 3.3, 2.0, 1.3, 0.0]
    field451_chi2 = [20688.5, 123.4, 31.0, 13.8, 5.0]
    v1012_fp = [62.24, 6.12, 2.55, 2.04, 1.02]
    v1012_chi2 = [31908.9, 542.3, 137.4, 61.2, 22.1]
    asassn_fp = [49.72, 2.76, 0.55, 0.55, 0.55]
    asassn_chi2 = [3298.5, 151.2, 39.2, 17.5, 6.3]

    fig = Figure(size=(760, 520))
    ax1 = Axis(fig[1, 1];
        xlabel="systematic_error_fraction (%)",
        ylabel="false-positive rate at chi2_threshold=10 (%)",
        title="Three real, independent ZTF fields: false-positive rate vs. systematic_error_fraction")
    ax2 = Axis(fig[2, 1];
        xlabel="systematic_error_fraction (%)",
        ylabel="reduced chi² (real confirmed variable)", yscale=log10,
        title="Each field's own real confirmed variable's reduced chi²")

    colors = (:tomato, :mediumspringgreen, :dodgerblue)
    labels = ("field 451 (152 stars)", "V1012 Mon (197 stars)", "ASASSN-V J072906.85 (182 stars)")

    for (fp, chi2, color, label) in zip(
        (field451_fp, v1012_fp, asassn_fp), (field451_chi2, v1012_chi2, asassn_chi2), colors, labels)
        lines!(ax1, floors, fp; color, linewidth=2, label)
        scatter!(ax1, floors, fp; color, markersize=10)
        lines!(ax2, floors, chi2; color, linewidth=2, label)
        scatter!(ax2, floors, chi2; color, markersize=10)
    end
    hlines!(ax2, [10.0]; color=:white, linestyle=:dash, linewidth=1)
    vlines!(ax1, [3.0]; color=:white, linestyle=:dot, linewidth=1)
    vlines!(ax2, [3.0]; color=:white, linestyle=:dot, linewidth=1)
    text!(ax1, 3.1, 55.0; text="chosen default (3%)", fontsize=12, color=:white)

    Legend(fig[3, 1], ax1; orientation=:horizontal, tellwidth=false)

    save(joinpath(ASSETS_DIR, "systematic-error-floor-sweep.png"), fig)
end

# --- IASC match_radius retuning: tracklet counts before/after ---------
# Real per-field tracklet counts from examples/iasc_demo.jl, before
# (match_radius converted from ZTF's pixel scale, ~10" tolerance) and
# after (match_radius from PS1's own PERROR, ~2" tolerance) retuning. The
# same 9 distinct known objects were recovered in both runs — this is a
# reduction in spurious duplicate tracklets, not lost detections. See
# "match_radius was too loose, by a measured, corrected amount" above.
let
    fields = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11"]
    before = [524, 776, 1817, 2619, 10422]
    after = [52, 132, 567, 731, 3478]

    fig = Figure(size=(720, 440))
    ax = Axis(fig[1, 1];
        xlabel="IASC practice field", ylabel="tracklets found",
        title="Retuning match_radius to PS1's own real astrometric precision (PERROR)",
        xticks=(1:length(fields), fields), yscale=log10)

    n = length(fields)
    x = 1:n
    barplot!(ax, x .- 0.15, before; width=0.3, color=:salmon, label="match_radius ≈ 10\" (ZTF's tolerance, reused as-is)")
    barplot!(ax, x .+ 0.15, after; width=0.3, color=:mediumspringgreen, label="match_radius = 2\" (10x PS1's real PERROR)")
    # A log axis has no true zero, so bars would otherwise auto-start from
    # whatever the smallest value happens to be (misleadingly making every
    # bar look like a similar-height floating block) — pin the bottom
    # near 1 instead, the closest a log axis gets to a real baseline.
    ylims!(ax, 1, 20000)

    axislegend(ax; position=:lt)
    save(joinpath(ASSETS_DIR, "iasc-match-radius-retuning.png"), fig)
end


# --- The pipeline at a glance ----------------------------------------------
# Not data: the order of the stages, for the Home and "How it works" pages.
let
    stages = [("Image set", ".zip, folder\nor FITS files"),
              ("Clean &\ncalibrate", "WCS refit to\nGaia DR3"),
              ("Detect", "on each frame\nminus the median"),
              ("Filter", "hot pixels,\nbright-star light"),
              ("Link", "straight lines,\nconstant rate"),
              ("Identify", "SkyBoT,\nMPC orbits"),
              ("Report", "MPC 80-column\nor ADES")]
    fig = Figure(size=(1100, 230))
    ax = Axis(fig[1, 1]; limits=(0.4, length(stages) + 0.6, 0, 1))
    hidedecorations!(ax); hidespines!(ax)
    for (i, (title, detail)) in enumerate(stages)
        poly!(ax, Rect(i - 0.42, 0.22, 0.84, 0.6); color=(:dodgerblue, 0.25), strokecolor=:dodgerblue, strokewidth=2)
        text!(ax, i, 0.62; text=title, align=(:center, :center), fontsize=15, font=:bold, color=:white)
        text!(ax, i, 0.36; text=detail, align=(:center, :center), fontsize=11, color=(:white, 0.8))
        i < length(stages) && arrows2d!(ax, [i + 0.43], [0.52], [0.14], [0.0]; color=:white, shaftwidth=2, tipwidth=8, tiplength=8)
    end
    save(joinpath(ASSETS_DIR, "pipeline-overview.png"), fig)
end

# --- Completeness, measured by injection --------------------------------------
# examples/injection_test.jl on a real Pan-STARRS1 set (XY54_p10), with the
# configuration in pipeline_config: 180 main-belt-like objects (V 18-23,
# 5-100"/h) and 80 fast, NEO-like ones (100-1000"/h), injected into the real
# frames with each frame's own PSF, zero point, noise and trailing.
let
    mag_edges = 18.0:0.5:23.0
    mag_rec = [13, 14, 11, 13, 9, 16, 16, 13, 3, 0]
    mag_n = [17, 16, 17, 16, 14, 18, 19, 18, 19, 26]
    speed_labels = ["5-10", "10-20", "20-40", "40-70", "70-100", "100-200", "200-400", "400-700", "700-1000"]
    speed_rec = [5, 13, 22, 40, 28, 6, 7, 7, 2]
    speed_n = [10, 21, 31, 67, 51, 14, 16, 22, 28]

    fig = Figure(size=(1000, 420))
    ax1 = Axis(fig[1, 1]; xlabel="V magnitude", ylabel="recovered (%)", title="By brightness (5-100\"/h)",
               xticks=18:23, limits=(nothing, (0, 105)))
    centres = collect(mag_edges[1:end-1]) .+ 0.25
    barplot!(ax1, centres, 100 .* mag_rec ./ mag_n; width=0.45, color=:mediumspringgreen)
    ax2 = Axis(fig[1, 2]; xlabel="rate (\"/h)", ylabel="recovered (%)", title="By rate of motion",
               xticks=(1:length(speed_labels), speed_labels), xticklabelrotation=π / 5, limits=(nothing, (0, 105)))
    barplot!(ax2, 1:length(speed_labels), 100 .* speed_rec ./ speed_n;
             color=[i <= 5 ? :mediumspringgreen : :orange for i in eachindex(speed_rec)])
    text!(ax2, 7.5, 92; text="NEO-like", align=(:center, :center), color=:orange, fontsize=13)
    save(joinpath(ASSETS_DIR, "completeness.png"), fig)
end

# --- Set quality: sources per reference star ------------------------------------
# set_quality's ratio (10σ sources per Gaia DR3 star in the field), the worst
# frame of each of ten real Pan-STARRS1 sets. One set (XY73_p00) is full of
# defective detector cells; the limit is 30.
let
    sets = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY54_p10",
            "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
    worst = [0.9, 0.7, 1.0, 0.9, 1.7, 1.0, 1.1, 5.4, 1126.4, 2.6]
    fig = Figure(size=(900, 400))
    ax = Axis(fig[1, 1]; ylabel="10σ sources per reference star (worst frame)", yscale=log10,
              xticks=(1:length(sets), sets), xticklabelrotation=π / 5, title="Set quality on ten real image sets",
              limits=(nothing, (0.3, 3000)))
    barplot!(ax, 1:length(sets), worst; color=[w > 30 ? :tomato : :mediumspringgreen for w in worst],
             fillto=0.3)
    hlines!(ax, [30.0]; color=:white, linestyle=:dash)
    text!(ax, 0.6, 40; text="limit: 30", color=:white, fontsize=12)
    save(joinpath(ASSETS_DIR, "set-quality.png"), fig)
end

# --- Bright-star light around a detection ------------------------------------------
# The median of an 8-20 px annulus around each detection, in the raw frame,
# above the sky in units of its noise: every real object measured on the
# Pan-STARRS1 sets against the three artifacts that had passed as candidates.
let
    fig = Figure(size=(800, 330))
    ax = Axis(fig[1, 1]; xlabel="background excess around the detection (σ)", yticks=(1:3, ["artifacts", "junk from a bad set", "real objects"]),
              title="What max_background_excess separates", limits=((-0.5, 3.0), (0.4, 3.6)))
    rangebars!(ax, [3], [-0.07], [0.12]; direction=:x, color=:mediumspringgreen, linewidth=6, whiskerwidth=14)
    rangebars!(ax, [2], [-0.2], [0.2]; direction=:x, color=:gray, linewidth=6, whiskerwidth=14)
    scatter!(ax, [1.72, 2.44], [1, 1]; color=:tomato, markersize=16)
    text!(ax, 1.72, 1.25; text="ghost", align=(:center, :bottom), color=:white, fontsize=12)
    text!(ax, 2.44, 1.25; text="diffraction spike", align=(:center, :bottom), color=:white, fontsize=12)
    vlines!(ax, [1.0]; color=:white, linestyle=:dash)
    text!(ax, 1.05, 3.4; text="cut: 1.0σ", color=:white, fontsize=12)
    save(joinpath(ASSETS_DIR, "background-excess.png"), fig)
end

# --- Astrometry against an independent reduction --------------------------------------
# An unknown object in a real Pan-STARRS1 set (XY54_p10, 2025-09-16), measured
# by the pipeline and, independently, in Astrometrica: pipeline minus
# Astrometrica, per frame (RA offsets as arcsec on the sky).
let
    dra = [-0.008, -0.007, -0.003, -0.009] .* 15 .* cosd(-10.73)   # seconds of RA -> arcsec
    ddec = [0.15, -0.10, -0.12, -0.05]
    fig = Figure(size=(460, 440))
    ax = Axis(fig[1, 1]; xlabel="ΔRA cos δ (\")", ylabel="ΔDec (\")", aspect=DataAspect(),
              limits=((-0.3, 0.3), (-0.3, 0.3)), title="Pipeline − Astrometrica, four frames")
    θ = range(0, 2π; length=200)
    lines!(ax, 0.2 .* cos.(θ), 0.2 .* sin.(θ); color=(:white, 0.4), linestyle=:dash)
    text!(ax, 0.0, 0.21; text="0.2\"", align=(:center, :bottom), color=(:white, 0.7), fontsize=11)
    scatter!(ax, dra, ddec; color=:dodgerblue, markersize=14)
    for (k, (x, y)) in enumerate(zip(dra, ddec))
        text!(ax, x + 0.015, y; text=string(k), color=:white, fontsize=12)
    end
    hlines!(ax, [0.0]; color=(:white, 0.2)); vlines!(ax, [0.0]; color=(:white, 0.2))
    save(joinpath(ASSETS_DIR, "astrometry-vs-astrometrica.png"), fig)
end

# --- From 5,548 tracklets to vet to 125 -------------------------------------------------
# Tracklets a person would have to vet, and known objects recovered, on the
# real Pan-STARRS1 sets: the first configuration (raw frames, threshold 8,
# header WCS) and the one in pipeline_config, on six sets; then on ten.
let
    labels = ["first configuration\n(6 sets)", "difference images\n(6 sets)", "+ artifact filters\n(10 sets)"]
    tracklets = [5548, 157, 125]
    fig = Figure(size=(760, 380))
    ax = Axis(fig[1, 1]; ylabel="tracklets to vet", yscale=log10, xticks=(1:3, labels),
              title="Tracklets a person has to look at", limits=(nothing, (10, 20000)))
    barplot!(ax, 1:3, tracklets; color=[:tomato, :orange, :mediumspringgreen], fillto=10,
             bar_labels=string.(tracklets), label_color=:white)
    save(joinpath(ASSETS_DIR, "tracklets-to-vet.png"), fig)
end

# --- Is a known object nearby the same one? ------------------------------------------------
# The geometry nearby_known_objects tests: a known object predicted near a
# candidate counts as that object only if it moves at the same rate in the
# same direction AND is displaced along the track (an orbit's error is mostly
# in time). Schematic: four positions 15 minutes apart, at 100"/h.
let
    fig = Figure(size=(1000, 360))
    xs = (0:3) .* 25.0
    for (col, (title, dx, dy, color, verdict)) in enumerate((
            ("Displaced along the track", 60.0, 0.0, :mediumspringgreen, "same object, imprecise orbit"),
            ("Displaced across the track", 0.0, 40.0, :tomato, "another object")))
        ax = Axis(fig[1, col]; title, aspect=DataAspect(), limits=((-20, 170), (-35, 70)),
                  xlabel="arcsec east-west", ylabel=col == 1 ? "arcsec north-south" : "")
        lines!(ax, xs, zeros(4); color=:yellow, linewidth=3)
        scatter!(ax, xs, zeros(4); color=:yellow, markersize=13)
        lines!(ax, xs .+ dx, fill(dy, 4); color, linewidth=2, linestyle=:dash)
        scatter!(ax, xs .+ dx, fill(dy, 4); color=:transparent, strokecolor=color, strokewidth=2, markersize=15)
        text!(ax, -15, -22; text="candidate (measured)", color=:yellow, fontsize=12)
        text!(ax, -15, 55; text="known object (predicted): $verdict", color, fontsize=12)
    end
    save(joinpath(ASSETS_DIR, "nearby-geometry.png"), fig)
end

println("Figures written to $ASSETS_DIR")
