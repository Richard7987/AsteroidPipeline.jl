#=
Generates the figures embedded in docs/src/investigation-log.md, from the
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

println("Figures written to $ASSETS_DIR")
