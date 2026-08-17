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
# Real measurements from sweeping systematic_error_fraction against (a)
# 152 real, matched, high-S/N stationary stars on real ZTF field 451
# (false-positive rate at chi2_threshold=10.0) and (b) the real,
# independently-confirmed variable ASASSN-V J183620.31's own real
# forced-photometry light curve (reduced chi2). See "Tuning
# find_variable_sources's systematic error floor..." in this log.
let
    floors = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]  # percent
    false_positive_rate = [13.2, 6.6, 3.3, 2.6, 2.0, 1.3, 0.0]  # percent, at chi2_threshold=10
    variable_reduced_chi2 = [20688.5, 484.4, 123.4, 55.1, 31.0, 13.8, 5.0]

    fig = Figure(size=(720, 460))
    ax1 = Axis(fig[1, 1];
        xlabel="systematic_error_fraction (%)",
        ylabel="false-positive rate at chi2_threshold=10 (%)",
        title="Real ZTF field 451 stationary stars (152) vs. a real confirmed variable")
    ax2 = Axis(fig[1, 1];
        ylabel="reduced chi² (real variable, ASASSN-V J183620.31)",
        yscale=log10, yaxisposition=:right, ygridvisible=false)
    hidespines!(ax2)
    hidexdecorations!(ax2)

    l1 = lines!(ax1, floors, false_positive_rate; color=:tomato, linewidth=2)
    scatter!(ax1, floors, false_positive_rate; color=:tomato, markersize=10)
    l2 = lines!(ax2, floors, variable_reduced_chi2; color=:dodgerblue, linewidth=2)
    scatter!(ax2, floors, variable_reduced_chi2; color=:dodgerblue, markersize=10)
    hlines!(ax2, [10.0]; color=:dodgerblue, linestyle=:dash, linewidth=1)
    vlines!(ax1, [2.0]; color=:white, linestyle=:dot, linewidth=1)
    text!(ax1, 2.05, 12.0; text="chosen default (2%)", fontsize=12, color=:white)

    Legend(fig[2, 1], [l1, l2],
           ["false-positive rate (left axis)",
            "real variable's reduced chi² (right axis, log scale; dashed = chi2_threshold=10)"];
           orientation=:horizontal, tellwidth=false)

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
