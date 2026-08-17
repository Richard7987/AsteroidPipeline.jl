using AsteroidPipeline
using Documenter
using DocumenterVitepress

# index.md, investigation-log.md, variable-star-validation.md,
# iasc-campaign-validation.md, and design-refinements.md are all real,
# permanent, committed pages under docs/src — the docs site's own
# content, not generated or copied in from README.md/INVESTIGATION_LOG.md
# (which no longer exist at the repo root; the short root README.md just
# links here instead).
DocMeta.setdocmeta!(AsteroidPipeline, :DocTestSetup, :(using AsteroidPipeline); recursive=true)

# Regenerates the figures investigation-log.md embeds, from the real
# numbers measured during this project's real-data investigations — run
# before makedocs so both local builds and CI (julia-actions/julia-docdeploy
# instantiates docs/Project.toml, which includes CairoMakie) produce them
# fresh rather than relying on committed images that could drift from the
# numbers in the text.
include("make_figures.jl")

makedocs(;
    modules=[AsteroidPipeline],
    authors="Alejandro",
    sitename="AsteroidPipeline.jl",
    # Needed explicitly (not just inside `format=`): the actual git remote
    # (Forgejo, self-hosted) isn't github.com, so Documenter's own "edit
    # this page" source-link auto-detection can't work here — this must
    # match the *GitHub mirror* deploydocs below actually deploys to.
    repo="github.com/Richard7987/AsteroidPipeline.jl.git",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/Richard7987/AsteroidPipeline.jl",
        devbranch="main",
        devurl="dev",
    ),
    pages=[
        "Home" => "index.md",
        "Investigation Log" => [
            "ZTF field 451 (the original investigation)" => "investigation-log.md",
            "Variable-star validation" => "variable-star-validation.md",
            "IASC / Pan-STARRS1 campaign validation" => "iasc-campaign-validation.md",
            "Design refinements" => "design-refinements.md",
        ],
        "API Reference" => [
            "Detection & Linking" => "api/detection.md",
            "Variable Stars" => "api/variables.md",
            "Astrometry & Plate-Solving" => "api/astrometry.md",
            "Cross-Matching" => "api/crossmatch.md",
            "Reference & ZOGY Differencing" => "api/reference-zogy.md",
            "Pipeline" => "api/pipeline.md",
            "Rotation Period" => "api/rotation.md",
        ],
    ],
)

DocumenterVitepress.deploydocs(;
    repo="github.com/Richard7987/AsteroidPipeline.jl",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
