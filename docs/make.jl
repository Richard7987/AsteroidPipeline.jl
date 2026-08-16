using AsteroidPipeline
using Documenter
using DocumenterVitepress

# README.md and INVESTIGATION_LOG.md are the actual source of truth for
# the project's purpose, status, and real-data findings (see
# INVESTIGATION_LOG.md itself for why they're kept separate). Copied here
# at build time, not committed under docs/src (see .gitignore), so the
# docs site can never drift out of sync with them.
const REPO_ROOT = joinpath(@__DIR__, "..")

index_content = read(joinpath(REPO_ROOT, "README.md"), String)
index_content = replace(index_content, "(INVESTIGATION_LOG.md)" => "(investigation-log.md)")
write(joinpath(@__DIR__, "src", "index.md"), index_content)

cp(joinpath(REPO_ROOT, "INVESTIGATION_LOG.md"), joinpath(@__DIR__, "src", "investigation-log.md"); force=true)

DocMeta.setdocmeta!(AsteroidPipeline, :DocTestSetup, :(using AsteroidPipeline); recursive=true)

makedocs(;
    modules=[AsteroidPipeline],
    authors="Alejandro",
    sitename="AsteroidPipeline.jl",
    # Needed explicitly (not just inside `format=`): this repo has no git
    # remote configured yet, and Documenter's own "edit this page" source
    # links can't be auto-detected from one that doesn't exist.
    repo="github.com/ale-bnes/AsteroidPipeline.jl.git",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/ale-bnes/AsteroidPipeline.jl",
        devbranch="main",
        devurl="dev",
    ),
    pages=[
        "Home" => "index.md",
        "Investigation Log" => "investigation-log.md",
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
    repo="github.com/ale-bnes/AsteroidPipeline.jl",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
