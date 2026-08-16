using AsteroidPipeline
using Documenter
using DocumenterVitepress

# README.md is the actual source of truth for the project's purpose and
# status. Copied here at build time, not committed under docs/src (see
# .gitignore), so the site's home page can never drift out of sync with
# it. INVESTIGATION_LOG.md lives permanently at docs/src/investigation-log.md
# instead (a real, committed page, not generated) — it's the docs site's
# own content, not duplicated from anywhere else.
const REPO_ROOT = joinpath(@__DIR__, "..")

index_content = read(joinpath(REPO_ROOT, "README.md"), String)
index_content = replace(index_content, "(docs/src/investigation-log.md)" => "(investigation-log.md)")
write(joinpath(@__DIR__, "src", "index.md"), index_content)

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
