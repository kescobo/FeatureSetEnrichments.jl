using FeatureSetEnrichments
using Documenter

DocMeta.setdocmeta!(FeatureSetEnrichments, :DocTestSetup, :(using FeatureSetEnrichments); recursive=true)

makedocs(;
    modules=[FeatureSetEnrichments],
    authors="Kevin Bonham, PhD <kbonham@wellesley.edu> and contributors",
    repo="https://github.com/kescobo/FeatureSetEnrichments.jl/blob/{commit}{path}#{line}",
    sitename="FeatureSetEnrichments.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kescobo.github.io/FeatureSetEnrichments.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kescobo/FeatureSetEnrichments.jl",
    devbranch="main",
)
