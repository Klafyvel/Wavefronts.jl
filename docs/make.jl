using Documenter, Wavefronts

makedocs(
    sitename="Wavefronts.jl documentation",
    pages = [
        "index.md",
        "example.md"
    ],
    modules = [Wavefronts]
)

deploydocs(
    repo = "github.com/Klafyvel/Wavefronts.jl.git",
)
