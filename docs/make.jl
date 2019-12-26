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
<<<<<<< Updated upstream
    repo = "Klafyvel/Wavefronts.jl",
=======
    repo = "github.com/Klafyvel/Wavefronts.jl.git",
>>>>>>> Stashed changes
)
