using Documenter, Katana

makedocs(
    modules = [Katana],
    format = :html,
    sitename = "Katana",
    pages = [
        "Home" => "index.md",
        "Manual" => [],
        "Library" => "library.md",
        "Developer" => []
    ]
)

deploydocs(
    deps = nothing,
    make = nothing,
    target = "build",
    repo = "github.com/lanl-ansi/Katana.jl.git",
    julia = "0.6"
)
