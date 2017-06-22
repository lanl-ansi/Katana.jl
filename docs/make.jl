using Documenter, Katana

makedocs(
    modules = [Katana],
    format = :html,
    sitename = "Katana",
    pages = [
        "Home" => "index.md",
        "Manual" => [],
        "Library" => [],
        "Developer" => []
    ]
)

