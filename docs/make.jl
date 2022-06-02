using Documenter
using iceFEM

pages = [
  "Home" => "index.md"
]

makedocs(
  sitename = "iceFEM.jl",
  format = Documenter.HTML(),
  modules = [iceFEM],
  pages = pages
)

deploydocs(
  repo = "https://github.com/Balaje/iceFEM.jl"
)
