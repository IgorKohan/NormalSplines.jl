using Documenter, NormalSplines

makedocs(
    sitename = "NormalSplines.jl",
	format = Documenter.HTML(),
	authors = "Igor Kohanovsky",	
    pages = [
				"Home" => "index.md",
				"Example Usage" => "Usage.md",
				"Public API" => "Public-API.md",
			]
)

deploydocs(
	     repo = "github.com/IgorKohan/NormalSplines.jl.git",
	    )
