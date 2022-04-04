## BINtools v0.1.1 (Release date: 2022-04-04)
* Moved rstantools from "suggests" to "imports"
* Updated all packages
* Received an error: "Files in the 'vignettes' directory but no files in 'inst/doc'". The solution was to run pkgdown::build_site(), delete /docs directory, and move the /doc directory to /inst. 

## BINtools v0.1.0 (Release date: 2021-06-08)
* This is the initial submission
