# SIMPLEGEN - instructions for contributors
### Bob Verity <r.verity@imperial.ac.uk>

## Package structure

The overall package structure follows the standard structure of any R package, with the following tweaks

- the "R_ignore" folder is setup to be ignored by R, but not by Github. This is where we can store files that are not part of a standard R package, but that we want included anyway, e.g. papers
- if you make an "ignore" folder in your local version of the package then this will be ignored by *both* R and github. This is where you can store local files that you may want to be inside the package for organisational reasons, but not make it into anyone elses version, e.g. test scripts
- Note, for both of the above folders, even though R claims to ignore them, in fact it copies and then deletes them when the package is built. Hence, avoid putting anything very large (e.g. whole genomes) in these folders as it will slow things down massively

