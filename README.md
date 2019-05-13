# SIMPLEGEN - instructions for contributors
Bob Verity <r.verity@imperial.ac.uk>


## Package structure

The overall package structure follows the standard structure of any R package, with the following tweaks:

- The "R_ignore" folder is setup to be ignored by R, but not by Github. This is where we can store files that are not part of a standard R package, but that we want included anyway, e.g. papers.
- If you make an "ignore" folder in your local version of the package then this will be ignored by *both* R and github. This is where you can store local files that you may want to be inside the package for organisational reasons, but not make it into anyone elses version, e.g. test scripts.
- Note, for both of the above folders, even though R claims to ignore them, in fact it copies and then deletes them when the package is built. Hence, avoid putting anything very large (e.g. whole genomes) in these folders as it will slow things down massively.
- Data can be stored in many places in an R package. For simplicity, here all data will be stored in inst/extdata. Objects in this folder will be loaded into R using a designated data loading function. This gives greater flexibility with file types - for example, data files can be stored in .csv making them easier to work with directly.


## Code structure and style

Contributors are asked to stick to the existing coding style, which follows Hadley Wickham's suggested style (http://r-pkgs.had.co.nz/r.html) with the following tweaks:

- single-line if statements, e.g. if (x < 10) y <- 5, are never OK. Always use {} over multiple lines.

Code should be structured and written in a way that is easy to pick up at a later date with minimal effort. Ease of reading takes precendence over code efficiency in most cases, as coding time is usually more costly than running time.

All R functions should be documented using the roxygen method (http://r-pkgs.had.co.nz/man.html). Both R and C++ functions should be documented clearly throughout.


## Github and branching

This package uses git and Github for version control. This has the advantage that it is impossible to completely break the code as we always have "savegames" going back through time.

We will use the git-flow branching pattern described here (https://nvie.com/posts/a-successful-git-branching-model/), please read this before working on the code. In simple terms:

- the master branch should always be the most recent stable release. Do not work directly on this branch.
- main development work should be on the develop branch.
- specific features should be on feature branches, but don't be shy of contributing to develop if that is more appropriate.
- once development is complete for a given release, we will create a release branch and at the same time merge changes back to master. After this point, the only changes should be bug-fixes.
- version numbers will be used to keep track of releases using the three-part X.Y.Z format, where X = major change, Y = small change e.g. added feature, Z = patch/bug fix. The development version will be 0.Y.Z, the first non-development release will be 1.0.0.


## Continuous integration

This package is setup for testing and continuous integration with both travis and appveyor. If you look at the main Github page (https://github.com/mrc-ide/SIMPLEGEN) you should see badges for both platforms showing that the current build is passing.

*If you make a change to the code that causes one of these builds to fail, please fix the error before continuing*


## C++

This package uses C++ through the Rcpp package. Unfortunately, however, debugging and profiling tools are not yet up to scratch for Rcpp. So, where possible within C++ code, hash-defines may be used to allow us to switch between different versions of the code that compile via different methods - for example within Rcpp vs. directly in Xcode on Mac.
