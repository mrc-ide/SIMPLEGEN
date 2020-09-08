# SIMPLEGEN - instructions for contributors
Bob Verity <r.verity@imperial.ac.uk>


## Package structure

The overall package structure follows the standard structure of any R package, with the following idiosyncrasies:

- The "R_ignore" folder is setup to be ignored by R, but not by Github. This is where we can store files that are not part of a standard R package, but that we want included anyway, e.g. papers.
- If you make an "ignore" folder in your local version of the package then this will be ignored by *both* R and github. This is where you can store local files that you may want to be inside the package for organisational reasons, but not make it into anyone elses version, e.g. test scripts.
- Note, for both of the above folders, even though R claims to ignore them, in fact it copies and then deletes them when the package is built. Hence, avoid putting anything very large (e.g. whole genomes) in these folders as it will slow things down massively.


## Code structure and style

Contributors are asked to stick to the existing coding style, which follows Hadley Wickham's suggested style (http://r-pkgs.had.co.nz/r.html) with the following tweaks:

- single-line if statements, e.g. if (x < 10) y <- 5, are never OK. Always use {} over multiple lines.

Code should be structured and written in a way that is easy to pick up at a later date with minimal effort. Ease of reading takes precendence over code efficiency in most cases, as coding time is usually more costly than running time.

All R functions should be documented using the roxygen method (http://r-pkgs.had.co.nz/man.html). Both R and C++ functions should be documented clearly throughout.


## Continuous integration

This package is setup for testing and continuous integration with both travis and appveyor. If you look at the main Github page (https://github.com/mrc-ide/SIMPLEGEN) you should see badges for both platforms showing that the current build is passing. It is crucial that certain branches pass all checks at all times (see below for details).


## Github and branching

This package uses git and Github for version control. This has the advantage that it is impossible to completely break the code as we always have saves going back through time.

We will use the git-flow branching pattern described here (https://nvie.com/posts/a-successful-git-branching-model/), please read this before working on the code. In simple terms:

- the master branch is the *outward-facing stable branch*. It should always represent the most recent official release, and therefore should never break. You should never work directly on this branch.
- the develop branch is the *inward-facing stable branch*. Similar to master, this branch should always pass all checks. The difference is that develop will continually move forward as new features are merged in, whereas master is frozen in time at official release points. You should never work directly on this branch.
- all code changes should be done through feature branches. These can be very small and short-lived, or more extensive. You are free to break and fix code as much as you like on feature branches. Once you are happy, these should be merged into develop using a Github pull request - not directly in the console. As long as the PR passes all checks you are free to complete the merge, or you can nominate a reviewer if you would like someone else to look over your changes. Always going through PRs in this way ensures that develop remains stable. We will periodically delete old feature branches once merged.
- version numbers will be used to keep track of releases using the three-part X.Y.Z format, where X = major change, Y = small change e.g. added feature, Z = patch/bug fix. The development version will be 0.Y.Z, the first non-development release will be 1.0.0.


## C++

This package uses C++ through the Rcpp package. Unfortunately, however, debugging and profiling tools are not yet up to scratch for Rcpp. So, where possible within C++ code, hash-defines may be used to allow us to switch between different versions of the code that compile via different methods - for example within Rcpp vs. directly in Xcode on Mac.
