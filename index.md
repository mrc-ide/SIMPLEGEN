[![Travis build status](https://travis-ci.org/mrc-ide/SIMPLEGEN.svg?branch=develop)](https://travis-ci.org/mrc-ide/SIMPLEGEN)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/SIMPLEGEN?branch=develop&svg=true)](https://ci.appveyor.com/project/mrc-ide/SIMPLEGEN)
[![Coverage status](https://codecov.io/gh/mrc-ide/SIMPLEGEN/branch/develop/graph/badge.svg)](https://codecov.io/github/mrc-ide/SIMPLEGEN?branch=develop)

# SIMPLEGEN

SIMPLEGEN is an R package for **sim**ulating **Pl**asmodium **e**pidemiological and **gen**etic data. The rationalle behind SIMPLEGEN is that computationally intensive forwards-simulation of *P.falciparum* genetic data can be avoided by splitting the process into multiple stages. These stages make up the SIMPLEGEN pipeline, shown below:

<img src="https://raw.githubusercontent.com/mrc-ide/SIMPLEGEN/master/R_ignore/images/pipeline.png" height="185px" width="800px" />

Notice that genotypes are not actually generated until near the end of the pipeline, at which point we have already discarded many unimportant events. For example, we never need to model the genetics of infections that eventually clear, or infections in hosts that are not connected to our final sample. By discarding unimportant events like this we get significant gains in speed, which in turn means we can model large spatial scales. For example, it would not be unusual to simulate a population of thousands of hosts distributed across hundreds of connected demes, from which we can sample whole *P.falciparum* genomes in a matter of seconds.

A second advantage of breaking simulation into stages is that users are free to enter the pipeline at any stage. For example, an entirely different epidemiological model could be used in the first stage, just as long as it can produce a transmission record in a standardised format. Simiilarly, if users wanted to keep the in-built epidemiological model but explore a different model of parasite genetics then they could enter at a later stage.

To get started, take a look at the [installation instructions](https://mrc-ide.github.io/SIMPLEGEN/articles/installation.html), followed by a [basic tutorial](https://mrc-ide.github.io/SIMPLEGEN/articles/basic_tutorial.html) on running the program.



