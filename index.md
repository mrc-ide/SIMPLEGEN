
<img src="https://raw.githubusercontent.com/mrc-ide/SIMPLEGEN/master/R_ignore/images/simplegenlogo.png" height="123px" width="533px" />

## Watch this space...

**SIMPLEGEN is currently in development and is not yet ready for real analysis sorry!** We hope to have a working release version start of 2024.

SIMPLEGEN is an R package for **sim**ulating **Pl**asmodium **e**pidemiological and **gen**etic data. The rationalle behind SIMPLEGEN is that computationally intensive forwards-simulation of *P. falciparum* genetic data can be avoided by splitting the process into multiple stages. These stages make up the SIMPLEGEN pipeline:

<img src="https://raw.githubusercontent.com/mrc-ide/SIMPLEGEN/master/R_ignore/images/pipeline.png" height="185px" width="800px" />


Notice that genotypes are not actually generated until near the end of the pipeline, at which point we have already discarded many unimportant events. For example, we never need to model the genetics of infections that eventually clear, or infections in hosts that are not connected to our final sample. By discarding unimportant events like this we get significant gains in speed, which has a knock-on effect on the size and spatial scale that we can handle. For example, the final software is aiming to simulate a population of thousands of hosts distributed across hundreds of connected demes, and to sample whole *P. falciparum* genomes.

A second advantage of breaking simulation into stages is that the user is free to enter the pipeline at any stage. For example, an entirely different epidemiological model could be used in the first stage, just as long as it can produce a transmission record in the required format for the second stage of the pipeline. Similarly, if someone wanted to use the inbuilt epidemiological model but wanted to explore a different model of parasite genetics then they could enter at a later stage. 

To get started with SIMPLEGEN, take a look at the [installation instructions](https://mrc-ide.github.io/SIMPLEGEN/articles/installation.html), followed by a [basic tutorial](https://mrc-ide.github.io/SIMPLEGEN/articles/basic_tutorial.html) on how to run the program. Alternatively, you may want to read about the [basic epidemiological model](https://mrc-ide.github.io/SIMPLEGEN/articles/model_description.html) that comes with SIMPLEGEN.



