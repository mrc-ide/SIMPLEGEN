[![Travis build status](https://travis-ci.org/mrc-ide/SIMPLEGEN.svg?branch=master)](https://travis-ci.org/mrc-ide/SIMPLEGEN)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/SIMPLEGEN?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/SIMPLEGEN)
[![Coverage status](https://codecov.io/gh/mrc-ide/SIMPLEGEN/branch/master/graph/badge.svg)](https://codecov.io/github/mrc-ide/SIMPLEGEN?branch=master)

# SIMPLEGEN
Simulating Plasmodium Epidemiological and Genetic Data

Forwards-in time simulation of Plasmodium falciparum genetic data can be
computationally intensive, as many genotypes are tracked but ultimately lost.
SIMPLEGEN avoids this problem by splitting simulation into three steps - first,
simulating transmission under a simple individual-based model, second, pruning
the infection history to focus on nodes that contribute to the final sample, and
third, simulating genetic data. A secondary advantage is that any third-party
epidemiological simulator can be used for the first step as long as it outputs
in the correct format. The major limitation of SIMPLEGEN is that assumes the
epidemiological and genetic processes are seperable, and hence is limited to
neutral variation, and cannot model selection.

**This package is currently in development. Code contributors must read the
guidelines in the R_ignore/CODER_README.txt file before contributing**.
