[![Travis build status](https://travis-ci.org/mrc-ide/SIMPLEGEN.svg?branch=develop)](https://travis-ci.org/mrc-ide/SIMPLEGEN)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/SIMPLEGEN?branch=develop&svg=true)](https://ci.appveyor.com/project/mrc-ide/SIMPLEGEN)
[![Coverage status](https://codecov.io/gh/mrc-ide/SIMPLEGEN/branch/develop/graph/badge.svg)](https://codecov.io/github/mrc-ide/SIMPLEGEN?branch=develop)

# SIMPLEGEN
Simulating Plasmodium Epidemiological and Genetic Data

Forwards-in-time simulation of Plasmodium falciparum genetic data can be
computationally intensive, as many genotypes are tracked but ultimately lost
when the host clears infection. SIMPLEGEN avoids this problem by splitting
simulation into three steps - first, simulating transmission under a simple
individual-based model, second, pruning the transmission record down to focus on
hosts that contribute to the final sample, and third, simulating genetic data
from the pruned transmission record. The main advantage of this method is speed,
which can be several orders of magnitude faster than ordinary forwards-in-time
simulation. A secondary advantage is that any third-party epidemiological
simulator can be used for the first step as long as it outputs in the correct
format. The major limitation of SIMPLEGEN is that it assumes the epidemiological
and genetic processes are seperable, hence it cannot be used to model selection
(i.e. drug resistance).

*This package is currently in development. Code contributors **must** read the
guidelines [here](https://github.com/mrc-ide/SIMPLEGEN/tree/style_guide) before contributing.*

