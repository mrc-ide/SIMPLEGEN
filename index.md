[![Travis build status](https://travis-ci.org/mrc-ide/SIMPLEGEN.svg?branch=develop)](https://travis-ci.org/mrc-ide/SIMPLEGEN)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/SIMPLEGEN?branch=develop&svg=true)](https://ci.appveyor.com/project/mrc-ide/SIMPLEGEN)
[![Coverage status](https://codecov.io/gh/mrc-ide/SIMPLEGEN/branch/develop/graph/badge.svg)](https://codecov.io/github/mrc-ide/SIMPLEGEN?branch=develop)

# SIMPLEGEN

SIMPLEGEN is an R package for *sim*ulating *Pl*asmodium *e*pidemiological and
*gen*etic data. The rationalle behind it is that naive forwards simulation of
*P.falciparum* genetic data tends to be very computationally expensive, as many
genotypes are tracked but ultimately lost when the host clears infection.
SIMPLEGEN avoids this problem by splitting the process of simulation into
multiple stages.

<br/>
<br/>
<img src="https://raw.githubusercontent.com/mrc-ide/SIMPLEGEN/master/R_ignore/images/pipeline.png" height="93px" width="300px" />
<br/>

At each stage only the minimum amount of information is stored
and passed forwards, which greatly reduces the computational complexity of the
problem.

Key features of SIMPLEGEN are:

* Item1


