# Superflow

## Pre-Requisites

**Superflow** is a package that analyzes susyNt so it is assumed that the **SusyNtuple** package is already checked out in your working directory. 

**Superflow** can also be configured for making fake ntuples and computing charge-flip probability ~out-of-the-box and so requiries the necessary packages in order to be compiled. These packages may need to be updated in order to run on ROOT6 by removing dependencies/includes on Cint. The necessary packages can be checked out by running the following command (once you have checked out Superflow):

`source Superflow/scripts/get_superflow_dep.sh`

At the moment this checks out two packages:

1. **DileptonMatrixMethod**

    Package written by Davide Gerbaudo (Davide.Gerbaudo@cern.ch) to calculate fake efficiencies, etc... located at [github.com/gerbaudo/DileptonMatrixMethod](https://github.com/gerbaudo/DileptonMatrixMethod).

2. **ChargeFlip**

    Package holding charge-flip probability maps and methods to calculate the charge-flip probablity. This package is located at [https://svnweb.cern.ch/cern/wsvn/atlasphys/Physics/SUSY/Analyses/WeakProduction/ChargeFlip/](https://svnweb.cern.ch/cern/wsvn/atlasphys/Physics/SUSY/Analyses/WeakProduction/ChargeFlip/)
