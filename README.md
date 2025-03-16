# Thermopred
This repository contains the official data, algorithms, and ML models present in the paper "`AI-Enhanced Quantum Chemistry Dataset for Thermochemical Properties of API-Like Compounds and Their Degradants`".

## How to use

Download the repository manually or via git:

```shell
$ git clone https://github.com/jeffrichardchemistry/thermopred
```

Enter the `thermopred` directory and run the following command to install the python package:

```shell
$ cd thermopred

$ python3 setup.py install
```

Once this is done, it is now possible to use the package by simply importing the modules. Import the modules as described below and pass a smiles for prediction.

```python
from Thermopred.Enthalpie import EnthalpieEnergy
from Thermopred.GibbsEnergy import GibbsFreeEnergy

smiles='CN1C=CN(CCCN(c2cc(Cl)ccc2O)c2ccccc2S)CC1'

ee = EnthalpieEnergy()
result_enthalpie = ee.predict(smiles)

gfe = GibbsFreeEnergy()
gfe.predict(smiles=smiles)
```
