# Thermopred
This repository contains the official data, algorithms, and ML models present in the paper "`AI-Enhanced Quantum Chemistry Dataset for Thermochemical Properties of API-Like Compounds and Their Degradants`".

## How to use

Download the repository manually or via git:

```shell
$ git clone https://github.com/jeffrichardchemistry/thermopred
```
Enter the `thermopred` models directory and run the following command to unzip models:

```shell
$ cd ../thermopred/Thermopred/models
$ tar -xzvf MLP_gibbs.sav.tar.gz && tar -xzvf RF_gibbs.sav.tar.gz && tar -xzvf XGB_gibbs.sav.tar.gz
```

and
```shell
$ cd ../thermopred/Thermopred/models_entalpie
$ tar -xzvf MLP_gibbs_enthalpie.sav.tar.gz && tar -xzvf RF_gibbs_Enthalpie.sav.tar.gz && tar -xzvf XGB_gibbs_Enthalpie.sav.tar.gz
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
