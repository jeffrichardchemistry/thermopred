[![DOI](https://img.shields.io/badge/DOI-10.1021/acs.jcim.5c01320-orange)](https://doi.org/10.1021/acs.jcim.5c01320)

# Thermopred
This repository contains the official data, algorithms, and ML models present in the paper "`AI-Enhanced Quantum Chemistry Dataset for Thermochemical Properties of API-Like Compounds and Their Degradants`".

## We recommend using Python 3.10, as this version has been thoroughly tested with Thermopred.

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

You can also install this package directly from PyPI using:
```shell
pip install thermopred
```

# How to cite
- LaTeX
```
@article{doi:10.1021/acs.jcim.5c01320,
author = {Santos, Diullio P. and Dias-Silva, Jefferson R. and Júnior, Luiz H. K. Q. and de Oliveira, Heibbe C. B.},
title = {ThermoPred: AI-Enhanced Quantum Chemistry Data Set and ML Toolkit for Thermochemical Properties of API-Like Compounds and Their Degradants},
journal = {Journal of Chemical Information and Modeling},
volume = {0},
number = {0},
pages = {null},
year = {0},
doi = {10.1021/acs.jcim.5c01320},
note ={PMID: 41334815}
}
```

- Text
```
Santos, D. P.; Dias-Silva, J. R.; Júnior, L. H. K. Q.; de Oliveira, H. C. B. ThermoPred: AI-Enhanced Quantum Chemistry Data Set and ML Toolkit for Thermochemical Properties of API-Like Compounds and Their Degradants. J. Chem. Inf. Model. 2025
```

