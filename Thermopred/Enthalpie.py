# Needs Scikit-learn==1.3.2

import os
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem

ABSOLUT_PATH = os.path.dirname(os.path.realpath(__file__))


class EnthalpieEnergy():
    
    def predict(self, smiles:str):
        models = self.__loadModels()
        fingerprint = self.__prepareData(smiles=smiles)
        
        # Value predict by each model
        modelsIndividualContrib = {modelname: model.predict(fingerprint)[0] for modelname,model in models.items()}
        
        finalValue = np.array(list(modelsIndividualContrib.values())).mean()
        
        finalValue = self.__backReescale2Original(value=finalValue)
        
        return finalValue

    def __backReescale2Original(self,value:float):
        maxdata = 49.89374206197381
        mindata = 1.0
        denormalizedValue = (value*(maxdata - mindata) + mindata)

        originalValue = -56.474209*denormalizedValue

        return originalValue

    def __prepareData(self, smiles:str):
        mol = Chem.MolFromSmiles(smiles)
        fp = np.array([AllChem.GetMorganFingerprintAsBitVect(mol,radius=2,nBits=4096)])
        return fp

    def __loadModels(self):
        modelspath = ABSOLUT_PATH+'/models_enthalpie'
        pathfilenames = [modelspath+f'/{name}' for name in os.listdir(modelspath)]
        pathfilenames = [modelpath for modelpath in pathfilenames if modelpath.endswith('.sav')]
        

        # The names (Keys) are: 'XGBRegressor', 'MLPRegressor', 'RandomForestRegressor'
        modelsLoaded = {type(pickle.load(open(path,'rb'))).__name__ : pickle.load(open(path,'rb')) for path in pathfilenames}
        return modelsLoaded

if __name__ == '__main__':
    ee = EnthalpieEnergy()
    result = ee.predict(smiles='CN1C=CN(CCCN(c2cc(Cl)ccc2O)c2ccccc2S)CC1')
    print(result)