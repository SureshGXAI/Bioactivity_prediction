import rdkit 
from rdkit import Chem 
from rdkit.Chem import Descriptors
import numpy as np
import warnings 
import os, sys

warnings.filterwarnings('ignore')

def LigandEfficiency(smile):
    mol = Chem.MolFromSmiles(smile)
    HA = Descriptors.HeavyAtomCount(mol)
    LogP = Descriptors.MolLogP(mol)
    LE = - BE / HA
    LEscale = 0.873 * np.exp(-0.026 * HA) - 0.064
    LELP = LogP / LE
    FQ = LE / LEscale 
    return LE, LEscale, LELP, FQ


def BindingConstant(BE):
    T = 298.15
    R = 1.987
    t = np.float128(-1* BE / R*T)
    ki = np.exp(t)
    return ki


def main(smile, BE):
    LigandEfficiency(smile)
    BindingConstant(BE)
    return 

if __name__ == "__main__":
    if len(sys.argv) == 3:
        smile = sys.argv[1]
        BE = sys.argv[2]
    else:
        print("Usage: python3 bioactivity_prediction.py smile BE", file=sys.stderr)


#print(BioactivePrediction('O=C(N[C@H](C(=O)N[C@H](B(O)O)CC(C)C)Cc1ccccc1)c2nccnc2', -9.45))
