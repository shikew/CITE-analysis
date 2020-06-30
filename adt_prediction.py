import joblib
import pandas as pd
import numpy as np

def rf_multi():
    cbmc_rna = pd.read_csv('F://PythonCode/graduation project/cbmc_rna_HVGs_200.csv', header=0, index_col=0)
    cbmc_rna = pd.DataFrame(cbmc_rna)
    columns_adt = cbmc_rna.columns
    cbmc_rna = cbmc_rna.T

    model = joblib.load('F://PythonCode/graduation project/save/rf_multi.pkl')
    yhat = model.predict(cbmc_rna)
    data_adt = yhat.T
    data_adt = pd.DataFrame(data_adt)
    data_adt.index = ["CD3","CD4","CD8","CD45RA","CD56","CD16","CD11c","CD14","CD19","CD34"]
    data_adt.columns = columns_adt
    return data_adt

data_adt = rf_multi()
