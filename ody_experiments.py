import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# import lazypredict as lz

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def compute_descriptors(molecule):
    descriptors = {d[0]: d[1](molecule) for d in Descriptors.descList} 
    descriptors = pd.Series(descriptors) 
    return descriptors

def read_and_clean_data(path):
    df = pd.read_csv(path)
    df = df.dropna(subset=['r_avg_IC50', 'f_avg_IC50'], how='all')
    df = df.drop_duplicates(subset=['SMILES'], keep='first')
    df = df[df['f_avg_IC50']<=100] # removes extreme outliers
    dfs = df[['SMILES','CID','f_avg_IC50','r_avg_IC50']]
    return df, dfs

df, dfs = read_and_clean_data('raw_df_notest.csv')

# get all molecular descriptors
descriptors = dfs['SMILES'].apply(lambda smiles: compute_descriptors(Chem.MolFromSmiles(smiles)))
descriptors['CID'] = dfs['CID'].tolist()
# get Morgan fingerprints
dfs['FP'] = dfs['SMILES'].apply(lambda smiles: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 3, 2048))
dfs['FP'] = dfs['FP'].apply(lambda fp: list(fp))












# # UMAP plot of Morgan fingerprints
# import umap
# reducer = umap.UMAP()
# embedding = reducer.fit_transform(dfs['FP'].tolist())
# plt.scatter(embedding[:, 0], embedding[:, 1], c=dfs['f_avg_IC50'], cmap='Spectral', s=5)
# plt.colorbar()
# plt.title('UMAP projection of Morgan fingerprints')
# plt.show()





# from molflux.datasets import featurise_dataset
# from molflux.features import load_from_dicts as load_representations_from_dicts

# featuriser = load_representations_from_dicts(
#     [
#         {"name": "morgan"},
#         {"name": "maccs_rdkit"},
#     ]
# )

# featurised_dataset = featurise_dataset(dfs, column="SMILES", representations=featuriser)

# print(featurised_dataset)



# X = dfs['FP'].to_numpy()
# y = dfs['f_avg_IC50'].to_numpy()
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# lz.Supervised.LazyRegressor(verbose=0, ignore_warnings=True, custom_metric=None).fit(X_train, X_test, y_train, y_test)

# lz.models.test_models(X_train, X_test, y_train, y_test)

# # regress umap embeddings against IC50 values
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LinearRegression
# from sklearn.metrics import mean_squared_error, r2_score
# from sklearn.preprocessing import StandardScsaler

# X = embedding
# y = dfs['f_avg_IC50'].tolist()
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
# scaler = StandardScaler()
# X_train = scaler.fit_transform(X_train)
# X_test = scaler.transform(X_test)
# reg = LinearRegression().fit(X_train, y_train)
# y_pred = reg.predict(X_test)
# print('MSE: ', mean_squared_error(y_test, y_pred))

