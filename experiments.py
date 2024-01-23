import pandas as pd
import lazypredict as lz
import rdkit.Chem as Chem

df = pd.read_csv('covid_submissions_all_info.csv')
df = df.dropna(subset=['r_avg_IC50', 'f_avg_IC50'], how='all')

df['r_avg_IC50'].isna().sum() # this is larger -> go with fluorescence for now
df['f_avg_IC50'].isna().sum()

# plot scatterplot of r_avg_IC50 vs f_avg_IC50