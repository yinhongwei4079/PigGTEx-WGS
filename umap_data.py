import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import requests
import seaborn as sns
import umap
import sys

from MulticoreTSNE import MulticoreTSNE as mTSNE
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

df = pd.read_csv(sys.argv[1],index_col = 0, sep='\t')
print("done")
df.index = df_index

def reduce_dim(X, algorithm='PCA', n_components=3):
    if algorithm == 'PCA':
        pca = PCA(n_components=n_components).fit(X)
        X_red = pca.transform(X)
        # merge your data into the same table as the reference data
        df_red = pd.DataFrame(X_red, 
                              index=df_index, 
                              columns=['component1', 'component2', 'component3'])
    elif algorithm == 'TSNE':
        # TSNE, Barnes-Hut have dim <= 3
        if n_components > 3:
            print('The Barnes-Hut method requires the dimensionaility to be <= 3')
            return None
        else:
            X_red = mTSNE(n_components=n_components, n_jobs=4).fit_transform(X)
            df_red = pd.DataFrame(X_red, 
                                  index=df_index, 
                                  columns=['component1', 'component2', 'component3'])
    elif algorithm == 'UMAP':
        umap_ = umap.UMAP(n_components=n_components).fit(X)
        X_red = umap_.transform(X)
        # merge your data into the same table as the reference data
        df_red = pd.DataFrame(X_red, 
                              index=df_index, 
                              columns=['component1', 'component2', 'component3'])
    else:
        return None
    
    return df_red

for i in range(3,20):
    for j in ['PCA','TSNE','UMAP']:
        df_umap = reduce_dim(df.values,algorithm=j,n_components=i)
        df_umap.to_csv("/vol3/agis/likui_group/yinhongwei/gte_reference/06_prun/pca/GTE.1602." + str(i) + "." + j + "." + ".csv")
