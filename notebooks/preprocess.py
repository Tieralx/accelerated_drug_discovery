import pandas as pd
import numpy as np
from pfeature import pcp_wp
from pfeature import ctd_wp
from pfeature import rri_wp
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re

def get_seqs (file_path):
    records_cd_hit = list(SeqIO.parse(file_path, "fasta"))
    seqs = []
    for records in records_cd_hit:
        seqs.append(re.sub(r'[Seq\(\'][\']','',str(records.seq)))
        
    return seqs

def pcp(input):
  a = input.rstrip('txt')
  output = a + 'dpc.csv'
  df_out = pcp_wp(input, output)
  df_in = pd.read_csv(output)
  return df_in

def ctd(input):
  a = input.rstrip('txt')
  output = a + 'ctd.csv'
  df_out = ctd_wp(input, output)
  df_in = pd.read_csv(output)
  return df_in

def rri(input):
  a = input.rstrip('txt')
  output = a + 'rri.csv'
  df_out = rri_wp(input, output)
  df_in = pd.read_csv(output)
  return df_in

def aplicar_pca(data, threshold = 0.90):
    print("Aplicando operación: PCA")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(data)
    # Aplicar PCA
    pca = PCA().fit(X_scaled)
    # Calcular la varianza explicada acumulada
    explained_variance = np.cumsum(pca.explained_variance_ratio_)
    num_components = 18 #np.argmax(explained_variance >= threshold) + 1
    pca = PCA(n_components=num_components)  # Número de componentes principales
    X_pca = pca.fit_transform(X_scaled)
    pca_columns = [f'PC{i+1}' for i in range(X_pca.shape[1])]  # Nombrar los componentes como PC1, PC2, etc.
    df_pca = pd.DataFrame(X_pca, columns=pca_columns)
    return df_pca

def escalado(data):
  # Crear una instancia del StandardScaler
  scaler = StandardScaler()
  # Ajustar y transformar los datos
  scaled_data = scaler.fit_transform(data)
  # Convertir los datos escalados a un DataFrame
  df_scaled = pd.DataFrame(scaled_data, columns=data.columns)
  # Mostrar los datos escalados
  return  df_scaled

def pfeature_process(cd_hit_path, file_path):
    feature_pcp = pcp(cd_hit_path)
    feature_ctd = ctd(cd_hit_path)
    feature_rri = rri(cd_hit_path)
    df_pfeatures = feature_pcp.merge(feature_ctd, how='inner', left_index=True, right_index=True)
    df_pfeatures = df_pfeatures.merge(feature_rri, how='inner', left_index=True, right_index=True)
    df_pfeatures = escalado(df_pfeatures)
    seqs = get_seqs(cd_hit_path) 
    df_pfeatures.index = seqs
    df_pfeatures.to_csv(file_path, index=True, index_label='Sequence')
    
    return df_pfeatures
    