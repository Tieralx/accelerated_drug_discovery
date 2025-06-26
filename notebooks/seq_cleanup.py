import pandas as pd
import re

def quitar_x_inicio(secuencia):
    return re.sub(r'^X+', '', secuencia)

def quitar_x_final(secuencia):
    return re.sub(r'X+$', '', secuencia)

def clean_seq(df):
    aa_pfeature = ['A','C','D','E','F','G','H','I','K','L','M','N','P',
                   'Q','R','S','T','V','W','Y']

    # 1. Validar que existe la columna 'sequences'
    if 'sequences' not in df.columns:
        raise KeyError("La columna esperada 'sequences' no se encuentra en el DataFrame.")

    # 2. Poner en mayúsculas todas las secuencias
    df['sequences'] = df['sequences'].str.upper()

    # 3. Dividir entre secuencias con y sin X
    aa_x = df[df['sequences'].str.contains('X')].copy()
    no_x = df[~df['sequences'].str.contains('X')].copy()

    # 4. Limpiar X al inicio y final
    aa_x['sequences'] = aa_x['sequences'].apply(quitar_x_inicio)
    aa_x['sequences'] = aa_x['sequences'].apply(quitar_x_final)

    # 5. Filtrar secuencias válidas
    valid_chars = lambda seq: all(c in aa_pfeature for c in seq)
    aa_x = aa_x[aa_x['sequences'].apply(valid_chars)]
    no_x = no_x[no_x['sequences'].apply(valid_chars)]

    # 6. Combinar todo
    filtered_df = pd.concat([no_x, aa_x], ignore_index=True)

    return filtered_df
