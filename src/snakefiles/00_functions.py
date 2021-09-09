"""functions.py: Accessory functions"""

def get_df_column(filepath, col_name, prefix, extension) :
    df = pd.read_csv(filepath)
    
    row = []
    for i in range(len(df)):
        suffix = df.loc[i, col_name]
        row += [prefix + suffix + extension]
        
    return(row)
