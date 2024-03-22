import pandas as pd
from scipy.sparse import coo_matrix

def dataframe_to_adjacency_matrix(df, row_column_name='gene_name',
                                  col_column_name='trait_reported',
                                  values_col_name='L2G',
                                  threshold=None,
                                  sparse_output=False):
    """
    Convert a DataFrame into an adjacency matrix.

    Parameters:
        df (pandas.DataFrame): Input DataFrame.
        row_column_name (str): Name of the column to use as row indices.
        col_column_name (str): Name of the column to use as column indices.
        values_col_name (str): Name of the column to use for values in the matrix.
        threshold (float or None): Threshold value for creating a binary matrix. If provided,
            values greater than the threshold will be set to 1, and values less than or equal
            to the threshold will be set to 0. Default is None.
        sparse_output (bool): If True, return a sparse matrix (scipy.sparse.coo_matrix),
            otherwise return a dense matrix (pandas.DataFrame). Default is False.

    Returns:
        scipy.sparse.coo_matrix or pandas.DataFrame: An adjacency matrix representation of the input DataFrame.
            If sparse_output is True, a sparse matrix is returned; otherwise, a dense matrix is returned.

    """
    # Select columns from the DataFrame
    selected_cols = [row_column_name, col_column_name, values_col_name]
    selected_df = df[selected_cols]
    
    # Pivot the DataFrame
    pivot_df = selected_df.pivot_table(index=row_column_name,
                                       columns=col_column_name,
                                       values=values_col_name,
                                       aggfunc='first',
                                       fill_value=0)
    
    if threshold is not None:
        # Threshold the values to create a binary matrix
        binary_matrix = pivot_df.apply(lambda x: (x > threshold).astype(int))
        
        if sparse_output:
            # Convert the binary matrix to a sparse matrix
            adjacency_matrix = coo_matrix(binary_matrix.values)
        else:
            adjacency_matrix = binary_matrix
    else:
        if sparse_output:
            # Convert the pivot table DataFrame to a sparse matrix
            adjacency_matrix = coo_matrix(pivot_df.values)
        else:
            adjacency_matrix = pivot_df
    
    return adjacency_matrix