
import pandas as pd

def filter_significant_variants(df, pval_threshold=5e-8):
    """
    Filters GWAS DataFrame for genome-wide significant variants.
    Returns the filtered DataFrame.
    """
    if 'p_value' not in df.columns:
        raise ValueError("DataFrame must contain a 'p_value' column.")

    filtered_df = df[df['p_value'] <= pval_threshold].copy()
    return filtered_df
