
import pandas as pd

def identify_independent_loci(df, window_size=500_000):
    """
    Identifies independent risk loci ±window_size apart per chromosome.
    Returns a filtered DataFrame of independent loci.
    """
    df = df.sort_values(by=['chromosome', 'p_value'])
    distinct_loci = []
    selected_positions = {}

    for _, row in df.iterrows():
        chrom = str(row['chromosome'])
        pos = row['base_pair_location']

        if chrom not in selected_positions:
            selected_positions[chrom] = []

        too_close = any(abs(pos - prev_pos) <= window_size for prev_pos in selected_positions[chrom])

        if not too_close:
            distinct_loci.append(row)
            selected_positions[chrom].append(pos)

    return pd.DataFrame(distinct_loci)
