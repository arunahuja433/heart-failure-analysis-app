
import pandas as pd

def analyze_expression(gene_df, hopkins_df):
    risk_genes = set(gene_df['gene_id'].str.upper())
    hopkins_df['geneid'] = hopkins_df['geneid'].str.upper()

    hopkins_df['padjpef'] = pd.to_numeric(hopkins_df['padjpef'], errors='coerce')
    hopkins_df['padjref'] = pd.to_numeric(hopkins_df['padjref'], errors='coerce')

    sig_pef = hopkins_df[((hopkins_df['padjpef'] < 0.05) | (hopkins_df['padjpef'] == 0.0)) & 
                         (hopkins_df['geneid'].isin(risk_genes))].copy()
    sig_pef['Direction'] = sig_pef['l2fcpef'].apply(lambda x: 'Up' if x > 0 else 'Down')

    sig_ref = hopkins_df[((hopkins_df['padjref'] < 0.05) | (hopkins_df['padjref'] == 0.0)) & 
                         (hopkins_df['geneid'].isin(risk_genes))].copy()
    sig_ref['Direction'] = sig_ref['l2fcref'].apply(lambda x: 'Up' if x > 0 else 'Down')

    geneids_pef = set(sig_pef['geneid'])
    geneids_ref = set(sig_ref['geneid'])

    exclusive_pef_ids = geneids_pef - geneids_ref
    exclusive_ref_ids = geneids_ref - geneids_pef

    sig_pef['Group'] = sig_pef['geneid'].apply(lambda x: 'Exclusive to HFpEF' if x in exclusive_pef_ids else 'Shared')
    sig_ref['Group'] = sig_ref['geneid'].apply(lambda x: 'Exclusive to HFrEF' if x in exclusive_ref_ids else 'Shared')

    return sig_pef, sig_ref
