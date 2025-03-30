
import requests

ENSEMBL_API = "https://rest.ensembl.org/overlap/region/human/"

def get_nearest_gene(chromosome, position):
    url = f"{ENSEMBL_API}{chromosome}:{position}-{position}?feature=gene;content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200 and response.json():
        results = response.json()
        gene_ids = [gene.get('gene_id', '') for gene in results]
        gene_names = [gene.get('external_name', '') for gene in results]
        return gene_ids[0], gene_names[0] if gene_ids else ("No gene found", "No name found")
    return "No gene found", "No name found"

def annotate_genes(df):
    genes = [get_nearest_gene(row['chromosome'], row['base_pair_location']) for _, row in df.iterrows()]
    gene_ids, gene_names = zip(*genes)
    df['gene_id'] = gene_ids
    df['external_name'] = gene_names
    return df
