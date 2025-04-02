import streamlit as st
import pandas as pd
import requests

# File path to your filtered GWAS .csv file (already filtered for significant SNPs)
gwas_file_path = "/Users/Arun/Desktop/Northwestern/Research/ShahCardio/hf_gwas_filtered.csv"
# Load the file to check the columns
df = pd.read_csv("/Users/Arun/Desktop/Northwestern/Research/ShahCardio/hf_gwas_filtered.csv")
print(df.columns)  # This will print the column names

# Ensembl API for fetching gene location
ENSEMBL_API = "https://rest.ensembl.org/lookup/symbol/human/"

def get_gene_location(gene_id_or_name):
    """ Fetch gene location (chromosome and base pair location) from Ensembl API. """
    if gene_id_or_name.startswith('ENSG'):
        # Gene ID format
        url = f"https://rest.ensembl.org/lookup/id/{gene_id_or_name}?content-type=application/json"
    else:
        # Gene Symbol format
        url = f"https://rest.ensembl.org/lookup/symbol/human/{gene_id_or_name}?content-type=application/json"

    response = requests.get(url)

    # Print the response for debugging purposes
    if response.status_code != 200:
        st.error(f"Error: {response.status_code} - {response.text}")
    else:
        # Only print to the console for debugging, not Streamlit UI
        print(f"API Response for {gene_id_or_name}: {response.json()}")

    if response.status_code == 200 and response.json():
        results = response.json()
        if 'seq_region_name' in results and 'start' in results:
            chromosome = results['seq_region_name']
            start_position = results['start']
            return chromosome, start_position

    return None, None

def filter_gwas_variants(chromosome, base_pair_location, pval_threshold=5e-8, window_size=500_000):
    """ Filter GWAS variants based on location and p-value from the filtered CSV file """
    chromosome = int(chromosome)
    try:
        # Load the filtered GWAS CSV file instead of the large .dta file
        df = pd.read_csv(gwas_file_path)
        df_filtered = df[(df['p_value'] <= pval_threshold) &
                         (df['chromosome'] == chromosome) &
                         (df['base_pair_location'] >= base_pair_location - window_size) &
                         (df['base_pair_location'] <= base_pair_location + window_size)]
        return df_filtered
    except Exception as e:
        st.error(f"An error occurred while filtering variants: {e}")
        return pd.DataFrame()  # Return an empty DataFrame in case of error

def get_gene_id_from_symbol(gene_symbol):
    """Fetch gene ID from gene symbol using the Ensembl API."""
    url = f"{ENSEMBL_API}lookup/symbol/human/{gene_symbol}?content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200 and response.json():
        results = response.json()
        if 'id' in results:
            return results['id']
    return None

def check_gene_association_in_hopkins(gene_input):
    """ Check if gene is associated with HFpEF or HFrEF in the Hopkins dataset """
    hopkins_file_path = "/Users/Arun/Desktop/Northwestern/Research/ShahCardio/hopkins.dta"
    hopkins_df = pd.read_stata(hopkins_file_path)
    hopkins_df['geneid'] = hopkins_df['geneid'].str.upper()

    hopkins_df['padjpef'] = pd.to_numeric(hopkins_df['padjpef'], errors='coerce')
    hopkins_df['padjref'] = pd.to_numeric(hopkins_df['padjref'], errors='coerce')

    # If input is gene symbol, convert it to gene ID
    gene_id = gene_input
    if not gene_id.startswith('ENSG'):  # If it's not a gene ID
        gene_id = get_gene_id_from_symbol(gene_input)

    if gene_id is None:
        st.error(f"Gene ID for {gene_input} not found.")
        return pd.DataFrame(), pd.DataFrame()  # Return empty DataFrames

    sig_pef = hopkins_df[((hopkins_df['padjpef'] < 0.05) | (hopkins_df['padjpef'] == 0.0)) &
                         (hopkins_df['geneid'].str.upper() == gene_input.upper())]
    sig_ref = hopkins_df[((hopkins_df['padjref'] < 0.05) | (hopkins_df['padjref'] == 0.0)) &
                         (hopkins_df['geneid'].str.upper() == gene_input.upper())]

    return sig_pef, sig_ref


# Handling the user input and response
if __name__ == "__main__":
    gene_input = st.text_input("Enter Gene ID or Gene Name")

    if gene_input:
        chromosome, base_pair_location = get_gene_location(gene_input)

        if chromosome and base_pair_location:
            st.write(f"Gene {gene_input} is located on chromosome {chromosome} at base pair {base_pair_location}.")

            # Checking association with HFpEF or HFrEF
            sig_pef, sig_ref = check_gene_association_in_hopkins(gene_input)

            if not sig_pef.empty:
                st.write(f"Gene {gene_input} is associated with HFpEF (p-adjusted < 0.05).")
                st.dataframe(sig_pef)
            if not sig_ref.empty:
                st.write(f"Gene {gene_input} is associated with HFrEF (p-adjusted < 0.05).")
                st.dataframe(sig_ref)

            # Filtering GWAS variants based on location
            st.write(f"Searching for variants within 500 kb of {gene_input} on chromosome {chromosome}...")
            filtered_variants = filter_gwas_variants(chromosome, base_pair_location)

            if not filtered_variants.empty:
                st.write(f"Found {len(filtered_variants)} variants near {gene_input}.")
                st.dataframe(filtered_variants[['variant_id', 'chromosome', 'base_pair_location', 'p_value']])

                # Allow user to download filtered variants
                st.download_button(
                    label="Download Filtered Variants",
                    data=filtered_variants.to_csv(index=False),
                    file_name=f"{gene_input}_filtered_variants.csv",
                    mime="text/csv"
                )

            else:
                st.write(f"No variants found within the specified window.")

     