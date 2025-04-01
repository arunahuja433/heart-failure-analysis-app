import streamlit as st
import pandas as pd
import os
import time

# Import processing functions from modular scripts
from gwas_filter import filter_significant_variants
from loci_identifier import identify_independent_loci
from gene_annotator import annotate_genes
from expression_analyzer import analyze_expression
from gene_snp_analysis import get_gene_location, filter_gwas_variants, check_gene_association_in_hopkins

# Set Streamlit app configuration
st.set_page_config(page_title="Heart Failure Analysis App", layout="wide")

# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Select a page", ("Home", "Gene/SNP Analysis"))

# Main page
if page == "Home":
    # Title and initial instructions for Home page
    st.title("Heart Failure Analysis App")
    st.write("Upload your `.dta` file (≤ 2 GB) or enter the file path (> 2 GB) to begin.")

    # -----------------------------
    # FILE UPLOAD OR PATH INPUT
    # -----------------------------

    # Option 1: Upload a .dta file directly (for files under 2GB)
    uploaded_file = st.file_uploader("Upload a `.dta` file (≤ 2 GB)", type=["dta"])

    # Option 2: Enter full file path (for very large files already on server)
    file_path = st.text_input("Enter the full path to your `.dta` file (> 2 GB):", "")

    df = None

    # -----------------------------
    # LOAD FILE: Handle uploaded file
    # -----------------------------
    if uploaded_file is not None:
        try:
            with open(uploaded_file.name, "wb") as f:
                f.write(uploaded_file.getbuffer())
            df = pd.read_stata(uploaded_file.name)
            st.success("File uploaded and read successfully!")
            os.remove(uploaded_file.name)
        except Exception as e:
            st.error(f"Error reading uploaded file: {e}")
    elif file_path:
        if os.path.exists(file_path):
            df = pd.read_stata(file_path)
            st.success("File loaded successfully from path!")
        else:
            st.error("The specified file path does not exist. Please enter a valid path.")

    # -----------------------------
    # STEP 0: PREVIEW DATA
    # -----------------------------
    if df is not None:
        st.subheader("Data Preview")
        st.dataframe(df.head())  # Show the first few rows of the uploaded data

        # -----------------------------
        # STEP 1: FILTER GWAS HITS
        # -----------------------------
        st.subheader("Step 1: Filter Genome-Wide Significant Variants")
        if 'p_value' not in df.columns:
            st.error("The dataset does not contain a 'p_value' column.")
        else:
            # Filter for genome-wide significant hits (p ≤ 5×10⁻⁸)
            filtered_df = filter_significant_variants(df)
            st.write(f"✅ Total genome-wide significant variants: {len(filtered_df)}")
            st.dataframe(filtered_df.head())  # Show preview of filtered variants

            # Allow user to download the filtered GWAS results
            st.download_button(
                label="Download Filtered GWAS Results",
                data=filtered_df.to_csv(index=False),
                file_name="filtered_gwas.csv",
                mime="text/csv"
            )

        # -----------------------------
        # STEP 2: IDENTIFY INDEPENDENT LOCI
        # -----------------------------
        st.subheader("Step 2: Identify Independent Risk Loci")
        if filtered_df is not None and not filtered_df.empty:
            independent_loci_df = identify_independent_loci(filtered_df)
            st.write(f"✅ Total independent loci identified: {len(independent_loci_df)}")
            st.dataframe(independent_loci_df.head())  # Show a preview of independent loci data

            # Allow user to download the independent loci results
            st.download_button(
                label="Download Independent Loci Results",
                data=independent_loci_df.to_csv(index=False),
                file_name="independent_loci.csv",
                mime="text/csv"
            )
        else:
            st.warning("Please filter the significant variants first.")

        # -----------------------------
        # STEP 3: ANNOTATE GENES
        # -----------------------------
        st.subheader("Step 3: Annotate with Nearest Genes")
        if independent_loci_df is not None and not independent_loci_df.empty:
            annotated_df = annotate_genes(independent_loci_df)
            st.write(f"✅ Total genes annotated: {len(annotated_df)}")
            st.dataframe(annotated_df.head())  # Show a preview of the annotated data

            # Allow user to download the annotated gene results
            st.download_button(
                label="Download Annotated Gene Results",
                data=annotated_df.to_csv(index=False),
                file_name="annotated_genes.csv",
                mime="text/csv"
            )
        else:
            st.warning("Please identify independent loci first.")

        # File path for the Hopkins dataset
        hopkins_file_path = "/Users/Arun/Desktop/Northwestern/Research/ShahCardio/hopkins.dta"

        # Step 4: Compare HFpEF vs HFrEF Expression Data
        st.subheader("Step 4: Compare HFpEF vs HFrEF Expression Data")

        # Check if the Hopkins file exists
        if os.path.exists(hopkins_file_path):
            try:
                # Read the Hopkins dataset
                hopkins_df = pd.read_stata(hopkins_file_path)
                st.success("Hopkins data loaded successfully!")

                # Pass both the annotated dataset and Hopkins dataset to analyze_expression
                expression_comparison_df, additional_info = analyze_expression(annotated_df, hopkins_df)

                # Display the comparison results
                st.write("HFpEF vs HFrEF expression comparison:")
                st.dataframe(expression_comparison_df)

                # Optionally, allow the user to download the comparison data
                st.download_button(
                    label="Download HFpEF vs HFrEF Comparison",
                    data=expression_comparison_df.to_csv(index=False),
                    file_name="hfpef_vs_hfref_comparison.csv",
                    mime="text/csv"
                )
            except Exception as e:
                st.error(f"Error loading Hopkins data: {e}")
        else:
            st.error("The specified Hopkins file does not exist. Please check the path.")



elif page == "Gene/SNP Analysis":
    # Gene/SNP Analysis Page
    # Gene/SNP Analysis Page
    st.title("Gene/SNP Analysis")
    gene_input = st.text_input("Enter Gene ID or Gene Name", "")

    if gene_input:
        # Step 1: Get Gene Location (Chromosome, Base Pair Location)
        chromosome, base_pair_location = get_gene_location(gene_input)

        if chromosome and base_pair_location:
            st.write(f"Gene {gene_input} is located on chromosome {chromosome} at base pair {base_pair_location}.")

            # Step 2: Check for association in Hopkins dataset (HFpEF or HFrEF)
            sig_pef, sig_ref = check_gene_association_in_hopkins(gene_input)

            if not sig_pef.empty:
                st.write(f"Gene {gene_input} is associated with HFpEF (p-adjusted < 0.05).")
                st.dataframe(sig_pef)
            if not sig_ref.empty:
                st.write(f"Gene {gene_input} is associated with HFrEF (p-adjusted < 0.05).")
                st.dataframe(sig_ref)

            # Step 3: Filter variants based on gene location
            st.write(f"Searching for variants within 500 kb of {gene_input} on chromosome {chromosome}...")

            try:
                # Start the progress bar
                progress_bar = st.progress(0)
                # Optionally, add a spinner to show that it's processing
                with st.spinner("Filtering GWAS variants..."):
                    # Simulating processing time (remove in production)
                    time.sleep(1)  # Just to simulate time passing

                    filtered_variants = filter_gwas_variants(chromosome, base_pair_location)

                    # Update progress to 50%
                    progress_bar.progress(50)

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
                        st.write("No variants found within the specified window.")

                    # Complete the progress bar
                    progress_bar.progress(100)

            except Exception as e:
                # Catch any errors or timeouts that may occur
                st.error(f"An error occurred while filtering variants: {e}")

        else:
            st.error("Gene location not found. Please check the Gene ID or Name.")
