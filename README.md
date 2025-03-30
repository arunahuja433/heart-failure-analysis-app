# Heart Failure Analysis App

This is a Streamlit-based web application designed for analyzing heart failure genetic data. The app allows users to upload GWAS `.dta` files and provides several analytical tools such as filtering significant variants, identifying independent loci, annotating genes, and comparing gene expression between Heart Failure with Preserved Ejection Fraction (HFpEF) and Heart Failure with Reduced Ejection Fraction (HFrEF).

## Features
- **GWAS Filtering**: Filters genome-wide significant variants based on p-value threshold.
- **Independent Loci Identification**: Identifies independent loci within ±500kb of significant variants.
- **Gene Annotation**: Annotates loci with the nearest genes using Ensembl data.
- **HFpEF vs HFrEF Comparison**: Compares gene expression data from two types of heart failure (HFpEF vs. HFrEF).

## Technologies
- **Streamlit**: For building the interactive web app.
- **Pandas**: For handling data processing.
- **PyArrow**: For efficiently handling large datasets.
- **Ensembl API**: For annotating genes.
- **Python**: The programming language used for the app.

## Installation

### Prerequisites:
- Python 3.7+
- Git (optional, for cloning the repo)

### Steps:
1. Clone this repository:
    ```bash
    git clone https://github.com/arunahuja433/heart-failure-analysis-app.git
    cd heart-failure-analysis-app
    ```

2. Create a virtual environment and activate it:
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # For MacOS/Linux
    venv\Scripts\activate  # For Windows
    ```

3. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

4. Run the app:
    ```bash
    streamlit run StreamLit.py
    ```

The app will open in your default web browser, and you can begin uploading `.dta` files for analysis.

## How to Use
1. Upload a `.dta` file or provide the full file path (for large files).
2. The app will filter significant GWAS variants, identify independent loci, and annotate genes.
3. You can then compare gene expression data between HFpEF and HFrEF groups.
4. Results can be downloaded as CSV files for further analysis.

## Contributing
1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -am 'Add new feature'`).
4. Push to the branch (`git push origin feature-branch`).
5. Open a pull request with a description of the changes.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

