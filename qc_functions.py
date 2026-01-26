import scanpy as sc
import anndata as ad
import pandas as pd
from pathlib import Path as _Path
from datetime import datetime
import json


MT_TO_SPECIES = {'MT-ND1': "Homo sapiens",
                "mt-Nd1":"Mus musculus", 
                "mt-nd1":"Danio rerio", 
                "mt:ND1":"Drosophila melanogaster", 
                "Mt-nd1":"Rattus Norvegicus"
                }

SPECIES_TO_MT = {
                "Homo sapiens":"MT-",
                "Mus musculus":"mt-",
                "Danio rerio":"mt-",
                "Rattus norvegicus":"Mt-",
                "Drosophila melanogaster":"mt:"
                }

# load dataset files create AnnData object
def read_10x_mtx(workspace_path):
    try:
        print("Trying scanpy 10x read")
        adata = sc.read_10x_mtx(
        workspace_path,  # the directory with the `.mtx` file
        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
        cache=False,  # write a cache file for faster subsequent reading
        )
        return adata
    except:
        print("Falling back to manual MTX read")
        adata = sc.read_mtx(f'{workspace_path}/matrix.mtx.gz')
        adata_bc=pd.read_csv(f'{workspace_path}/barcodes.tsv.gz',header=None)
        adata_features=pd.read_csv(f'{workspace_path}/features.tsv.gz',header=None, delimiter = '\t')
        adata= adata.T
        adata.obs['cell_id']= adata_bc
        adata.var['gene_name'] = adata_features[1].values
        adata.var.index= adata.var['gene_name']

        return adata
    
def calculate_qc_metrics(mt: str, adata: ad.AnnData):
    min_genes=200
    min_cells=3
    sc.pp.filter_cells(adata, min_genes)
    sc.pp.filter_genes(adata, min_cells)
    # The scanpy function {func}`~scanpy.pp.calculate_qc_metrics` calculates common quality control (QC) metrics, which are largely based on `calculateQCMetrics` from scater {cite}`McCarthy2017`. One can pass specific gene population to {func}`~scanpy.pp.calculate_qc_metrics` in order to calculate proportions of counts for these populations. Mitochondrial, ribosomal and hemoglobin genes are defined by distinct prefixes as listed below. 

    # mitochondrial genes, "MT-" for human, "mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith(mt)
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )


def find_species(features_path):
    """Read a `features.tsv.gz` file and return a pandas DataFrame.

    Parameters
    - features_path: path to the gzipped TSV (str or Path)

    Returns
    - DataFrame with columns [feature_id, gene_name, ...] depending on the file.
    """

    p = _Path(features_path)
    print("path to features file:", p)
    if not p.exists():
        raise FileNotFoundError(f"features file not found: {p}")

    # Many 10x feature files are gzipped TSVs with 2 or 3 columns.
    try:
        print("Trying to read features file")
        df = pd.read_csv(p, sep='\t', header=None, compression='infer')
        print('here is if it is there: ', 'MT-ND1' in df[1].values)
        
        for mt_gene, species in MT_TO_SPECIES.items():
            if mt_gene in df[1].values:
                print(f"Detected species: {species} (found {mt_gene})")
                return species
        
        raise ValueError("Mitochondrial gene not found in features file.")
    except Exception as e:
        raise RuntimeError(f"Failed to read features file '{p}': {e}")