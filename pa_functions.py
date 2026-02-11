import scanpy as sc
import json

principle_components=50

def read_adata(adata_path):
    print("------ read_adata begins -------")
    try:
        adata = sc.read_h5ad(adata_path)
        print(f"Successfully read adata from {adata_path}")
        return adata
    except Exception as e:
        raise RuntimeError(f"Failed to read adata file '{adata_path}': {e}")
    
def normalize(adata):
    print("-----normalize begins----")
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    print ("----normalize completed-----")

def feature_selection(adata):
    # ## Feature selection
    # 
    # As a next step, we want to reduce the dimensionality of the dataset and only include the most informative genes. This step is commonly known as feature selection. The scanpy function `pp.highly_variable_genes` annotates highly variable genes by reproducing the implementations of Seurat {cite}`Satija2015`, Cell Ranger {cite}`Zheng2017`, and Seurat v3 {cite}`stuart2019comprehensive` depending on the chosen `flavor`. 
    print("-------- feature selection begins --------")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata)

def dimentionality_reduction(adata):
    global principle_components
    # ## Dimensionality Reduction
    # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
    print("------- dimentionality reduction begins ------")
    sc.tl.pca(adata)

    #find ideal number of principle components
    df = adata.uns["pca"]['variance_ratio'].cumsum(axis=0)
    if df[-1] < 0.95:
        print("use 50")
    else:
        for num in range(len(df)):
            if df[num] >=0.95:
                principle_components = num+1
                break
                
    # Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function {func}`~scanpy.tl.leiden` or {func}`~scanpy.tl.tsne`. In our experience, there does not seem to be signifigant downside to overestimating the numer of principal components.
    sc.pl.pca_variance_ratio(adata, n_pcs=principle_components, log=True, save=".png")

    # You can also plot the principal components to see if there are any potentially undesired features (e.g. batch, QC metrics) driving signifigant variation in this dataset. In this case, there isn't anything too alarming, but it's a good idea to explore this.
    sc.pl.pca(
        adata,
        color=["pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (1, 2)],
        ncols=2,
        size=2,
        save=".png",
    )
    
    print("----- nearest_neighbor_graph begins -----")
    sc.pp.neighbors(adata, n_pcs=principle_components)
    # This graph can then be embedded in two dimensions for visualiztion with UMAP (McInnes et al., 2018):
    sc.tl.umap(adata)
    # We can now visualize the UMAP according to the `sample`. 
    sc.pl.umap(
        adata,
        # Setting a smaller point size to get prevent overlap
        size=2,
        save="1.png"
    )


def init_project(selected_project_path, adata):
    sc.settings.figdir = f"{selected_project_path}/cellborg-cli/figures"

    #TODO: concat datasets
    print("TODO: Concate datasets")
    normalize(adata)
    print('Successfully normalized concatonated dataset')
    feature_selection(adata)
    print('Successfully selected features')
    dimentionality_reduction(adata)
    print('Successfully reduced dimentionality')

    numpcs = {
            "num_pcs":principle_components,
            "gene_list": adata.var_names.to_list()
    }
    with open(f"{selected_project_path}/cellborg-cli/project_values.json", "w") as outputfile:
            json.dump(numpcs, outputfile)
    print("created project_values.json file")