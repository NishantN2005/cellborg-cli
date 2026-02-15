import scanpy as sc
import pandas as pd
import json
import os

principle_components=50
resolution_global = 0.5
project_path = None
sc.settings.autoshow = False

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


def clustering(project_path, adata, resolution):
    global resolution_global
# ## Clustering
# 
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) {cite}`traag2019louvain`. Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    print("------- clustering begins --------")
    sc.tl.leiden(adata, resolution = resolution, flavor="igraph", n_iterations=-1)
    #sc.pl.umap(adata, color=["leiden"], save="2.png")

    # **** creating JSON for clustering ******
    print("test X_umap")
    print(adata.obsm)
    umap_df = pd.DataFrame(
    adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
    )
    umap_df["cluster"] = adata.obs["leiden"]


    umap_dict = umap_df.to_dict(orient="index")

    if not os.path.exists(f"{project_path}/cellborg-cli/figures/umap"):
        os.makedirs(f"{project_path}/cellborg-cli/figures/umap")
    umap_path = f"{project_path}/cellborg-cli/figures/umap/umap_clustering_res_{int(resolution*100)}.json"

    print('creating umap json file...')
    with open(umap_path, "w") as f:
        json.dump(umap_dict, f)


    resolution_global = resolution
    return (
        adata.obs["leiden"].cat.categories
    )

def cell_type_annotation(adata, annotations):
    adata.obs["cell_type_lvl1"] = adata.obs["leiden"].map(annotations)

def init_project(selected_project_path, adata):
    sc.settings.figdir = f"{selected_project_path}/cellborg-cli/figures"
    global project_path
    project_path = selected_project_path
    print("------ init_project begins -------")

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

def do_clustering(adata, resolution):
    global resolution_global
    global project_path
    try:
        clusters = clustering(project_path, adata, resolution)

        #add clustering resolution
        with open(f"{project_path}/cellborg-cli/project_values.json", "r+") as f:
            data = json.load(f)
            data['clust_resolution'] = resolution
            #write over file
            f.seek(0)
            json.dump(data, f)
            f.truncate()
        print("Added clustering resolution to project values")
        

        clusters = clusters.to_list()
        print(f"cluster: {type(clusters)}", clusters)

    except Exception as err:
        print('ERROR: ',str(err))


def gene_expression(adata, gene_list):
    try:
        print('----Starting Gene Expression----')
        gene_mask = adata.var_names.isin(gene_list)
        selected_genes = adata.var_names[gene_mask]
        expression_df = pd.DataFrame(
        adata[:, gene_mask].X.toarray(),  # Convert to dense matrix if needed
        columns=selected_genes,
        index=adata.obs.index
        )
        expression_df["cluster"] = adata.obs["leiden"]
        umap_df = pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
        )
    
        combined_df = pd.concat([umap_df, expression_df], axis=1)
        combined_dict = combined_df.to_dict(orient="index")

        print("saving gene expression to local...")
        gene_expression_path = f"{project_path}/cellborg-cli/figures/gene_expression_per_cell_with_clusters.json"
        with open(gene_expression_path, "w") as f:
            json.dump(combined_dict, f)
        
        print("Saved gene expression data to JSON file.")
    except Exception as err:
        print('ERROR: ',str(err))

def annotations(adata, annotations_dict):
    global project_path
    try:
        print('----Starting Annotations----')
        cell_type_annotation(adata, annotations_dict)
        print("creating annotations json...")
        adata.obs["cell_type_lvl1"] = adata.obs["leiden"].map(annotations_dict)
        annotations_df = pd.DataFrame(adata.obs["cell_type_lvl1"])
        annotations_df["cluster"] = adata.obs["leiden"]
        umap_df = pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
        )
    
        combined_df = pd.concat([umap_df, annotations_df], axis=1)
        combined_dict = combined_df.to_dict(orient="index")

        
        annotations_path = f"{project_path}/cellborg-cli/figures/umap/umap_annotations.json"
        print("saving annotations json...")
        with open(annotations_path, "w") as f:
            json.dump(combined_dict, f)
        
        print("Saved annotations to JSON file.")
    except Exception as err:
        print('ERROR: ',str(err))