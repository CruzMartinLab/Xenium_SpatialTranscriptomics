import scanpy as sc
import squidpy as sq
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from spatialdata_io import xenium
#Set seed for reproducability
import random
import torch
np.random.seed(123)
random.seed(123)
torch.manual_seed(123)

#Load Xenium Data
xenium_path = "C:\\Users\\cmneu\\OneDrive\\Desktop\\Xenium_Seurat\\MT_sag"
sdata = xenium(xenium_path)
adata = sdata.tables["table"]

#QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
adata.obs["neg_probe_pct"] = (
    adata.obs["control_probe_counts"] / adata.obs["total_counts"] * 100
)
adata.obs["neg_decoding_pct"] = (
    adata.obs["control_codeword_counts"] / adata.obs["total_counts"] * 100
)

#Plot QC Metrics
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

sns.histplot(adata.obs["total_counts"], kde=True, ax=axs[0])
axs[0].set_title("Total Transcripts per Cell")
axs[0].set_xlabel("Total Counts")

sns.histplot(adata.obs["n_genes_by_counts"], kde=True, ax=axs[1])
axs[1].set_title("Unique Genes per Cell")
axs[1].set_xlabel("Genes Detected")

sns.histplot(adata.obs["neg_probe_pct"], kde=True, ax=axs[2])
axs[2].set_title("Negative Probe %")
axs[2].set_xlabel("Percentage")

plt.tight_layout()
plt.show()

#Cell Filtering
#Remove cells with:
# - Very low total transcripts (bottom 5%)
# - Very high total transcripts (top 99%) (doublet removal)
# - Very low genes detected
# - High mitochondrial fraction (suggests dying cells)

#Mitochondrial Fraction
mito_genes = adata.var_names.str.startswith("MT-")
mito_counts = adata[:, mito_genes].X.sum(axis=1)
if hasattr(mito_counts, "A1"):
    mito_counts = mito_counts.A1
else:
    mito_counts = np.asarray(mito_counts).flatten()
adata.obs["mito_frac"] = mito_counts / adata.obs["total_counts"]

#Visualization of Mitochrondrial Dist.
adata.obs["mito_frac"].hist(bins=50)
plt.xlabel("Mitochondrial Fraction")
plt.ylabel("Number of Cells")
plt.title("Mitochondrial Content Distribution")
plt.show()

#Actual Filter
adata = adata[adata.obs["mito_frac"] < 0.1].copy()

#Transcript Count
low_counts = np.percentile(adata.obs["total_counts"], 5)
high_counts = np.percentile(adata.obs["total_counts"], 99)

sc.pp.filter_cells(adata, min_counts=low_counts)
sc.pp.filter_cells(adata, max_counts=high_counts)

# Check minimum and maximum counts after filtering
print(f"Min count after filtering: {np.min(adata.obs['total_counts'])}")
print(f"Max count after filtering: {np.max(adata.obs['total_counts'])}")

#Low Gene Detection
low_genes = np.percentile(adata.obs["n_genes_by_counts"], 5)
sc.pp.filter_genes(adata, min_counts=low_genes)

#Gene Filtering
# - Remove genes detected in fewer than threshold cells
# - Remove genes with extremely low expression

#Gene expression visualization
gene_counts = (adata.X > 0).sum(axis=0).A1
plt.hist(gene_counts, bins=100)
plt.xlabel("Number of Cells a Gene is Expressed In")
plt.ylabel("Number of Genes")
plt.title("Gene Detection Across Cells")
plt.show()

#Gene filter
min_cells = int(0.005 * adata.n_obs) # approx 735 cells
sc.pp.filter_genes(adata, min_cells = min_cells)
sc.pp.filter_genes(adata, min_counts = 100)

#Detect and Remove Doublets
import scrublet as scr

scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs["doublet_scores"] = doublet_scores
adata.obs["predicted_doublets"] = predicted_doublets

#Remove detected doublets
adata = adata[~adata.obs["predicted_doublets"], :]

#Normalization
adata.layers["counts"] = adata.X.copy() #Preserve raw counts
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)#highly variable gene selection
adata = adata[:, adata.var.highly_variable]
adata.X = adata.X + 1e-6 #Deal with zero counts that don't actually exist
sc.pp.normalize_total(adata, target_sum = 1e4)

#PCA, Neighbors, Clustering
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver="arpack") #computes pca with eigenvalues - good for sparse matrix
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["doublet_scores", "mito_frac"])
sc.tl.leiden(adata)

#Spatial Analysis
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
sq.gr.centrality_scores(adata, cluster_key="leiden")

# Plot spatial clusters
sc.pl.spatial(adata, color="leiden", spot_size=20)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "leiden"], wspace=0.4)
##FOR NORMAL LEIDEN CLUSTERING UNCOMMENT ALL LINES WITH ## AND COMMENT OUT EVERYTHING BELOW THIS

#Trying SpatialLeiden Clustering
import spatialleiden as sl
seed = 123
sc.pp.log1p(adata)
sc.pp.pca(adata, random_state=seed)
sc.pp.neighbors(adata, random_state=seed)

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs = 10)

adata.obsp["spatial_connectivities"] = sl.distance2connectivity(
    adata.obsp["spatial_distances"]
)
sc.tl.leiden(adata, directed = False, random_state=seed)
sl.spatialleiden(adata, layer_ratio=1.8, directed=(False, True), seed=seed)
sc.pl.embedding(adata, basis="spatial", color=["leiden", "spatialleiden"])

## Change to tissue name before running this line + wherever you want to save it
adata.write("C:/Users/cmneu/OneDrive/Desktop/Xenium_Seurat/MT_sag/MT_sag_QC_filtered.h5ad")
