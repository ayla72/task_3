import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define the Alzheimer’s-related genes of interest
genes_of_interest = [
    "APP", "PSEN1", "MAPT", "SLC2A1",  # AD and energy metabolism
    "GFAP", "HSPA1A", "HSP90AA1",      # Stress response
    "SNAP25", "SYN1", "GRIN1"          # Synaptic health
]

# Load the dataset into an AnnData object
adata = sc.read_h5ad("Al.h5ad")

# Inspect adata.var to locate the column for gene names
print("Inspecting the first few rows of adata.var to locate gene names:")
print(adata.var.head())

# Use the 'feature_name' column to filter for Alzheimer’s-related genes if available
if 'feature_name' in adata.var.columns:
    adata_filtered = adata[:, adata.var['feature_name'].isin(genes_of_interest)].copy()
else:
    adata_filtered = adata[:, [gene for gene in genes_of_interest if gene in adata.var_names]].copy()

# Split data by disease status for AD and normal subjects
adata_AD = adata_filtered[adata_filtered.obs['disease'] == 'Alzheimer disease']
adata_normal = adata_filtered[adata_filtered.obs['disease'] == 'normal']

# Calculate mean expression for genes of interest in AD and normal cohorts
mean_expression_AD = adata_AD.to_df().mean()
mean_expression_normal = adata_normal.to_df().mean()

# Combine results into a DataFrame for easy comparison
expression_comparison = pd.DataFrame({
    'Gene': adata_filtered.var['feature_name'].values,
    'Mean Expression (AD)': mean_expression_AD.values,
    'Mean Expression (Normal)': mean_expression_normal.values
})

# Save mean expression comparison to CSV
expression_comparison.to_csv("mean_expression_comparison_AD_vs_Normal.csv", index=False)
print("Mean expression comparison saved to 'mean_expression_comparison_AD_vs_Normal.csv'")

# Plot mean expression levels of genes in AD vs. normal aging
plt.figure(figsize=(10, 6))
sns.barplot(x="Gene", y="value", hue="variable",
            data=pd.melt(expression_comparison, id_vars=["Gene"], value_vars=["Mean Expression (AD)", "Mean Expression (Normal)"]))
plt.title("Mean Gene Expression Levels in AD vs. Normal Aging")
plt.ylabel("Mean Expression")
plt.xlabel("Genes")
plt.xticks(rotation=45)
plt.legend(title="Condition")
plt.savefig("mean_expression_AD_vs_Normal.png", dpi=300, bbox_inches="tight")  # Save the mean expression plot
plt.show()

# Additional: Calculate variance for genes in AD and normal groups
variance_expression_AD = adata_AD.to_df().var()
variance_expression_normal = adata_normal.to_df().var()

variance_comparison = pd.DataFrame({
    'Gene': adata_filtered.var['feature_name'].values,
    'Variance (AD)': variance_expression_AD.values,
    'Variance (Normal)': variance_expression_normal.values
})

# Save variance comparison to CSV
variance_comparison.to_csv("variance_expression_comparison_AD_vs_Normal.csv", index=False)
print("Variance expression comparison saved to 'variance_expression_comparison_AD_vs_Normal.csv'")

# Print variance comparison
print("\nVariance in Gene Expression (AD vs. Normal):")
print(variance_comparison)

# Visualize variance in gene expression in AD vs. normal aging
plt.figure(figsize=(10, 6))
sns.barplot(x="Gene", y="value", hue="variable",
            data=pd.melt(variance_comparison, id_vars=["Gene"], value_vars=["Variance (AD)", "Variance (Normal)"]))
plt.title("Variance in Gene Expression Levels in AD vs. Normal Aging")
plt.ylabel("Variance in Expression")
plt.xlabel("Genes")
plt.xticks(rotation=45)
plt.legend(title="Condition")
plt.savefig("variance_expression_AD_vs_Normal.png", dpi=300, bbox_inches="tight")  # Save the variance plot
plt.show()