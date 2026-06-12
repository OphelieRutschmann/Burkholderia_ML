

def extract_gene_mutations(feature_meta, geno_df):
"""See in which gene each sample has mutations
requires some sort of data cleaning before performing this step (i.e. variants that are very often -1 might need to be removed because they are low confidence.)"""
    meta_df = pd.DataFrame(feature_meta).T # rows: var_id. columns: (gene, impact, effect)

    # filter gene_df to keep only columns (variants) that have a gene annotation
    gene_df = geno_df.loc[:, meta_df["gene"].notna()]  

    # group gene_df by gene (works because meta_df and gene_df have the same indices)
    # take the max (the variants are grouped by gene. If any of the variants for a specific gene has a 1 (i.e. there is a mutation in that gene), return 1)
    gene_df = gene_df.groupby(meta_df["gene"].dropna(), axis=1).max()

    # transform -1 (allele absent) to 
    gene_df = gene_df.clip(lower=0)

    gene_df.columns = [f"snp_gene__{g}" for g in gene_df.columns]