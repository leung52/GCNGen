import pandas as pd


def generate_node_file(corr_file: str, go_terms_file: str, output_file: str) -> None:
	edges = pd.read_csv(corr_file, sep='\t', usecols=['GENE1', 'GENE2'])
	nodes = set(edges['GENE1']).union(set(edges['GENE2']))
	go_terms_df = pd.read_csv(go_terms_file)
	filtered_df = go_terms_df[go_terms_df['gene'].isin(nodes)]
	filtered_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
	corr = 'processed_data/filtered_tmu_genes_kendall_b_genes_exp_5_cells_p_005_corr_99_percentile.tsv'
	go = 'Tmu_scRNASeq_objects/Tmu_gene_annotations_interpro_GO_EC_Pfam.csv'
	out = 'cytoscape/tmu_gen_nodes_kendall_b_p_0.05_99th_percentile.tsv'

	generate_node_file(corr, go, out)

	edge = 'cytoscape/tmu_gen_edges_kendall_b_p_0.05_99th_percentile.tsv'
	corr_df = pd.read_csv(corr, sep='\t')
	corr_df.to_csv(edge, sep='\t', index=False)
