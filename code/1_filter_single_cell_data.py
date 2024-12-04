import pandas as pd
import multiprocessing as mp



def _init_gene_cutoff(cutoff: float) -> None:
	if 0 < cutoff < 1:
		global _gene_cutoff
		_gene_cutoff = cutoff
	else:
		raise ValueError("Cutoff must be between 0 and 1.")


def _init_cell_cutoff(cutoff: int) -> None:
	if cutoff > 0:
		global _cell_cutoff
		_cell_cutoff = cutoff
	else:
		raise ValueError("Cutoff must be greater than 0.")

def remove_genes(df: pd.DataFrame) -> pd.DataFrame:
	"""Remove gene rows where gene is expressed in less than _gene_cutoff of cells."""
	n = round(_gene_cutoff * len(df.columns))
	return df[(df != 0).sum(axis=1) >= n]

def remove_cells(df: pd.DataFrame) -> pd.DataFrame:
	return df.loc[:, (df != 0).sum(axis=0) >= _cell_cutoff]


if __name__ == "__main__":
	_init_gene_cutoff(0.05)
	with mp.Pool(processes=mp.cpu_count()) as pool:
		iter_df = pd.read_csv('processed_data/tmu_genes_kendall_b_genes_in_more_than_5_percent_of_cells.tsv', sep='\t', chunksize=1000)
		header_written = False
		output_file = 'processed_data/filtered_for_genes_expressed_in_more_than_5_percent_of_cells.tsv'
		for filtered_df in pool.imap_unordered(remove_genes, iter_df):
			if not filtered_df.empty:
				filtered_df.to_csv(output_file, sep='\t', header=not header_written, index=False, mode='a')
				header_written=True
