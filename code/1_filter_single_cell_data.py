import pandas as pd
import multiprocessing as mp


def filter_genes_by_cell_expr(percent_cutoff: float, data_file: str, save_to: str) -> None:
	_set_gene_cutoff(percent_cutoff)
	_use_mp_to_filter(_remove_genes, data_file, save_to)

def filter_cells_by_genes_expr(n_genes: int, data_file: str, save_to: str) -> None:
	_set_cell_cutoff(n_genes)
	_use_mp_to_filter(_remove_cells, data_file, save_to)

def _use_mp_to_filter(filter, data_file: str, save_to: str) -> None:
	with mp.Pool(processes=mp.cpu_count()) as pool:
		iterable_dfs = pd.read_csv(data_file, sep='\t', chunksize=1000)
		header_written = False
		for filtered_df in pool.imap_unordered(filter, iterable_df):
			if not filtered_df.empty:
				filtered_df.to_csv(save_to, sep='\t', header=not header_written, index=False, mode='a')
				header_written=True

def _remove_genes(df: pd.DataFrame) -> pd.DataFrame:
	"""Remove gene rows where gene is expressed in less than _gene_cutoff of cells."""
	n = round(_gene_cutoff * len(df.columns))
	return df[(df != 0).sum(axis=1) >= n]

def _remove_cells(df: pd.DataFrame) -> pd.DataFrame:
	return df.loc[:, (df != 0).sum(axis=0) >= _cell_cutoff]

def _set_gene_cutoff(cutoff: float) -> None:
	if 0 < cutoff < 1:
		global _gene_cutoff
		_gene_cutoff = cutoff
	else:
		raise ValueError("Cutoff must be between 0 and 1.")

def _set_cell_cutoff(cutoff: int) -> None:
	if cutoff > 0:
		global _cell_cutoff
		_cell_cutoff = cutoff
	else:
		raise ValueError("Cutoff must be greater than 0.")


if __name__ == "__main__":
	pass
