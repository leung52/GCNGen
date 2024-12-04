import pandas as pd
from itertools import combinations
import scipy.stats as ss
import multiprocessing as mp
from tqdm import tqdm
import os


_CORR_FUNCTIONS = {
	'kendall': ss.kendalltau
}
_gene_exp_df = None
_corr_method = None

def load_single_cell_expression_data(file_path: str) -> pd.DataFrame:
	print('Loading gene expression data...')
	df = pd.read_csv(file_path, sep='\t', index_col=0)
	df = df.T  # Transpose: columns become rows and vice versa
	print('done.')
	return df

def list_gene_pairs(gene_exp_df: pd.DataFrame) -> list:
	expected_pairs = len(gene_exp_df.columns) * (len(gene_exp_df.columns) - 1) // 2
	gene_pairs_ls = list(tqdm(combinations(gene_exp_df.columns, 2), total=expected_pairs, desc='Gathering gene pairs'))
	return gene_pairs_ls

def chunkit(ls: list, chunksize: int) -> list:
	print(f'Chunking list...')
	chunked_ls = [ls[i:i + chunksize] for i in range(0, len(ls), chunksize)]
	print('done.')
	return chunked_ls

def generate_interaction_file(gene_exp_df, gene_pairs_chunked_ls, corr_method, output_file) -> None:
	if _init_worker(gene_exp_df, corr_method):
		header_written = False
		with mp.Pool(processes=mp.cpu_count()) as pool:
			for corr_df in tqdm(pool.imap_unordered(_generate_corr_df, gene_pairs_chunked_ls), total=len(gene_pairs_chunked_ls)):
				corr_df.to_csv(output_file, mode='a', sep='\t', index=False, header=not header_written)
				header_written = True

def _init_worker(df: pd.DataFrame, corr_method: str) -> bool:
	"""Initialize the global variable for each worker."""
	if _check_corr_method_ok(corr_method):
		global _corr_method
		_corr_method = _CORR_FUNCTIONS[corr_method]
		global _gene_exp_df
		_gene_exp_df = df
		return True
	else:
		return False

def _check_corr_method_ok(corr_method: str) -> bool:
	if corr_method in _CORR_FUNCTIONS:
		return True
	else:
		print(f"'{corr_method}' is not a valid argument. Please use one of the following:")
		for key in _CORR_FUNCTIONS.keys():
			print(f"        - '{key}'")
		return False

def _generate_corr_df(gene_pairs_ls: list) -> pd.DataFrame:
	corr_dict = {
		'GENE1': [],
		'GENE2': [],
		'Correlation': [],
		'p-value': []
	}

	for gene1, gene2 in gene_pairs_ls:
		t, p = _corr_method(_gene_exp_df[gene1], _gene_exp_df[gene2])
		corr_dict['GENE1'].append(gene1)
		corr_dict['GENE2'].append(gene2)
		corr_dict['Correlation'].append(t)
		corr_dict['p-value'].append(p)

	return pd.DataFrame(corr_dict)

def combine_files(file_names: list, output_file: str) -> None:
	with open(output_file, 'w') as out_f:
		first = True
		for file in file_names:
			with open(file, 'r') as f:
				if not first:
					next(f)
				for line in f:
					out_f.write(line)
			os.remove(file)
			first = False

def calculate_correlation(gene_exp_file: str, out_file: str, corr_method: str) -> None:
	gene_exp_df = load_single_cell_expression_data(gene_exp_file)
	gene_pairs_ls = list_gene_pairs(gene_exp_df)
	gene_pairs_ls = chunkit(gene_pairs_ls, 250000)
	x = len(gene_pairs_ls)
	file_base = 'out/temporary'
	files = []

	for n, ls in enumerate(gene_pairs_ls):
		n += 1
		ls = chunkit(ls, 1000)
		out_f = file_base + str(n) + '.tsv'
		print(f'File {n}/{x} being generated')
		generate_interaction_file(gene_exp_df, ls, corr_method, out_f)
		files.append(out_f)

	combine_files(files, out_file)


if __name__ == "__main__":
	file1 = 'processed_data/DS50_DS51_clean_DGE_combined_6000.tsv'
	file2 = 'processed_data/filtered_for_genes_expressed_in_more_than_5_percent_of_cells.tsv'
	file3 = 'processed_data/filtered_for_genes_expressed_in_more_than_10_percent_of_cells.tsv'

	out = 'processed_data/tmu_genes_kendall_b_genes_in_more_than_10_percent_of_cells.tsv'

	calculate_correlation(file3, out, 'kendall')
