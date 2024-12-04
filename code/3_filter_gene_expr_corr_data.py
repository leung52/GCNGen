import multiprocessing as mp
import pandas as pd
import numpy as np


def _init_p(p_value: float) -> None:
	if 0 < p_value <= 1:
		global _p
		_p = p_value
	else:
		raise ValueError("P-value must be smaller than or equal to 1 or greater than to 0")

def _init_min_corr(min_corr_value: float) -> None:
	if 0 <= min_corr_value <= 1:
		global _min_corr
		_min_corr = min_corr_value
	else:
		raise ValueError("Minimun correlation coefficient value must be between 0 and 1 inclusive.")

def filter_p_value(df: pd.DataFrame) -> pd.DataFrame:
	return df[df['p-value'] < _p]

def filter_correlation_coefficient(df: pd.DataFrame) -> pd.DataFrame:
	return df[df['Correlation'] >= _min_corr]

def filter_p_and_corr(df: pd.DataFrame) -> pd.DataFrame:
	return filter_p_value(filter_correlation_coefficient(df))

def filter_with_mp_and_save_to_tsv(input_file: str, output_file: str, filter_func) -> None:
	# Have to initialise global variables
	with mp.Pool(processes=mp.cpu_count()) as pool:
		iter_df = pd.read_csv(input_file, sep='\t', chunksize=5000)
		header_written = False
		for filtered_df in pool.imap_unordered(filter_func, iter_df):
			if not filtered_df.empty:
				filtered_df.to_csv(output_file, sep='\t', header=not header_written, index=False, mode='a')
				header_written=True


if __name__ == "__main__":
	file = "processed_data/filtered_tmu_genes_kendall_b_genes_exp_5_cells_p_005.tsv"
	out = "processed_data/filtered_tmu_genes_kendall_b_genes_exp_5_cells_p_005_corr_99_percentile.tsv"
	corr_coef = np.array(pd.read_csv(file, sep='\t', usecols = ['Correlation']))
	cutoff = np.quantile(corr_coef, 0.99)
	del corr_coef
	print(cutoff)

	_init_min_corr(cutoff)
	filter = filter_correlation_coefficient

	filter_with_mp_and_save_to_tsv(file, out, filter)
