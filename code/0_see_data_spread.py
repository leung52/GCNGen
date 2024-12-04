import pandas as pd
import matplotlib.pyplot as plt


def load_gene_expr_data(file: str) -> pd.DataFrame:
	return pd.read_csv(file, sep='\t', index_col='GENE')

def all_data(df: pd.DataFrame) -> list:
	return list(df.values.flatten())

def count_genes_expressed_in_cell(df: pd.DataFrame) -> list:
	return list(df.sum(axis=0))

def count_cells_expressing_gene(df: pd.DataFrame) -> list:
	return list(df.sum(axis=1))

def load_correlation_df(file: str) -> pd.DataFrame:
	return pd.read_csv(file, sep='\t')

def get_correlation_coefficients(df: pd.DataFrame) -> list:
	return list(df['Correlation'])

def plot(data: list, file_name: str, log_scale: bool) -> None:
	plt.hist(data, bins=30, log=log_scale, color='xkcd:sky blue')
	plt.savefig(file_name)


if __name__ == "__main__":
	files = []

	for x in ["005", "001", "0001"]:
		file = f'processed_data/filtered_tmu_genes_kendall_b_genes_exp_5_cells_p_{x}.tsv'
		files.append(file)

	for file in files:
		ls = get_correlation_coefficients(load_correlation_df(file))
		line = (file.split('/')[-1]).split('.')[0]
		f = f'figures/log_{line}.png'
		plot(ls, f, True)
