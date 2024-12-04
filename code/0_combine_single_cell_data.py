from functools import reduce
import pandas as pd


def merge_sc_data_files(files: list(str), save_to: str) -> None:
	print('Loading dataframes...')
	dfs = []
	for file in files:
		dfs.append(pd.read_csv(file, sep='\t', index_col=0))
	print('done.')

	print('Merging dataframes into one...')
	df_merged = reduce(lambda left,right: pd.merge(left,right,on=['GENE'], how='outer'), dfs)
	del dfs
	print('done.')

	print('Replacing NaN with int(0)...')
	merge_df = merge_df.fillna(0)
	print('done.')

	print('Ensuring data type is int...')
	merge_df = merge_df.astype(int)
	print('done.')

	print('Saving dataframe to tsv...')
	merge_df.to_csv(save_to, sep='\t')
	print('done.')


if __name__ == "__main__":
	tmu_file_1 = 'Tmu_scRNASeq_objects/DS50_clean_DGE_3000.txt' # tab separated
	tmu_file_2 = 'Tmu_scRNASeq_objects/DS51_clean_DGE_3000.txt'
	output_file = 'processed_data/DS50_DS51_clean_DGE_combined_6000.tsv'
	merge_sc_data_files(tmu_file_1, tmu_file_2, output_file)

