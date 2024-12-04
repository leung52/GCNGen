
import pandas as pd


def merge_single_cell_gene_expression_data_files(file1: str, file2: str, output_file: str) -> pd.DataFrame:
	print('Loading dataframes...')
	df1 = pd.read_csv(file1, sep='\t', index_col=0)
	df2 = pd.read_csv(file2, sep='\t', index_col=0)
	print('done.')
	print('Merging dataframes into one...')
	merge_df = df1.merge(df2, how='outer', on='GENE') # merge dataframes on genes
	print('done.')
	print('Replacing NaN with int(0)...')
	merge_df = merge_df.fillna(0) # convert NaN to 0
	print('done.')
	print('Ensuring data type is int...')
	merge_df = merge_df.astype(int) # ensure data type is int
	print('done.')
	print('Saving dataframe to tsv...')
	merge_df.to_csv(output_file, sep='\t') # Save to tsv file
	print('done.')
	return merge_df


if __name__ == "__main__":
	tmu_file_1 = 'Tmu_scRNASeq_objects/DS50_clean_DGE_3000.txt' # tab separated
	tmu_file_2 = 'Tmu_scRNASeq_objects/DS51_clean_DGE_3000.txt'
	output_file = 'processed_data/DS50_DS51_clean_DGE_combined_6000.tsv'
	merge_single_cell_gene_expression_data_files(tmu_file_1, tmu_file_2, output_file)
