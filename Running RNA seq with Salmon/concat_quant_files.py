import os
import pandas as pd

#############################
###########INPUTS############
#############################

path_to_results = '/media/New_Volume/Justin/salmon_results_3_6_ALLTRANSCRIPTS'
path_to_output_df = '/media/New_Volume/Justin'
name_of_output_df = 'prostateLiang_blood_df.txt'


#############################
############Code#############
#############################

#load dataframe with transcripts other data (gene name, chrom, start, end)
data_df = pd.read_csv('transcriptome_data.txt', delimiter='\t')
data_df.columns = ['Name', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position']
print(data_df[:20])
list_of_reads_results = []

#creates list of filenames from the data results from salmon
for result_file in os.listdir(path_to_results):
	#skips annoying '.' files that are added in background, for example '.DS_Store'
	if result_file[0:1] == '.':
		pass
	else:
		list_of_reads_results.append(result_file)
list_of_reads_results.sort()
# for i in list_of_reads_results:
# 	print(i)
#print(list_of_reads_results)

#creates df's from quant.sf files in separtate folders in the path files.
#Loop will change directories to proper directory, merge dfs, and label data with
#last 10 characters of folder name in path.  Finally, outputs df to tab separated txt file
#in designated path of output location.
count = 0
for data_set in list_of_reads_results:
	#creates first df to initialize loop
	if count == 0:
		new_path = path_to_results + '/' + data_set
		os.chdir(new_path)
		ID_of_sample = data_set
		temp_df=pd.read_csv('quant.sf', delimiter='\t')
		results_df = temp_df[['Name', 'NumReads']].copy()
		results_df.columns = ['Name', ID_of_sample]
		results_df = pd.merge(data_df, results_df, on='Name')
		print(results_df[:20])
		count += 1
		os.chdir(path_to_results)

	# elif count >5:
	# 	break

	else:
		new_path = path_to_results + '/' + data_set
		os.chdir(new_path)
		ID_of_sample = data_set
		temp_df=pd.read_csv('quant.sf', delimiter='\t')
		numReads_df = temp_df[['Name', 'NumReads']].copy()
		numReads_df.columns = ['Name', ID_of_sample]
		#print(numReads_df[:20])
		results_df = pd.merge(results_df, numReads_df, on='Name')
		print('Sucessful Merge of: ' + data_set)

		count += 1
		os.chdir(path_to_results)
#print(results_df[:20])
os.chdir(path_to_output_df)
print('Outputting dataframe to .txt file...')
results_df.to_csv('%s' %name_of_output_df, header = True, index=False, sep='\t')
print('Sucessfully created .txt file. \n Program Terminate')
