#!/usr/bin/bash
#This script allows for all fastq files to be read in a specfic path.
#I have used this for 6 fastq files to 500 fastq files.  
#paths may need to be modified for datasets
#code linking fastq files together may need to be modified based on file names

import os, subprocess

#paths to data, RNA program, and where python script is run from.

data_path = '/media/New_Volume/Justin/C202SC18031469/raw_data'
salmon_path = '/media/New_Volume/Justin/RNA_Programs/salmon/bin'
origin_path_of_program = '/media/New_Volume/Justin/RNA_Programs'
read_path = '/media/New_Volume/Justin/C202SC18031469/raw_data'
list_of_reads = []
pairs_of_reads = []

#obtain all data file names from path and adds to a list
for read_file in os.listdir(data_path):
	list_of_reads.append(read_file)
list_of_reads.sort()

#Creates a printout showing which files will be run with which files
#So user can make sure they are lining up correctly
for i in range(0, len(list_of_reads), 2):
	#start at 0 and move by 2's so the last iteration is the len -3
	if i == len(list_of_reads) - 3:
		#print(i)
		pairs_of_reads.append(list_of_reads[i]+ ': ' + list_of_reads[i+1])
		break
	else:
		#print(i)
		pairs_of_reads.append(list_of_reads[i] + ': ' +list_of_reads[i+1])
		i += 1
for i in pairs_of_reads:
	print(i)

for i in range(0, len(list_of_reads), 2):
	#start at 0 and move by 2's so the last iteration is the len -3
	#print(list_of_reads[i][0:22],list_of_reads[i+1][0:22])
	#print(list_of_reads[i].split('.')[1][-5:], list_of_reads[i+1].split('.')[1][-5:])
	if i == len(list_of_reads) - 3:
		#print(i)

		#in this case the file names were 22 characters long.  Checks if they match
		#then checks if they are read1 and read2
		if list_of_reads[i][0:22] == list_of_reads[i+1][0:22]:
			if list_of_reads[i].split('.')[1][-2:] == 'R1' and list_of_reads[i+1].split('.')[1][-2:]== 'R2':
				print(i,True)
			else:print(i, False)
		else:
			print(i, False)
		break
	else:
		if list_of_reads[i][0:22] == list_of_reads[i+1][0:22]:
			if list_of_reads[i].split('.')[1][-2:] == 'R1' and list_of_reads[i+1].split('.')[1][-2:]== 'R2':
				print(i,True)
			else:print(i, False)
		else:
			print(i, False)

# for i in pairs_of_reads:
# 	print(i)


#creates file name output for salmon function.
#Salmon must be run as a bash exectible.  Not a sh which is the default for os and running_subprocess
#which is why in the subprocess we call exectuable bash

#This portion tells the computer to run all the fastq files through salmon so we don't have to do them one at a time
#because for some reason people waste their day doing that

num_remaining = len(list_of_reads)/2
for i in range(0, len(list_of_reads), 2):

	partial_path_read_one = read_path + '/%s' %list_of_reads[i]
	partial_path_read_two = read_path + '/%s' %list_of_reads[i+1]
	full_path_read_one = '"' + partial_path_read_one + '"'
	full_path_read_two = '"' + partial_path_read_two + '"'
	# print(full_path_read_one)
	# print(full_path_read_two)
	temp_output_name = list_of_reads[i].replace('.','_').split('_')[0:4]
	output_name = '_'.join(temp_output_name)
	print('Currently running Sample: ' + output_name + '\n')
	command = './salmon quant -i transcripts_index_mouse_salmon_keepDuplicates -l A -1 <(gunzip -c %s) -2 <(gunzip -c %s) -o /media/New_Volume/Justin/salmon_mouse_5_4_2018/%s' %(full_path_read_one, full_path_read_two, output_name)
	os.chdir(salmon_path)
	print('running_subprocess')
	subprocess.call([command], shell=True, executable='/bin/bash')
	num_remaining = num_remaining - 1
	print('\nNumber of data files remaining: ' + str(num_remaining) + '\n')

os.chdir(origin_path_of_program)




