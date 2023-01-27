
#!/usr/bin/env python3

import os
import csv
from numpy import NAN
import pandas as pd
import glob
from pathlib import Path
import run_out_constants as c

'''

This script will scrape the dragen metrics files and output a csv file for each seq type

usage: python3 dragen_XXXXmetrics_scraper.py

'''

# search through the xxxx_metrics directories for the metrics files
def file_finder(pattern, path):
    
    search_results = []
    
    for root, dirs, files in os.walk(path):
        for name in files:
            if pattern in name:
                search_results.append(os.path.join(root,name))
    
    return search_results

seq_type = [
            'genome_seq', 
            'exome_seq'
            ]

# gather the files from the xxxx_metrics directories
exome_coverage = file_finder(c.QC_COV_EXT, c.EXO_XXXX_COV) + file_finder(c.QC_COV_EXT2, c.EXO_XXXX_COV)
exome_mapping = file_finder(c.MAPPING_EXT, c.EXO_XXXX_MAP) 
genome_coverage = file_finder(c.COV_EXT, c.WGS_XXXX_COV)
genome_mapping = file_finder(c.MAPPING_EXT, c.WGS_XXXX_MAP)

#make a dict of the files by seq type
seq_type_dict = {
                'genome_seq': {'coverage' : genome_coverage, 'mapping' : genome_mapping}, 
                'exome_seq': {'coverage' : exome_coverage, 'mapping' : exome_mapping}
                }


#iterate through the seq types
for seq in seq_type:

    with open (f'{seq}_mapping.csv', 'w') as f1, open (f'{seq}_coverage.csv', 'w') as f2:
        
        #make the fieldnames    
        mapping_fieldnames_dict = c.XXXX_DICT_FIELDNAME_REF[seq]['MAPPING/ALIGNING SUMMARY']
        coverage_fieldnames_dict = c.XXXX_DICT_FIELDNAME_REF[seq]['COVERAGE SUMMARY']

        #make 2 dict writer objects each time
        #write fieldnames 
        mapping_writer = csv.DictWriter(f1, fieldnames = dict.fromkeys(c.DICT_FIELDNAME_REF[seq]['MAPPING/ALIGNING SUMMARY'].values()), delimiter=',')
        coverage_writer = csv.DictWriter(f2, fieldnames = dict.fromkeys(c.DICT_FIELDNAME_REF[seq]['COVERAGE SUMMARY'].values()), delimiter=',')

        #write header row
        mapping_writer.writeheader() 
        coverage_writer.writeheader() 

        #iterate through the reports
        for type, report in seq_type_dict.items():

            for key, files in report.items():

                if key == 'coverage':

                    files = seq_type_dict[seq]['coverage']

                    #loop through the files in the seq_type_dict
                    for file in files:
                        
                        #open the file
                        with open(file, 'r') as f:
                            
                            #get the sample name
                            sample = os.path.basename(file).split('.')[-3]
                            #get the work order
                            # work_order = os.path.basename(file).split('/')[-2]

                            file_data = {}
                            
                            file_data.update({'COVERAGE SUMMARY': {sample: dict.fromkeys(list(coverage_fieldnames_dict.values()), NAN)}})
                            file_data['COVERAGE SUMMARY'][sample]['Sample'] = sample
                            # file_data['COVERAGE SUMMARY'][sample]['Work Order'] = work_order

                            #iterate through the lines in the file
                            for l in f:

                                #split the line into a list
                                line = l.split(',')

                                #strip the whitespace from the 2nd index position of the list
                                line[2] = line[2].replace(' ', '')
                                line[-1] = line[-1].strip()
                                
                                if len(line) == 4:
                                    
                                    value = line[2]
                                    
                                    value_fieldname = coverage_fieldnames_dict.get(value)
                                    file_data[line[0]][sample][value_fieldname] = line[3]

                                    continue

                                if len(line) == 5:

                                    value = line[2].strip()            
                                    value_fieldname = coverage_fieldnames_dict.get(value)
                                    
                                    if value_fieldname == '% Aligned bases genome' or value_fieldname == '% Aligned reads genome':

                                        file_data['COVERAGE SUMMARY'][sample][value_fieldname] = float(line[4])
                                        continue
                                    
                                    continue

                        # print('coverage_data: ', file_data)

                if key == 'mapping':

                    files = seq_type_dict[seq]['mapping']

                    #loop through the files in the seq_type_dict
                    for file in files:
                        
                        #open the file
                        with open(file, 'r') as f:
                            
                            f_reader = csv.reader(f, delimiter=',')
                            #get the sample name
                            sample = os.path.basename(file).split('.')[-3]
                            #get the work order
                            # work_order = os.path.basename(file).split('/')[-2]

                            file_data = {}
                            
                            file_data.update({'MAPPING/ALIGNING SUMMARY': {sample: dict.fromkeys(list(mapping_fieldnames_dict.values()), NAN)}})
                            file_data['MAPPING/ALIGNING SUMMARY'][sample]['Sample'] = sample
                            # file_data['MAPPING/ALIGNING SUMMARY'][sample]['Work Order'] = work_order

                            #iterate through the lines in the file
                            for line in f_reader:
                                
                                if  line[0] == 'MAPPING/ALIGNING SUMMARY':
                                    
                                    #split the line into a list
                                    # line = l.split(',')

                                    #strip the whitespace from the 2nd index position of the list
                                    line[2] = line[2].replace(' ', '')
                                    line[-1] = line[-1].strip()
                                
                                    # print('line: ', line)
                                   

                                    if len(line) == 4:
                                    
                                        
                                        value = line[2]
                                        
                                        value_fieldname = mapping_fieldnames_dict.get(value)
                                        file_data[line[0]][sample][value_fieldname] = line[3]

                                        continue

                                    if len(line) == 5:


                                        value = line[2]   
                  
                                        value_fieldname = mapping_fieldnames_dict.get(value)
                                        percent_fieldname  = mapping_fieldnames_dict.get(value_fieldname)
                                        
                                        if not percent_fieldname:
                                            print('filename', f)
                                            print('value_fieldname', value_fieldname)
                                            print('percent_lookup', percent_fieldname)
                                            exit()

                                        file_data[line[0]][sample][value_fieldname] = line[3]
                                        file_data[line[0]][sample][percent_fieldname] = line[4]

                                        continue

                        # print('mapping_data: ', file_data)
                        
                        # print('dictionary as gollows')
                        # for k,v in file_data['MAPPING/ALIGNING SUMMARY'][sample].items():
                        #     print(k,v)
                            

    # write rows
    # print('Writing rows to csv files...')
    
    # for sample in combined_data['MAPPING/ALIGNING SUMMARY']:
    #     mapping_writer.writerow(combined_data['MAPPING/ALIGNING SUMMARY'][sample])

    # for sample in combined_data['COVERAGE SUMMARY']:
    #     coverage_writer.writerow(combined_data['COVERAGE SUMMARY'][sample])
