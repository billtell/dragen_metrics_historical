#!/usr/bin/env python3

import os
import csv
from numpy import NAN
import pandas as pd
import glob
from pathlib import Path
import run_out_constants as c


'''
This script will parse the dragen run.out files and create a csv file with the metrics.

usage: python3 dragen_run_out.py 

'''

TEST_DIR = '/home/gteller/run_test'

#search for the run.out files
def get_out_files(dirpath):

    return list(Path(dirpath).rglob('*.out'))


seq_type = [
            'genome_seq', 
            'exome_seq'
            ]


#iterate through the seq types
for seq in seq_type:

    with open (f'{seq}_mapping.csv', 'w') as f1, open (f'{seq}_coverage.csv', 'w') as f2:

        #make the fieldnames    
        mapping_fieldnames_dict = c.DICT_FIELDNAME_REF[seq]['MAPPING/ALIGNING SUMMARY']
        coverage_fieldnames_dict = c.DICT_FIELDNAME_REF[seq]['COVERAGE SUMMARY']
        
        #make 2 dict writer objects each time
        #write fieldnames 
        mapping_writer = csv.DictWriter(f1, fieldnames = dict.fromkeys(c.DICT_FIELDNAME_REF[seq]['MAPPING/ALIGNING SUMMARY'].values()), delimiter=',')
        coverage_writer = csv.DictWriter(f2, fieldnames = dict.fromkeys(c.DICT_FIELDNAME_REF[seq]['COVERAGE SUMMARY'].values()), delimiter=',')
        
        #write header row
        mapping_writer.writeheader() 
        coverage_writer.writeheader() 

        #get the files
        files = get_out_files(os.path.join(TEST_DIR, seq))
        
        #write the data
        for file in files:
            
            #open the file
            with open (file, 'r') as f:
                
                #sample name
                #user thing. maybe
                sample = os.path.basename(file).split('.')[-3]
                
                if seq == 'genome_seq':
                    work_order = os.path.basename(file).split('.')[1]
                    if work_order == 'genome_seq' or work_order == 'DRAGEN':
                        work_order = os.path.basename(file).split('.')[2]
                    if work_order == 'DRAGEN' or work_order == 'somatic' or work_order == 'GERMLINE' or work_order == 'complete':
                        work_order = os.path.basename(file).split('.')[0]

                if seq == 'exome_seq':
                    work_order == os.path.basename(file).split('.')[1]
                    if work_order == 'exome_seq' or work_order == 'original' or work_order == 'all_reads' or work_order == 'germline' or work_order == 'ALL_READS':
                        work_order = os.path.basename(file).split('.')[0]
                
                data = {}
                data.update({'MAPPING/ALIGNING SUMMARY': {sample: dict.fromkeys(list(mapping_fieldnames_dict.values()), NAN)}, 
                             'COVERAGE SUMMARY': {sample: dict.fromkeys(list(coverage_fieldnames_dict.values()), NAN)}})

                data['MAPPING/ALIGNING SUMMARY'][sample]['Sample'] = sample
                data['MAPPING/ALIGNING SUMMARY'][sample]['Work Order'] = work_order
                data['COVERAGE SUMMARY'][sample]['Sample'] = sample
                data['COVERAGE SUMMARY'][sample]['Work Order'] = work_order
                
                dragen_metrics_found = False

                #capture the version line
                for l in f:

                    if 'DRAGEN Host Software Version' in l:
                        version = l.split(' ')[4]
                        data['MAPPING/ALIGNING SUMMARY'][sample]['Version'] = version
                        data['COVERAGE SUMMARY'][sample]['Version'] = version
                        continue
                        
                    #only read in the lines after line containing 'DRAGEN METRICS'
                    if l.startswith('DRAGEN METRICS'):
                        dragen_metrics_found = True
                        continue

                    if dragen_metrics_found:

                        if 'MAPPING/ALIGNING SUMMARY' in l:
                            space_count = 0
                            string_cat = False
                            strings = []
                            string_word = ''

                            for i, string in enumerate(l, 1):

                                if string == ' ':
                                    space_count += 1

                                if string != ' ' and space_count == 0:
                                    string_cat = True
                                    string_word += string
                                    space_count = 0

                                if string != ' ' and space_count == 1:
                                    string_cat = True
                                    string_word += f' {string}'
                                    space_count = 0

                                if space_count == 2:
                                    if string_word:
                                        strings.append(string_word)

                                    string_cat = False
                                    string_word = ''
                                    space_count = 0

                                if i == len(l): 
                                    if string_word:
                                        if 'NA\n' in string_word or ' \n' in string_word:
                                            string_word = NAN
                                        if string_word == '\n':
                                            string_word = NAN

                                        strings.append(string_word)

                            #strip the new line character from strings
                            #strip the newline character
                            for i, val in enumerate(strings):
                                if isinstance(val, str) and val.endswith('\n'):
                                    strings[i] = float(val.strip())
                                    #strip the whitespace
                                if isinstance(val, str) and val.startswith(' '):
                                    strings[i] = str(val.strip())

                            #add to dictionary data{'MAPPING/ALIGNING SUMMARY': {metric: {'value': value, 'percent': percent}
                            # print(strings)

                            if len(strings) == 3 or len(strings) == 4:
                            
                                value = strings[1]
                                
                                value_fieldname = mapping_fieldnames_dict.get(value)
                                data['MAPPING/ALIGNING SUMMARY'][sample][value_fieldname] = float(strings[2])
                            
                            if len(strings) == 5:
                    
                                value = strings[1]

                                value_fieldname = mapping_fieldnames_dict.get(value)
                                percent_fieldname = mapping_fieldnames_dict.get(value_fieldname)
                                
                                #example of what it is grabbing, these correspond to the header row
                                
                                #"Number of duplicate marked reads" : "Duplicate reads",
                                    # "Duplicate reads" : "%Duplicate reads"
                                
                                data['MAPPING/ALIGNING SUMMARY'][sample][value_fieldname] = float(strings[2])
                                data['MAPPING/ALIGNING SUMMARY'][sample][percent_fieldname] = float(strings[3])  
                                 
                        if 'COVERAGE SUMMARY' in l:
                            space_count = 0
                            string_cat = False
                            strings = []
                            string_word = ''

                            if 'PCT of genome with coverage' in l  or 'PCT of target region with coverage' in l:

                                #split the line into a list
                                l = l.split(' ')
                                #remove the newline character
                                l = [i.strip() for i in l]
                                #remove the empty strings
                                l = [i for i in l if i]
                                #concatenate the strings except for index position 0,1,-1
                                l = [l[0] + ' ' + l[1],' '.join(l[2:-1]), l[-1]]
                                
                                cov_header = l[1]
                                #add to dictionary 
                                data['COVERAGE SUMMARY'][sample][coverage_fieldnames_dict.get(cov_header)] = l[2]

                                continue

                            for i, string in enumerate(l, 1):

                                if string == ' ':
                                    space_count += 1

                                if string != ' ' and space_count == 0:
                                    string_cat = True
                                    string_word += string
                                    space_count = 0

                                if string != ' ' and space_count == 1:
                                    string_cat = True
                                    string_word += f' {string}'
                                    space_count = 0

                                if space_count == 2:
                                    if string_word:
                                        strings.append(string_word)
                                    string_cat = False
                                    string_word = ''
                                    space_count = 0

                                if i == len(l): 
                                    if string_word:
                                        if 'NA\n' in string_word or ' \n' in string_word:
                                            string_word = NAN
                                        if string_word == '\n':
                                            string_word = NAN

                                        strings.append(string_word)

                            #add to dictionary data{'COVERAGE SUMMARY': {metric: {'value': value, 'percent': percent}}}
                            if len(strings) == 3:
                                
                                value = strings[1].strip()
                                data['COVERAGE SUMMARY'][sample][coverage_fieldnames_dict.get(value)] = float(strings[2])
                                #if the result of the .get is None, then print the value    
                                continue

                            if len(strings) >= 4:
                                
                                value = strings[1].strip()
                                
                                value_fieldname = coverage_fieldnames_dict.get(value)
                                
                                if value_fieldname == '% Aligned bases genome' or value_fieldname == '% Aligned reads genome':

                                    data['COVERAGE SUMMARY'][sample][value_fieldname] = float(strings[3])
                                    continue
                                
                                data['COVERAGE SUMMARY'][sample][value_fieldname] = float(strings[2])
                                
                                continue
            
        # print('run_out data:', data)

        #write the data to the csv files
        # for sample in data['MAPPING/ALIGNING SUMMARY']:
        #     mapping_writer.writerow(data['MAPPING/ALIGNING SUMMARY'][sample])

        # for sample in data['COVERAGE SUMMARY']:
        #     coverage_writer.writerow(data['COVERAGE SUMMARY'][sample])
    
        #pull metrics from files in the seq_type_dict
        # #loop through the seq_type_dict


