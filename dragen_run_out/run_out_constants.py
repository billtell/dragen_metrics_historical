
#pathways
EXO_XXXX_MAP = '/storage1/fs1/gtac-mgi/Active/Bioinformatics_analysis/exome_seq/xxxx_metrics/mapping'
EXO_XXXX_COV = '/storage1/fs1/gtac-mgi/Active/Bioinformatics_analysis/exome_seq/xxxx_metrics/coverage'

WGS_XXXX_MAP = '/storage1/fs1/gtac-mgi/Active/Bioinformatics_analysis/genome_seq/xxxx_metrics/mapping'
WGS_XXXX_COV = '/storage1/fs1/gtac-mgi/Active/Bioinformatics_analysis/genome_seq/xxxx_metrics/wgs_coverage'


#extensions
COV_EXT = '.wgs_coverage_metrics.csv'
QC_COV_EXT = '.qc-coverage-region-1_coverage_metrics.csv' 
QC_COV_EXT2 =  '.qc-coverage-region-2_coverage_metrics.csv'
MAPPING_EXT = '.mapping_metrics.csv'

#header dicts and column ordering lists
COV_HEADING_v2 = {
    "Sample" : "Sample",
    "version" : "Version",
    "Work Order" : "Work Order",
    "Aligned bases" : "Aligned bases",
    "Aligned bases in genome" :  "% Aligned bases genome",
    "Average alignment coverage over genome" : "Average coverage genome",
    "Uniformity of coverage (PCT > 0.2*mean) over genome" : "Uniformity of cov genome (PCT > 0.2*mean)",
    "Uniformity of coverage (PCT > 0.4*mean) over genome" : "Uniformity of cov genome (PCT > 0.4*mean)",
    "PCT of genome with coverage [1500x: inf)" : "PCT of genome [1500x: inf)",                              
    "PCT of genome with coverage [1000x: inf)" : "PCT of genome [1000x: inf)",                            
    "PCT of genome with coverage [ 500x: inf)" : "PCT of genome [500x: inf)",                                 
    "PCT of genome with coverage [ 100x: inf)" : "PCT of genome [100x: inf)",                           
    "PCT of genome with coverage [ 50x: inf)" : "PCT of genome [50x: inf)",                              
    "PCT of genome with coverage [ 20x: inf)" : "PCT of genome [20x: inf)",                           
    "PCT of genome with coverage [ 15x: inf)" : "PCT of genome [15x: inf)",                                
    "PCT of genome with coverage [ 10x: inf)" : "PCT of genome [10x: inf)",                                   
    "PCT of genome with coverage [ 3x: inf)" : "PCT of genome [3x: inf)",                                 
    "PCT of genome with coverage [ 1x: inf)" : "PCT of genome [1x: inf)",                                 
    "PCT of genome with coverage [ 0x: inf)" : "PCT of genome [0x: inf)",                                 
    "PCT of genome with coverage [1000x:1500x)" : "PCT of genome [1000x:1500x)",                                
    "PCT of genome with coverage [ 500x:1000x)" : "PCT of genome [500x:1000x)",                               
    "PCT of genome with coverage [ 100x: 500x)" : "PCT of genome [100x:500x)",                                  
    "PCT of genome with coverage [ 50x: 100x)" : "PCT of genome [50x:100x)",                                 
    "PCT of genome with coverage [ 20x: 50x)" : "PCT of genome [20x:50x)",                             
    "PCT of genome with coverage [ 15x: 20x)" : "PCT of genome [15x:20x)",                                 
    "PCT of genome with coverage [ 10x: 15x)" : "PCT of genome [10x:15x)",                                    
    "PCT of genome with coverage [ 3x: 10x)" : "PCT of genome [3x:10x)",                                
    "PCT of genome with coverage [ 1x: 3x)" : "PCT of genome [1x:3x)",                                   
    "PCT of genome with coverage [ 0x: 1x)" : "PCT of genome [0x:1x)",                              
    "Average chr X coverage over genome" : "Ave. chrX coverage genome",
    "Average chr Y coverage over genome" : "Ave. chrY coverage genome",
    "Average mitochondrial coverage over genome" : "Ave. chrM coverage genome",
    "Average autosomal coverage over genome" : "Ave. autosomal coverage genome",
    "Median autosomal coverage over genome" : "Median autosomal coverage genome",
    "Mean/Median autosomal coverage ratio over genome" : "Mean/Median autosomal coverage ratio",
    "Aligned reads" : "Aligned reads",
    "Aligned reads in genome" : "% Aligned reads genome",
    }

#coverage_exome
QC_COV_HEADINGS = {
    "Sample" : "Sample",
    "version" : "Version",
    "Work Order" : "Work Order",
    'Aligned bases':  'Aligned bases', 
    'Aligned bases in QC coverage region': 'Aligned bases in QC region',
                'Aligned bases in QC region' : '% Aligned bases in QC region', 
    'Average alignment coverage over QC coverage region':  'Average alignment coverage',
    'Uniformity of coverage (PCT > 0.2*mean) over QC coverage region': 'Uniformity of coverage (PCT > 0.2*mean)',
    'Uniformity of coverage (PCT > 0.4*mean) over QC coverage region': 'Uniformity of coverage (PCT > 0.4*mean)',
    'PCT of QC coverage region with coverage [1500x: inf)': 'PCT of QC region [1500X: inf)',
    'PCT of QC coverage region with coverage [1000x: inf)': 'PCT of QC region [1000X: inf)',
    'PCT of QC coverage region with coverage [ 500x: inf)': 'PCT of QC region [ 500X: inf)',
    'PCT of QC coverage region with coverage [ 100x: inf)': 'PCT of QC region [ 100X: inf)',
    'PCT of QC coverage region with coverage [  50x: inf)': 'PCT of QC region [  50X: inf)',
    'PCT of QC coverage region with coverage [  20x: inf)': 'PCT of QC region [  20X: inf)',
    'PCT of QC coverage region with coverage [  15x: inf)': 'PCT of QC region [  15X: inf)',
    'PCT of QC coverage region with coverage [  10x: inf)': 'PCT of QC region [  10X: inf)',
    'PCT of QC coverage region with coverage [   3x: inf)': 'PCT of QC region [   3X: inf)',
    'PCT of QC coverage region with coverage [   1x: inf)': 'PCT of QC region [   1X: inf)',
    'PCT of QC coverage region with coverage [   0x: inf)': 'PCT of QC region [   0X: inf)',
    'Average chr X coverage over QC coverage region': 'Ave. chrX coverage QC region',
    'Average chr Y coverage over QC coverage region': 'Ave. chrY coverage QC region',
    'Average mitochondrial coverage over QC coverage region': 'Ave. chrM coverage QC region',
    'Average autosomal coverage over QC coverage region': 'Ave. autosomal coverage QC region',
    'Median autosomal coverage over QC coverage region': 'Median autosomal coverage QC region',
    'Mean/Median autosomal coverage ratio over QC coverage region': 'Mean/Median autosomal coverage ratio',
    'Aligned reads': 'Aligned reads',
    'Aligned reads in QC coverage region': 'Aligned reads in QC region',
                'Aligned reads in QC region' : '% On Target',
    }

MAPPING_HEADINGS_v2 = {
    "Sample" : "Sample",
    "version" : "Version",
    "Work Order" : "Work Order",
    "Total input reads" : "Total reads",
                "Total reads" : "%Total reads",
    "Number of duplicate marked reads" : "Duplicate reads",
                "Duplicate reads" : "%Duplicate reads",
    "Number of duplicate marked and mate reads removed" : "Duplicates removed",
    "Number of unique reads (excl. duplicate marked reads)" : "Unique reads",
                "Unique reads" : "%Unique reads",
    "Reads with mate sequenced" : "Reads w/mate",
                "Reads w/mate" : "%Reads w/mate",
    "Reads without mate sequenced" : "Reads wo/mate",
                "Reads wo/mate" : "%Reads wo/mate", 
    "QC-failed reads" : "QC failed reads",
                "QC failed reads" : "%QC-failed reads",
    "Mapped reads" : "Mapped-reads",
                "Mapped-reads" : "%Mapped-reads",
    "Mapped reads adjusted for filtered mapping" : "Mapped reads (adjusted)",
                "Mapped reads (adjusted)" : "%Mapped adjusted",
    "Mapped reads R1" : "Mapped R1",
                "Mapped R1" : "%Mapped R1",
    "Mapped reads R2" : "Mapped R2",
                "Mapped R2" : "%Mapped R2", 
    "Number of unique & mapped reads (excl. duplicate marked reads)" : "Unique mapped reads",
                "Unique mapped reads" : "%Unique mapped reads",
    "Unmapped reads" : "Unmapped-reads",
                "Unmapped-reads" : "%Unmapped Reads",
    "Unmapped reads adjusted for filtered mapping" : "Unmapped reads (adjusted)",
                "Unmapped reads (adjusted)" : "%Unmapped adjusted",
    "Adjustment of reads matching non-reference decoys" : "Reads matching non-ref decoys",
                "Reads matching non-ref decoys" : "%Reads matching non-ref decoys",
    "Singleton reads (itself mapped; mate unmapped)" : "Mapped singletons",
                "Mapped singletons" : "%Mapped singletons",
    "Paired reads (itself & mate mapped)" : "Read pairs mapped",
                "Read pairs mapped" : "%Read pairs mapped",
    "Properly paired reads" : "Properly Paired reads",
                "Properly Paired reads" : "%Properly paired reads",
    "Not properly paired reads (discordant)" : "Discordant pairs",
                "Discordant pairs" :  "%Discordant pairs",
    "Paired reads mapped to different chromosomes" : "Pairs mapped to diff chr",
                "Pairs mapped to diff chr" : "%Pairs mapped to diff chr",
    "Paired reads mapped to different chromosomes (MAPQ>=10)" : "Pairs mapped to diff chr(MAPQ>=10)",
                "Pairs mapped to diff chr(MAPQ>=10)" : "%Pairs mapped to diff chr(MAPQ>=10)",
    "Reads with MAPQ [40:inf)" : "Reads w/ MAPQ [40:+)",
                "Reads w/ MAPQ [40:+)" : "%Reads w/ MAPQ [40:+)",
    "Reads with MAPQ [30:40)" : "Reads w/ MAPQ [30:40)",
                "Reads w/ MAPQ [30:40)" : "%Reads w/ MAPQ [30:40)",
    "Reads with MAPQ [20:30)" : "Reads w/ MAPQ [20:30)",
                "Reads w/ MAPQ [20:30)" : "%Reads w/ MAPQ [20:30)",
    "Reads with MAPQ [10:20)" : "Reads w/ MAPQ [10:20)",
                "Reads w/ MAPQ [10:20)" : "%Reads w/ MAPQ [10:20)",
    "Reads with MAPQ [ 0:10)" : "Reads w/ MAPQ [ 0:10)",
                "Reads w/ MAPQ [ 0:10)" : "%Reads w/ MAPQ [ 0:10)",
    "Reads with MAPQ NA (Unmapped reads)" : "MAPQ NA Unmapped reads",
                "MAPQ NA Unmapped reads" : "%MAPQ NA Unmapped reads",
    "Reads with indel R1" : "Reads w/indel R1",
                "Reads w/indel R1" : "%Reads w/indel R1",
    "Reads with indel R2" : "Reads w/indel R2",
                "Reads w/indel R2" : "%Reads w/indel R2",
    "Total bases" : "Total bases.",
    "Total bases R1" : "R1 bases",
    "Total bases R2" : "R2 bases",
    "Mapped bases" : "Mapped bases.",
    "Mapped bases R1" : "Mapped R1 bases",
    "Mapped bases R2" : "Mapped R2 bases",
    "Soft-clipped bases" : "S-clipped bases",
                "S-clipped bases" : "%S-clipped bases",
    "Soft-clipped bases R1" : "S-clipped bases R1",
                "S-clipped bases R1" : "%S-clipped bases R1",
    "Soft-clipped bases R2" : "S-clipped bases R2",
                "S-clipped bases R2" : "%S-clipped bases R2",
    "Hard-clipped bases" : "H-clipped bases",
                "H-clipped bases" : "%H-clipped bases", 
    "Hard-clipped bases R1" : "H-clipped bases R1",
                "H-clipped bases R1" : "%H-clipped bases R1",
    "Hard-clipped bases R2" : "H-clipped bases R2",
                "H-clipped bases R2" : "%H-clipped bases R2",
    "Mismatched bases R1" : "Mismatched-bases R1",
                "Mismatched-bases R1" : "%Mismatched-bases R1",
    "Mismatched bases R2" : "Mismatched-bases R2",
                "Mismatched-bases R2" : "%Mismatched-bases R2",
    "Mismatched bases R1 (excl. indels)" : "Mismatched bases R1 (wo/indels)",
                "Mismatched bases R1 (wo/indels)" : "%Mismatched bases R1 (wo/indels)",
    "Mismatched bases R2 (excl. indels)" : "Mismatched bases R2 (wo/indels)",
                "Mismatched bases R2 (wo/indels)" : "%Mismatched bases R2 (wo/indels)",
    "Q30 bases" : "Q30-bases",
                "Q30-bases" : "%Q30-bases", 
    "Q30 bases R1" : "Q30-bases R1",
                "Q30-bases R1" : "%Q30-bases R1",
    "Q30 bases R2" : "Q30-bases R2",
                "Q30-bases R2" : "%Q30-bases R2",
    "Q30 bases (excl. dups & clipped bases)" : "Q30 bases - dups & clipped bases)",
    "Total alignments" : "total alignments",
    "Secondary alignments" : "Sec. alignments",
    "Supplementary (chimeric) alignments" : "Chim. alignments",
    "Estimated read length" : "Est. read length",
    "Bases in reference genome" : "Bp in ref. genome",
    "Bases in target bed [% of genome]" : "Bases in target [% of genome]",
    "Average sequenced coverage over genome" : "Average coverage: genome",
    "Insert length: mean" : "Ins. length: mean",
    "Insert length: median" : "Ins. length: median",
    "Insert length: standard deviation" : "Insert length: stdev",
    "Provided sex chromosome ploidy" : "Provided sex chr ploidy",
    "Estimated sample contamination" : "Estimated contamination",
    "Estimated sample contamination standard error" : "Estimated cont. std error",
    "DRAGEN mapping rate [mil. reads/second]" : "DRAGEN mapping rate [M reads/sec]",
    }

COV_HEADING= {
    'Sample' : 'Sample',
    'version' : 'Version',
    'WorkOrder' : 'Work Order',
    'Alignedbases' : 'Aligned bases',
    'Alignedbasesingenome': '% Aligned bases genome',
    'Averagealignmentcoverageovergenome' : 'Average coverage genome',
    'Uniformityofcoverage(PCT>0.2*mean)overgenome' : 'Uniformity of cov genome (PCT > 0.2*mean)',
    'Uniformityofcoverage(PCT>0.4*mean)overgenome' : 'Uniformity of cov genome (PCT > 0.4*mean)',
    'PCTofgenomewithcoverage[1500x:inf)' : 'PCT of genome [1500x: inf)',
    'PCTofgenomewithcoverage[1000x:inf)' : 'PCT of genome [1000x: inf)',
    'PCTofgenomewithcoverage[500x:inf)' : 'PCT of genome [500x: inf)',
    'PCTofgenomewithcoverage[100x:inf)' : 'PCT of genome [100x: inf)',
    'PCTofgenomewithcoverage[50x:inf)' : 'PCT of genome [50x: inf)',
    'PCTofgenomewithcoverage[20x:inf)' : 'PCT of genome [20x: inf)',
    'PCTofgenomewithcoverage[15x:inf)' : 'PCT of genome [15x: inf)',
    'PCTofgenomewithcoverage[10x:inf)' : 'PCT of genome [10x: inf)',
    'PCTofgenomewithcoverage[3x:inf)' : 'PCT of genome [3x: inf)',
    'PCTofgenomewithcoverage[1x:inf)' : 'PCT of genome [1x: inf)',
    'PCTofgenomewithcoverage[0x:inf)' : 'PCT of genome [0x: inf)',
    'PCTofgenomewithcoverage[1000x:1500x)' : 'PCT of genome [1000x:1500x)',
    'PCTofgenomewithcoverage[500x:1000x)' : 'PCT of genome [500x:1000x)',
    'PCTofgenomewithcoverage[100x:500x)' : 'PCT of genome [100x:500x)',
    'PCTofgenomewithcoverage[50x:100x)' : 'PCT of genome [50x:100x)',
    'PCTofgenomewithcoverage[20x:50x)' : 'PCT of genome [20x:50x)',
    'PCTofgenomewithcoverage[15x:20x)' : 'PCT of genome [15x:20x)',
    'PCTofgenomewithcoverage[10x:15x)' : 'PCT of genome [10x:15x)',
    'PCTofgenomewithcoverage[3x:10x)' : 'PCT of genome [3x:10x)',
    'PCTofgenomewithcoverage[1x:3x)' : 'PCT of genome [1x:3x)',
    'PCTofgenomewithcoverage[0x:1x)' : 'PCT of genome [0x:1x)',
    'AveragechrXcoverageovergenome' : 'Ave. chrX coverage genome',
    'AveragechrYcoverageovergenome' : 'Ave. chrY coverage genome',
    'Averagemitochondrialcoverageovergenome' : 'Ave. chrM coverage genome',
    'Averageautosomalcoverageovergenome' : 'Ave. autosomal coverage genome',
    'Medianautosomalcoverageovergenome' : 'Median autosomal coverage genome',
    'Mean/Medianautosomalcoverageratioovergenome' : 'Mean/Median autosomal coverage ratio',
    'Alignedreads' : 'Aligned reads',
    'Alignedreadsingenome' : '% Aligned reads genome'
    }

QC_COV_HEADINGS_v2 = {
    'Sample' : 'Sample',
    'version' : 'Version',
    'WorkOrder' : 'Work Order',
    'Alignedbases' : 'Aligned bases',
    'AlignedbasesinQCcoverageregion' : 'Aligned bases in QC region',
    'Aligned bases in QC region' : '% Aligned bases in QC region',
    'AveragealignmentcoverageoverQCcoverageregion' : 'Average alignment coverage',
    'Uniformityofcoverage(PCT>0.2*mean)overQCcoverageregion' : 'Uniformity of coverage (PCT > 0.2*mean)',
    'Uniformityofcoverage(PCT>0.4*mean)overQCcoverageregion' : 'Uniformity of coverage (PCT > 0.4*mean)',
    'PCTofQCcoverageregionwithcoverage[1500x:inf)' : 'PCT of QC region [1500X: inf)',
    'PCTofQCcoverageregionwithcoverage[1000x:inf)' : 'PCT of QC region [1000X: inf)',
    'PCTofQCcoverageregionwithcoverage[500x:inf)' : 'PCT of QC region [ 500X: inf)',
    'PCTofQCcoverageregionwithcoverage[100x:inf)' : 'PCT of QC region [ 100X: inf)',
    'PCTofQCcoverageregionwithcoverage[50x:inf)' : 'PCT of QC region [  50X: inf)',
    'PCTofQCcoverageregionwithcoverage[20x:inf)' : 'PCT of QC region [  20X: inf)',
    'PCTofQCcoverageregionwithcoverage[15x:inf)' : 'PCT of QC region [  15X: inf)',
    'PCTofQCcoverageregionwithcoverage[10x:inf)' : 'PCT of QC region [  10X: inf)',
    'PCTofQCcoverageregionwithcoverage[3x:inf)' : 'PCT of QC region [   3X: inf)',
    'PCTofQCcoverageregionwithcoverage[1x:inf)' : 'PCT of QC region [   1X: inf)',
    'PCTofQCcoverageregionwithcoverage[0x:inf)' : 'PCT of QC region [   0X: inf)',
    'AveragechrXcoverageoverQCcoverageregion' : 'Ave. chrX coverage QC region',
    'AveragechrYcoverageoverQCcoverageregion' : 'Ave. chrY coverage QC region',
    'AveragemitochondrialcoverageoverQCcoverageregion' : 'Ave. chrM coverage QC region',
    'AverageautosomalcoverageoverQCcoverageregion' : 'Ave. autosomal coverage QC region',
    'MedianautosomalcoverageoverQCcoverageregion' : 'Median autosomal coverage QC region',
    'Mean/MedianautosomalcoverageratiooverQCcoverageregion' : 'Mean/Median autosomal coverage ratio',
    'Alignedreads' : 'Aligned reads',
    'AlignedreadsinQCcoverageregion' : 'Aligned reads in QC region',
    'Aligned reads in QC region' : '% On Target'
    }

#mapping_metrics
MAPPING_HEADINGS_v1 = {
    "Sample" : "Sample",
    "version" : "Version",
    "WorkOrder" : "Work Order",
    "Totalinputreads" : "Total reads",
        "Total reads" : "%Total reads",
    "Numberofduplicatemarkedreads" : "Duplicate reads",
        "Duplicate reads" : "%Duplicate reads",
    "Numberofduplicatemarkedandmatereadsremoved" : "Duplicates removed",
        "Duplicates removed" : "%Duplicates removed",
    "Numberofuniquereads(excl.duplicatemarkedreads)" : "Unique reads",
        "Unique reads" : "%Unique reads",
    "Readswithmatesequenced" : "Reads w/mate",
        "Reads w/mate" : "%Reads w/mate",
    "Readswithoutmatesequenced" : "Reads wo/mate",
        "Reads wo/mate" : "%Reads wo/mate",
    "QC-failedreads" : "QC failed reads",
        "QC failed reads" : "%QC-failed reads",
    "Mappedreads" : "Mapped-reads",
        "Mapped-reads" : "%Mapped-reads",
    "Mappedreadsadjustedforfilteredmapping" : "Mapped reads (adjusted)",
        "Mapped reads (adjusted)" : "%Mapped adjusted",
    "MappedreadsR1" : "Mapped R1",
        "Mapped R1" : "%Mapped R1",
    "MappedreadsR2" : "Mapped R2",
        "Mapped R2" : "%Mapped R2",
    "Numberofunique&mappedreads(excl.duplicatemarkedreads)" : "Unique mapped reads",
        "Unique mapped reads" : "%Unique mapped reads",
    "Unmappedreads" : "Unmapped-reads",
        "Unmapped-reads" : "%Unmapped Reads",
    "Unmappedreadsadjustedforfilteredmapping" : "Unmapped reads (adjusted)",
        "Unmapped reads (adjusted)" : "%Unmapped adjusted",
    "Adjustmentofreadsmatchingnon-referencedecoys" : "Reads matching non-ref decoys",
        "Reads matching non-ref decoys" : "%Reads matching non-ref decoys",
    "Singletonreads(itselfmapped;mateunmapped)" : "Mapped singletons",
        "Mapped singletons" : "%Mapped singletons",
    "Pairedreads(itself&matemapped)" : "Read pairs mapped",
        "Read pairs mapped" : "%Read pairs mapped",
    "Properlypairedreads" : "Properly Paired reads",
        "Properly Paired reads" : "%Properly paired reads",
    "Notproperlypairedreads(discordant)" : "Discordant pairs",
        "Discordant pairs" : "%Discordant pairs",
    "Pairedreadsmappedtodifferentchromosomes" : "Pairs mapped to diff chr",
        "Pairs mapped to diff chr" : "%Pairs mapped to diff chr",
    "Pairedreadsmappedtodifferentchromosomes(MAPQ>=10)" : "Pairs mapped to diff chr(MAPQ>=10)",
        "Pairs mapped to diff chr(MAPQ>=10)" : "%Pairs mapped to diff chr(MAPQ>=10)",
    "ReadswithMAPQ[40:inf)" : "Reads w/ MAPQ [40:+)",
        "Reads w/ MAPQ [40:+)" : "%Reads w/ MAPQ [40:+)",
    "ReadswithMAPQ[30:40)" : "Reads w/ MAPQ [30:40)",
        "Reads w/ MAPQ [30:40)" : "%Reads w/ MAPQ [30:40)",
    "ReadswithMAPQ[20:30)" : "Reads w/ MAPQ [20:30)",
        "Reads w/ MAPQ [20:30)" : "%Reads w/ MAPQ [20:30)",
    "ReadswithMAPQ[10:20)" : "Reads w/ MAPQ [10:20)",
        "Reads w/ MAPQ [10:20)" : "%Reads w/ MAPQ [10:20)",
    "ReadswithMAPQ[0:10)" : "Reads w/ MAPQ [ 0:10)",
        "Reads w/ MAPQ [ 0:10)" : "%Reads w/ MAPQ [ 0:10)",
    "ReadswithMAPQNA(Unmappedreads)" : "MAPQ NA Unmapped reads",
        "MAPQ NA Unmapped reads" : "%MAPQ NA Unmapped reads",
    "ReadswithindelR1" : "Reads w/indel R1",
        "Reads w/indel R1" : "%Reads w/indel R1",
    "ReadswithindelR2" : "Reads w/indel R2",
        "Reads w/indel R2" : "%Reads w/indel R2",
    "Totalbases" : "Total bases.",
    "TotalbasesR1" : "R1 bases",
    "TotalbasesR2" : "R2 bases",
    "Mappedbases" : "Mapped bases.",
    "MappedbasesR1" : "Mapped R1 bases",
    "MappedbasesR2" : "Mapped R2 bases",
    "Soft-clippedbases" : "S-clipped bases",
        "S-clipped bases" : "%S-clipped bases",
    "Soft-clippedbasesR1" : "S-clipped bases R1",
        "S-clipped bases R1" : "%S-clipped bases R1",
    "Soft-clippedbasesR2" : "S-clipped bases R2",
        "S-clipped bases R2" : "%S-clipped bases R2",
    "Hard-clippedbases" : "H-clipped bases",
        "H-clipped bases" : "%H-clipped bases",
    "Hard-clippedbasesR1" : "H-clipped bases R1",
        "H-clipped bases R1" : "%H-clipped bases R1",
    "Hard-clippedbasesR2" : "H-clipped bases R2",
        "H-clipped bases R2" : "%H-clipped bases R2",
    "MismatchedbasesR1" : "Mismatched-bases R1",
        "Mismatched-bases R1" : "%Mismatched-bases R1",
    "MismatchedbasesR2" : "Mismatched-bases R2",
        "Mismatched-bases R2" : "%Mismatched-bases R2",
    "MismatchedbasesR1(excl.indels)" : "Mismatched bases R1 (wo/indels)",
        "Mismatched bases R1 (wo/indels)" : "%Mismatched bases R1 (wo/indels)",
    "MismatchedbasesR2(excl.indels)" : "Mismatched bases R2 (wo/indels)",
        "Mismatched bases R2 (wo/indels)" : "%Mismatched bases R2 (wo/indels)",
    "Q30bases" : "Q30-bases",
        "Q30-bases" : "%Q30-bases",
    "Q30basesR1" : "Q30-bases R1",
        "Q30-bases R1" : "%Q30-bases R1",
    "Q30basesR2" : "Q30-bases R2",
        "Q30-bases R2" : "%Q30-bases R2",
    "Q30bases(excl.dups&clippedbases)" : "Q30 bases - dups & clipped bases)",
    "Totalalignments" : "total alignments",
    "Secondaryalignments" : "Sec. alignments",
    "Supplementary(chimeric)alignments" : "Chim. alignments",
    "Estimatedreadlength" : "Est. read length",
    "Basesinreferencegenome" : "Bp in ref. genome",
    "Basesintargetbed[%ofgenome]" : "Bases in target [% of genome]",
        "Bases in target [% of genome]" : "%Bases in target [% of genome]",
    "Averagesequencedcoverageovergenome" : "Average coverage: genome",
    "Insertlength:mean" : "Ins. length: mean",
    "Insertlength:median" : "Ins. length: median",
    "Insertlength:standarddeviation" : "Insert length: stdev",
    "Providedsexchromosomeploidy" : "Provided sex chr ploidy",
    "Estimatedsamplecontamination" : "Estimated contamination",
    "Estimatedsamplecontaminationstandarderror" : "Estimated cont. std error",
    "DRAGENmappingrate[mil.reads/second]" : "DRAGEN mapping rate [M reads/sec]",
    }

DICT_FIELDNAME_REF = {
        'genome_seq': {'MAPPING/ALIGNING SUMMARY': MAPPING_HEADINGS_v2,
        'COVERAGE SUMMARY' : COV_HEADING_v2}, 
            
        'exome_seq': {'MAPPING/ALIGNING SUMMARY' : MAPPING_HEADINGS_v2,
        'COVERAGE SUMMARY' : QC_COV_HEADINGS}
        }

XXXX_DICT_FIELDNAME_REF = {
        'genome_seq': {'MAPPING/ALIGNING SUMMARY': MAPPING_HEADINGS_v1,
        'COVERAGE SUMMARY' : COV_HEADING},

        'exome_seq': {'MAPPING/ALIGNING SUMMARY' : MAPPING_HEADINGS_v1,
        'COVERAGE SUMMARY' : QC_COV_HEADINGS_v2}
        }