#CSAT is a program for finding all the splicing events in a mapped sample (bam)
#and write the normalized output to a text file for their comparison among conditions
#Created by Ivó Hernández Hernández (ivohh91@gmail.com)

#-------Import required packages-------#
from argparse import ArgumentParser, FileType
from pybiomart import Server
import multiprocessing as mp
from scipy import stats
import pandas as pd
import numpy as np
import pybedtools
import pysam
import glob
import time
import sys
import re


#-------------------------------------------------Parse command line argumens----------------------------------------------------------------#
parser = ArgumentParser(prog='CSAT', description='Extract and normalize splice juntions from bam files')
parser.add_argument('--bam-txt', dest='bam_files', nargs='?', type=FileType('r'), help='Input file with bam file names in independent lines')
parser.add_argument('--species', dest='species', action='store', help='Species of the bam files [mmusculus|hsapiens|rnorvegicus...]')
parser.add_argument('--output', dest='output', type=str, action='store', help='path/to/output/file')
parser.add_argument('--config-file', dest='config', action='store', help='Path to configuration file')
parser.add_argument('--known-sites', dest='known', nargs='?', type=FileType('r'), help='Path to known_sites.db')
parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use in multiprocessing steps')

args = parser.parse_args()

#Control that the number of arguments is ok
if len(sys.argv) < 11:
    print('Insufficient number of arguments\n') 
    parser.print_help()
    exit(1)
#--------------------------------------------------------------------------------------------------------------------------------------------#

#--------------Build known splicing sites dictionary------------#
def build_known(known=args.known):
    '''This function creates a dicionary with 
    the jumps as keys and strand as values'''
    sites_known = {}

    for line in known:

        splitted = line.rstrip().split("\t")

        Chrom    = splitted[0] 
        take_off = int(splitted[1])
        landing  = int(splitted[2])
        strand   = int(splitted[3])

        composition = Chrom+'_'+str(take_off)+'_'+str(landing)

        sites_known[composition]  = strand
    
    return sites_known

known = build_known(args.known)
#----------------------------------------------------------------#

#Get matrix
server  = Server(host='http://www.ensembl.org')
dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets[args.species+'_gene_ensembl'])
result  = dataset.query(attributes = ['chromosome_name','start_position','end_position','strand'])

matrix = result.to_numpy()
#---------------------------------------------------#

######//From here the functions for bam processing are defined\\#######
#Function to get quality of the reads
def get_quality(string):
    lista=[]
    for i in string:
        lista.append(ord(i)-33)
    return np.mean(lista)

#Defining functions
def build_blocks(fetched):
    blocks = []
    for line in fetched:
        blocks.append(line.get_blocks())
    return blocks

#Function to count the number of IR reads
def count_IR(blocks,take_off,strand):
    count = 0
    for line in blocks:
        for h in range(len(line)):
            if strand == 1 or strand == 'default':
                if take_off in range(line[h][0]+1,line[h][1]+1) and take_off+1 in range(line[h][0]+1,line[h][1]+1):
                    count +=1
                    break
            else:
                if take_off-1 in range(line[h][0]+1,line[h][1]+1) and take_off in range(line[h][0]+1,line[h][1]+1):
                    count +=1
                    break

    return count

#Function to recurrently merge output files
def merge_sd(lista):
    nss = []
    for i in lista:
        nss.append(pd.read_csv(i,header=0,sep='\t', low_memory=False))
        
    for i in range(len(nss)-1):
        left=nss[i]
        right=nss[i+1]
    
        nss[i+1] = pd.merge(left, right, how='outer', on=['Coordinates'])
    return nss[-1]


#Parse CIGAR strings
def parse_CIGAR(start,cigar_arr):
    loc_N = np.where(np.array(cigar_arr) == 'N')
    result = np.tile((start,start),(np.shape(loc_N)[1],1))
    for h, loc in enumerate(np.nditer(loc_N)):
        take_off = 0
        for i in range(loc-1):
            if cigar_arr[i+1] == 'M':
                take_off += int(cigar_arr[i])
            if cigar_arr[i+1] == 'I':
                take_off += int(cigar_arr[i])
            if cigar_arr[i+1] == 'N':
                take_off += int(cigar_arr[i])
            if cigar_arr[i+1] == 'D':
                take_off += int(cigar_arr[i])
        result[h][0] += take_off-1
        result[h][1] += take_off+int(cigar_arr[loc-1])
    
    return result

#Function to fill the dictionary
def fill_dictionary(known,jumps,Chrom,dicty,Avg_Q,matrix):
    
    for jump in jumps:
        
        take_off = jump[0]
        landing  = jump[1]
    
        full_comp = Chrom+'_'+str(take_off)+'_'+str(landing)

        if full_comp in known.keys():
            strand  = known[full_comp]
            ann = 'annotated'
        else: 
            ann = 'novel'
            try:
                strand  = matrix[(matrix[:,0] == Chrom) & (matrix[:,1] < take_off) & (matrix[:,2] > landing)][0,3]
            except:
                code = []
                if matrix[(matrix[:,0] == Chrom) & (matrix[:,1] < take_off) & (matrix[:,2] > take_off)].size != 0:
                    code.append(matrix[(matrix[:,0] == Chrom) & (matrix[:,1] < take_off) & (matrix[:,2] > take_off)][0,3])
                
                elif matrix[(matrix[:,0] == Chrom) & (matrix[:,1] < landing) & (matrix[:,2] > landing)].size != 0:
                    code.append(matrix[(matrix[:,0] == Chrom) & (matrix[:,1] < landing) & (matrix[:,2] > landing)][0,3])
                
                if len(np.unique(code)) == 1:
                    strand = code[0]
                else:
                    strand = 'default'


        if strand == -1:
            xchange  = landing
            landing  = take_off
            take_off = xchange
        
        composition = Chrom+'_'+str(take_off)

        if composition not in dicty:             
            dicty[composition]  = [{landing:1},0,{landing:[Avg_Q]},strand,[ann]] 

        else:
            if landing in dicty[composition][0].keys():
                dicty[composition][0][landing] += 1
                dicty[composition][2][landing].append(Avg_Q)
                
            else:
                dicty[composition][0][landing] = 1
                dicty[composition][2][landing] = [Avg_Q]
                dicty[composition][4].append(ann)


#Function to get bed files from joined dataframes
def bed_prep(joint):
    """This function prepares jumps dataframe for making
    the closest gene finding"""
    
    #Process the merged file to split the firs column
    new = joint['Coordinates'].str.split(pat=':',expand=True)
    new2 = new[1].str.split(pat='-', expand=True)
    joint.insert(1,'Chromosome',new[0])
    joint.insert(2,"5'-coord",new2[0])
    joint.insert(3,"3'-coord",new2[1])
    joint = joint.drop(labels='Coordinates', axis=1)
    
    #Insert a column with unique identifier with the format Id[number]
    Ids = []
    for i in range(1,len(joint['Chromosome'])+1):
        Ids.append('Id'+str(i))

    joint.insert(0,'Id',Ids)
    
    #Make bed files for the both the 5' coordinate and 3' coordinate
    take_off_bed = pybedtools.BedTool.from_dataframe(joint[['Chromosome',"5'-coord","5'-coord",'Id']]).sort()
    landing_bed = pybedtools.BedTool.from_dataframe(joint[['Chromosome',"3'-coord","3'-coord",'Id']]).sort()
    
    return joint, take_off_bed, landing_bed

def biomart_bed(species):
    """This function returns bed objects from biomart data
    for the exon endings and startings separately"""
    
    #Retrieve data from biomart
    server  = Server(host='http://www.ensembl.org')
    dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets[species+'_gene_ensembl'])
    result  = dataset.query(attributes = ['chromosome_name','exon_chrom_start','exon_chrom_end','external_gene_name','ensembl_gene_id'])
    result.columns = ['Chromosome', 'Exon_start','Exon_end','Gene_name', 'Ensmbl_ID']
    
    #Create both biomart df w/o duplicates (exon end and exon start)
    biomart_take_off = result[['Chromosome','Exon_end','Exon_end','Gene_name','Ensmbl_ID']].drop_duplicates()
    biomart_landing  = result[['Chromosome','Exon_start','Exon_start','Gene_name','Ensmbl_ID']].drop_duplicates()
    
    #Create biomart bed objects sorted (very important)
    biomart_take_off_bed = pybedtools.BedTool.from_dataframe(biomart_take_off).sort()
    biomart_landing_bed  = pybedtools.BedTool.from_dataframe(biomart_landing).sort()
    
    return biomart_take_off_bed, biomart_landing_bed


def closest_gene(joint, take_off_bed, landing_bed, biomart_take_off_bed, biomart_landing_bed):
    """This function retrieves the a csv with the complete list of
    splicing events and corresponding gene name and ID"""
    
    #Perform BedTool.closest
    closest_take_off = take_off_bed.closest(biomart_take_off_bed)
    closest_landing  = landing_bed.closest(biomart_landing_bed)
    
    #Convert bed to pandas dataframes
    closest_take_off_df = closest_take_off.to_dataframe(low_memory=False, header=None)
    closest_landing_df  = closest_landing.to_dataframe(low_memory=False, header=None)
    
    #Select only data of interest
    take_off_tomerge = closest_take_off_df[['name','thickEnd','itemRgb']]
    landing_tomerge  = closest_landing_df[['name','thickEnd','itemRgb']]
    
    #Rename the columns
    take_off_tomerge.columns = ['Id','left_Gene','left_Ensembl_Id']
    landing_tomerge.columns  = ['Id','right_Gene','right_Ensembl_Id']
    
    #Perform merging of the closest objects
    merge_1 = pd.merge(joint,take_off_tomerge,how='outer', on=['Id'])
    merge_2 = pd.merge(merge_1,landing_tomerge,how='outer',on=['Id']).drop('Id', axis=1)
    
    #Reorder columns for final output
    order = merge_2.columns[0:3].append(merge_2.columns[-4:]).append(merge_2.columns[3:-4]).tolist()
    final = merge_2[order]
    
    return final


#--------------------------------------------------Define main function---------------------------------------------------------#
def main(file, known=known, matrix=matrix):
    '''Steps of the function:
        1-Create a dictionary of found jumps {chr_take-off:[{landing:count},All reads in the region,[Rel abundance of each jump]}
        2-Find reads covering each jump
        3-Print to output the normalized frequency of each jump'''

    #1---------------------------Creating dictionary------------------------------#
    #Open the file to be analyzed
    with open(file+'.cols', "r+") as cols:

        sites_sample = {}

        #print('Building splicing dictionary for file:', file,'\n')
        for site in cols:

            #First split the line by spaces
            splitted    = site.rstrip().split("\t")

            #Define variables with each field
            Chrom       = splitted[0]
            Start_site  = int(splitted[1])
            CIGAR       = splitted[2]
            Q_string    = splitted[3]

            #Break down CIGAR string and get relevant features
            CIGAR_split = re.findall(r"[^\W\d_]+|\d+", CIGAR)

            #Calculate quality
            Avg_Q       = get_quality(Q_string)

            jumps       = parse_CIGAR(Start_site,CIGAR_split)
            fill_dictionary(known=known,jumps=jumps, Chrom=Chrom, dicty=sites_sample, Avg_Q=Avg_Q, matrix=matrix)    


    #2---------------------Finding reads covering each jump site-----------------------
    ###Find total reads overlapping region of interest
    #print('Normalizing jumps:', file,'\n')
    bam = pysam.AlignmentFile(file, 'rb')

    
    #3---------------------Writing normalized frequency of each jump---------------------
    with open(file+".sd", "w+") as final:
        final.write('Coordinates'+'\t'+'Strand_'+file+'\t'+'Annotation_'+file+'\t'+'Splice_cts_'+file+'\t'+'IR_reads_'+file+'\t'+'Norm_cts_'+file+'\t'+'Avg_qual_'+file+'\n')
        for item in sites_sample.items():
            chrom    = item[0].split('_')[0]
            take_off = int(item[0].split('_')[1])
            strand   = item[1][3] 

            if strand == -1:
                All  = bam.fetch(chrom,take_off-1,take_off)
            else:
                All  = bam.fetch(chrom,take_off,take_off+1)

            blocks   = build_blocks(All)
            IR_reads = count_IR(blocks=blocks,take_off=take_off,strand=strand)

            val_sum = np.sum(list(sites_sample[item[0]][0].values()))
            percent = np.asarray(list(sites_sample[item[0]][0].values()))/(IR_reads+val_sum)*100

            for val, entry in enumerate(item[1][0].keys()):
                annotation = item[1][4][val]
                avg_q   = np.mean(sites_sample[item[0]][2][entry])
                final.write(chrom+':'+str(take_off)+'-'+str(entry)+'\t'+str(strand)+'\t'+annotation+'\t'+str(sites_sample[item[0]][0][entry])+'\t'+str(IR_reads)+'\t'+str(percent[val])+'\t'+str(avg_q)+'\n')
    
#----------------------------------------------------------------------------------------#
   
####Commencing main computation    
#print(str(time.localtime()[3])+'h:'+str(time.localtime()[4])+"':"+str(time.localtime()[5])+"''")

files = [line.rstrip() for line in args.bam_files]

#-----------Parallelize main function------------#
pool = mp.Pool(args.cores)
pool.map(main, [file for file in files])
pool.close()
#------------------------------------------------#

#Python under 3.5 won't search recursively using glob
print('Preparing final files...', str(time.localtime()[3])+'h:'+str(time.localtime()[4])+"':"+str(time.localtime()[5])+"''")

#Prepare biomart data
d,e = biomart_bed(args.species)

#Parse config file
config = pd.read_csv(args.config, header=None, low_memory=False, sep=' ')
config.columns = ['Sample_type','Sample_name', 'Condition']

conditions = np.unique(config['Condition']) 


#Arrange data to an adequate analysis format
def format_output(i,config=config,d=d,e=d,output=args.output):
    
    lista  = np.array(config[config.Condition == i]['Sample_name'])
    joint  = merge_sd(lista)
    a,b,c  = bed_prep(joint)
    interm = closest_gene(a,b,c,d,e)
    interm = interm.fillna(0)

    interm.rename(columns = {interm.columns[7]:'Strand'}, inplace=True)
    interm.rename(columns = {interm.columns[8]:'Annotation'}, inplace=True)

    for z in range(13,len(interm.columns)-1,6):
        interm.loc[interm['Strand'] == 0, 'Strand'] = interm[interm.columns[z]]
        interm.loc[interm['Annotation'] == 0, 'Annotation'] = interm[interm.columns[z+1]]
    interm.drop(interm.columns[13:-5:6], axis =1, inplace=True) #Erasing strand
    interm.drop(interm.columns[13:-4:5], axis =1, inplace=True) #Erasing annotation
    
    #Create summary table
    df2 = pd.DataFrame(data=interm.iloc[:,0:9])
    gene_list = np.unique(interm['left_Gene'])[1:]
    start = 9
    
    unique, index = np.unique(config[config.Condition == i]['Sample_type'], return_index=True)
    types = unique[index.argsort()]
    
    for h in types:
        n_samples_h = config[(config.Condition == i) & (config.Sample_type == h)].shape[0]
        df2['Avg_cts_'+h] = np.mean(interm.iloc[:,start:start+4*n_samples_h:4].to_numpy(),axis=1)
        df2['Avg_IR_'+h]  = np.mean(interm.iloc[:,start+1:start+1+4*n_samples_h:4].to_numpy(),axis=1)
        df2['Norm_abund_'+h] = np.mean(interm.iloc[:,start+2:start+2+4*n_samples_h:4].to_numpy(),axis=1)
        df2['SEM_Norm1_'+h] = stats.sem(interm.iloc[:,start+2:start+2+4*n_samples_h:4].to_numpy(),axis=1)
        df2['Avg_Q_'+h] = np.mean(interm.iloc[:,start+3:start+3+4*n_samples_h:4].to_numpy(),axis=1)
        start = start+4*n_samples_h
            
    #Insert Columns for Norm2
    for g, name in enumerate(lista):
        interm.insert(12+5*g, 'Norm2_'+name, np.zeros(interm.shape[0]))
            
    #Big step of second normalization
    gene_list = np.unique(interm['left_Gene'])[1:]
    
    for gene in gene_list:
        a = np.array(np.sum(interm.loc[interm.left_Gene == gene, interm.columns[9:-1:5]], axis=0))
        b = interm.loc[interm.left_Gene == gene, interm.columns[9:-1:5]].to_numpy()
        interm.loc[interm.left_Gene == gene, interm.columns[12:-1:5]] = b/a*100

    interm.fillna(0, inplace=True)

    #Write table with individual values
    interm.to_csv(output+'/'+i+'.sf', index=False) 

    #Complete df2 with Norm2
    start_a = 12
    start_b = 13

    for t in range(len(types)):
        n_samples_t = config[(config.Condition == i) & (config.Sample_type == types[t])].shape[0]
        to_insert_a = np.mean(interm.iloc[:,start_a:start_a+5*n_samples_t:5], axis=1)
        to_insert_b = stats.sem(interm.iloc[:,start_a:start_a+5*n_samples_t:5], axis=1)
        df2.insert(start_b+7*t, 'Norm2_'+types[t], to_insert_a)
        df2.insert(start_b+1+7*t, 'SEM_Norm2_'+types[t], to_insert_b)
        start_a = start_a+5*n_samples_t
    
    if len(types) == 2:
        df2['Log2_F.C._Norm1'] = np.log2((df2.iloc[:,18]+1)/(df2.iloc[:,11]+1))
        df2['DRA_Norm1'] = df2.iloc[:,18]-df2.iloc[:,11]
        df2['Log2_F.C._Norm2'] = np.log2((df2.iloc[:,20]+1)/(df2.iloc[:,13]+1))
        df2['DRA_Norm2'] = df2.iloc[:,20]-df2.iloc[:,13]

    df2.to_csv(output+'/'+i+'.sff', index=False)
#----------------------------------------------------------------------------------------------------#

#Format and write results
#-----------Parallelize format function------------#
pool = mp.Pool(args.cores)
pool.map(format_output, [condition for condition in conditions])
pool.close()
#------------------------------------------------#

#print(str(time.localtime()[3])+'h:'+str(time.localtime()[4])+"':"+str(time.localtime()[5])+"''")


