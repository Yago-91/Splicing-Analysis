# *Deep Splicing Analysis tool*

This is a python 3 software dedicated to finding splicing events. It allows novel splicing detection including intergenic jumps. It is called from the command line utility DSAT.sh. 




## Usage


    Arguments:

    -b|--bam            File with bam files names in new lines.

    -p|--Nthreads       Number of parallel threads to use.

    -s|--species        The species of bam files [mmusculus|hsapiens...]

    -c|--cols           Files previously processed by samtools
                        to get the splicing reads.

    -f|--config-file    Path to configuration file final formatting.

    -g|--gzip           If used, the intermediate files will be 
                        compressed.
    
    -d|--delete         If used, intermediate files will be removed
                        at the end of the process.

    --known             Pass a .gtf file to create the database or the 
                        file name of the database to use.             

    -o|--outdir         The directory to where save output.

    -h|--help           Prints this help.  
    
<br></br> 
 
## Required packages

Python >=3.5

**Python modules**  
argparse >=1.4  
pybiomart >=0.2  
multiprocessing  
scipy >=1.4
pandas >=0.25
numpy >=1.18  
pybedtools >=0.8   
pysam >=0.15  
glob 
time   
sys  
re

**Unix tools**  
samtools >=1.9


```python

```
