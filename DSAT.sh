#!/bin/bash

###The goal of this api is to detect known and novel splicing events
###1) Create known-splicing-db using selected GTF
###2) Retrieve spliced reads
###3) Detect known events and push to vast-tools2
###4) Detect unknown events and build table with main characteristics: Ens_ID-Tx_ID-IGV*coord

function help() {
    cat<<EOF
    Basic arguments:

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
    
EOF
}


while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--Nthreads)
    CORES="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--species)
    species="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--cols)
    COLS='YES'
    shift # past argument
    ;;
    -f|--config-file)
    config="$2"
    shift # past argument
    shift # past value
    ;;
    -g|--gzip)
    gzip='YES'
    shift # past argument
    ;;
    -d|--delete)
    Del='YES'
    shift # past argument
    ;;
    --known)
    known="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outdir)
    OUT_dir="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
    help
    exit 0
    ;;
    *)    # unknown option
    echo "Unrecognized option"; help; exit 1
    ;;

esac
done

#Reads with matching pattern
if [[ ! -d $OUT_dir ]]; then
    echo "Output directory does not exist"
    echo "Creating $OUT_dir directory..."$'\n'  
    mkdir $OUT_dir
fi      

parallel(){
  samtools view -@ $CORES $bam | grep '[0-9][0-9]*M[0-9][0-9]*N[0-9][0-9]*M' | cut -d $'\t' -f3,4,6,11 > $bam.cols
}

if [[ $known = *.gtf ]];
then    
    echo "Building known splicing sites database $(date +"%x %r %Z")"$'\n'
    python ~/bin/Splicing_analysis/Extract_known_splicing_sites.py $known > known_sites.tmp
    cat known_sites.tmp | sed 's/-/-1/g' | sed 's/+/1/g' > known_sites.db
    rm -fr known_sites.tmp
fi

IN_BAM=$(cat $BAM)

if [[ $COLS = 'YES' ]];
then
    echo '.cols files found'$'\n'
    echo "Started splicing analysis...$(date +"%r")"$'\n'
    python ~/bin/Splicing_analysis/DSAT.py --bam-txt $BAM --species $species --output $OUT_dir --config-file $config --known known_sites.db --cores $CORES
    echo "Finished splicing analysis...$(date + "%r")"$'\n'
else
    echo "Generating .cols files $(date +"%r")"$'\n'

    for bam in ${IN_BAM//,/ }
        do
            parallel &

        done

    wait

    echo "Samtools step completed $(date +"%r")"$'\n'
    echo "Started splicing analysis...$(date +"%r")"$'\n'
    python ~/bin/Splicing_analysis/DSAT.py --bam-txt $BAM --species $species --output $OUT_dir --config-file $config --known known_sites.db --cores $CORES
    echo "Finished splicing analysis...$(date +"%r")"$'\n'
fi

mv -f *.cols $OUT_dir
mv -f *.sd $OUT_dir

if [[ $gzip = 'YES' ]]; then
    for i in $OUT_dir/*.cols
        do
            pigz -7 -p $CORES $i
        
        done
fi

if [[ $Del = 'YES' ]]; then
    rm -fr $OUT_dir/*.cols
    rm -fr $OUT_dir/*.sd
fi

echo "End of the yourney $(date +"%x %r %Z")"