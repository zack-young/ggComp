#!/bin/sh
set -eo pipefail
#

# Parse the command line and set variables to control logic
parseCommandLine() {
    # Special case that nothing was provided on the command line so print usage
    # - include this if it is desired to print usage by default
    if [ "$#" -eq 0 ]; then
        printUsage
        exit 0
    fi

    optstring="h::i:o:v"

    optstringLong="help::,input:,\
chr_lis:,config:,config2:,single_CNV:,pair_CNV:,\
output:,DP_low:,DP_high:,GQ_sample:,no_CNV:,plus_CNV:,\
bin_size:,LEVEL:,version,folder_lis:,processes:,train,niter:"

    # - error message will be printed if unrecognized option or missing parameter but status will be 0

    GETOPT_OUT=$(getopt --options $optstring --longoptions $optstringLong -- "$@")

    exitCode=$?
    if [ $exitCode -ne 0 ]; then
        #echo ""
        printUsage
        exit 1
    fi

    eval set -- "$GETOPT_OUT"
    #echo $GETOPT_OUT
    # Loop over the options
    # - the error handling will catch cases were argument is missing
    # - shift over the known number of options/arguments
    while true; do
        #echo "Command line option is $opt"
        case "$1" in
            -h|--help) # -h or --help  Print usage
            case "$2" in
                "")
                    printUsage
                    exit 0
                    ;;
                "CNV_detector"|"SNP_extractor"|"DSR_counter"|"SGR_PHR_definer"|"HMM_smoother")
                    [ "$2" = "CNV_detector" ] && printUsageCNV
                    [ "$2" = "SNP_extractor" ] && printUsageSNP
                    [ "$2" = "DSR_counter" ] && printUsageDSR
                    [ "$2" = "SGR_PHR_definer" ] && printUsageSGR
                    [ "$2" = "HMM_smoother" ] && printUsageHMMS
                    exit 0
                    ;;
                *)
                    printUsage
                    exit 0
                    ;;
            esac
            ;;
            -i|--input) # -i inputFile or --input-file inputFile  Specify the input file
                # Input file must be specified so $2 can be used
                input=$2
                shift 2
                ;;
            -o|--output) # -o outputFile or --output-file outputFile  Specify the output file
                output=$2
                shift 2  # Because output file is $2
                ;;
            --config)
                config=$2
                shift 2  # Because output file is $2
                ;;
            --config2) 
                config2=$2
                shift 2  # Because output file is $2
                ;;
            --chr_lis) 
                chr_lis=$2
                shift 2  # Because output file is $2
                ;;
            --single_CNV) 
                SINGLE=true
                config=$2
                shift 2  # Because output file is $2
                ;;            
            --no_CNV) 
                NC=true
                config=$2
                shift 2
                ;;
            --plus_CNV) 
                PC=true
                config=$2
                shift 2  # Because output file is $2
                ;;            
            --pair_CNV) 
                PAIR=true
                config=$2
                shift 2
                ;;
            --DP_low)
                DP_low=$2
                shift 2
                ;;
            --DP_high)
                DP_high=$2
                shift 2
                ;;
            --GQ_threshold)
                GQ_threshold=$2
                shift 2
                ;;
            --bin_size)
                bin_size=$2
                shift 2
                ;;
            --LEVEL) 
                LEVEL=$2
                shift 2
                ;;
            --folder_lis)
                FolderNameList=$2
                shift 2
                ;;
            --processes)
                Procs=$2
                shift 2
                ;;
            --train)
                TrainFlag=true
                shift 1
                ;;
            --niter)
                Niter=$2
                shift 2
                ;;
            -v|--version) # -v or --version  Print the version
                printVersion
                exit 0
                ;;
            --) # No more arguments , , , , 
                shift
                [ "$1" = "Pair_lis_generator" -o "$1" = "CNV_detector" -o "$1" = "SNP_extractor" -o "$1" = "DSR_counter" -o "$1" = "SGR_PHR_definer" -o "$1" = "HMM_smoother" ] && MODE=$1
                shift
                break
                ;;
            *) # Unknown option - will never get here because getopt catches up front
                echo ""
                echo "Invalid option $1." >&2
                printUsage
                exit 1
                ;;
        esac
    done
    # Get a list of all command line options that do not correspond to dash options.
    # - These are "non-option" arguments.
    # - For example, one or more file or folder names that need to be processed.
    # - If multiple values, they will be delimited by spaces.
    # - Command line * will result in expansion to matching files and folders.
    shift $((OPTIND-1))
    additionalOpts=$*
}

# Print the program usage
# - calling code needs to exit with the appropriate status
printUsage() {
    echo ""
    echo "Program: ggComp (A pairwised comparison method to identify similar genetic regions and shared CNV regions between accessions.)
Version: 1.0"
    echo ""
    echo "Usage:   ggComp [-v|--version] [-h|--help] <command> <argument>"
    echo ""
    echo "Commands:"
    echo ""
    echo "    CNV_detector            detect CNV regions"
    echo ""
    echo "    SNP_extractor           extract SNP information from VCF"
    echo ""
    echo "    DSR_counter             compute the DSR"
    echo ""
    echo "    SGR_PHR_definer         identify SGR and PHR"
    echo ""
    echo "    HMM_smoother            smoothing SGR and PHR result using HMM"

} 


printUsageCNV() {
    echo ""
    echo "About: detect CNV regions"
    echo "Usage: ggComp CNV_detector <--chr_lis <STRING>>  <--single_CNV <FILE> | --pair_CNV <FILE>>"
    echo "   --chr_lis <STRING>    chromosomes that listed in behind file "
    echo "                         e.g. 'chr1A chr1B chr1D'"
    echo ""
    echo "Options:"
    echo ""
    echo "   --single_CNV          detect CNV region of samples ::: input file contains chromosme"
    echo "                         list, path of BED file, path of BAM file and output file"
    echo "                         e.g. chromosome BED BAM output_dev"
    echo "                         produce three files: *.mask_CNV"
    echo "                                              *.mask_CNV_deletion" 
    echo "                                              *.mask_CNV_duplication"
    echo ""
    echo "   --pair_CNV            identify pairwise sample specific CNV regions and shared CNV"
    echo "                         regions::: input file contains sample names, path of CNV file" 
    echo "                         and path of output file directory"
    echo "                         e.g. sample1_name sample1_directory sample2_name sample2_directory pair_dircetory"
    echo ""
    echo "chromosome length in bed file shall not exceed the limit of bedtools (400Mb)"
    echo "see test/single_CNV.config and example/pair_CNV.config for more details"
    echo " (colunms separated by tab)"
}

printUsageSNP() {
    echo ""
    echo "About: extract SNP information from VCF for DSR calculation"
    echo "Usage: ggComp SNP_extractor <--config <FILE>>"
    echo ""
    echo "Options:"
    echo ""
    echo "   --config <FILE>       file contains path of vcf files, sample ID list in vcf for"
    echo "                         extracting(separated by ','), path of output file"
    echo "                         e.g. VCF_file    vcf_ID   output_file"
    echo ""
    echo "see test/SNP_extractor.config for more details (colunms separated by tab)"
}

printUsageDSR() {
    echo ""
    echo "About: compute the DSR"
    echo "Usage: ggComp DSR_counter <--config <file>>  [Options]"
    echo ""
    echo "   --config <FILE>       file contains：path of files contain SNP inforamtion;" 
    echo "                         path of output files; end position of chromosomes" 
    echo "                         (use 'defalut' to call built-in end position of each chromosomes)"
    echo "                         e.g. SNP_file output_file END\default "
    echo ""
    echo "Options:"
    echo ""
    echo "   --DP_low <INT>        exclude SNP with DP <= <int>. default 3"
    echo ""
    echo "   --DP_high <INT>       exclude SNP with DP >= <int>. default 99"
    echo ""
    echo "   --GQ_sample <INT>     exclude SNP with GQ <= <int>. default 8"
    echo ""
    echo "   --bin_size <INT>      bln size. default 1000000"
    echo ""
    echo "see test/DSR_counter.config for more details (colunms separated by tab) "



}
printUsageSGR() {
    echo ""
    echo "About: detect SGR and PHR between sample pairs"
    echo "Usage: ggComp SGR_PHR_definer <--noCNV <file> | --plus_CNV <file>> [Options]"
    echo ""
    echo "   --no_CNV <FILE>      ignore CNV"
    echo "                        file contains path of files contain DSR information"    
    echo "                        and path of output files."
    echo "                        e.g. DSR_file  output"
    echo ""
    echo "   --plus_CNV <FILE>    take CNV into consideration"  
    echo "                        file contains path of files contain DSR information"    
    echo "                        path of CNV file and path of output files."
    echo "                        e.g. DSR_file  CNV_file output"
    echo ""
    echo "Options"
    echo ""
    echo "   --LEVEL <INT>       threshold divide SGR and PHR. default 10"
    echo ""
    echo "see test/SGR_PHR_definer.config for more details (colunms separated by tab) "

}

printUsageHMMS() {
    echo ""
    echo "About: smoothing SGR and PHR result using HMM"
    echo "Usage: ggComp HMM_smoother [Options]"
    echo ""
    echo "Options:"
    echo ""
    echo "   --input <FOLDER>      folder path of SGR and PHR phasing results"
    echo ""
    echo "   --output <FOLDER>     output path"
    echo ""
    echo "   --folder_lis <FILE>   containing a list of folder names of SGR and PHR phasing"
    echo "                         results that are goining to be smoothed"
    echo "                         e.g. C2_C10"
    echo "                              C2_C10"
    echo ""
    echo "   --processes <INT>     (optional) maximum worker processes"
    echo "                         default: 2"
    echo ""
    echo "   --train               (optional) train the model using input data (otherwise using"
    echo "                         the model published in the article)"
    echo ""
    echo "   --niter <INT>         (optional) maximum number of iterations to perform in training"
    echo "                         default: 60"
}

# Print the program version
# - calling code needs to exit with the appropriate status
printVersion() {
    echo ""
    echo "ggComp 1.0 2021-04"
    echo "License Expat: The MIT license
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law."
    echo ""
}

# Initialize variables
sp=$(cd "$(dirname "$0")";pwd)
DP_low=3
DP_high=99
GQ_sample=8
bin_size=1000000
LEVEL=10
MODE="ALL"
input=""
output=""
SINGLE=false
PAIR=false
PC=false
NC=false
FolderNameList=""
Procs=2
TrainFlag=false
Niter=60
# Parse the command line options
# - pass all arguments to the function
parseCommandLine "$@"
# Print input and output file
par_num=`echo $@ | awk -F '\t' '{print NF}'`
par_lis=$@
if test "$par_num" = 1 ;then 
[ "$par_lis" = "Pair_lis_generator" ] && printUsageLis && exit 0
[ "$par_lis" = "BED_generato" ] && printUsageBED && exit 0
[ "$par_lis" = "CNV_detector" ] && printUsageCNV && exit 0
[ "$par_lis" = "SNP_extractor" ] && printUsageSNP && exit 0
[ "$par_lis" = "DSR_counter" ] && printUsageDSR && exit 0
[ "$par_lis" = "SGR_PHR_definer" ] && printUsageSGR && exit 0
[ "$par_lis" = "HMM_smoother" ] && printUsageHMMS && exit 0
fi
#echo $input $output




if test $MODE = "CNV_detector" ;then    # call CNV
    if ($SINGLE) ;then
        if [[ -z "${config}" ]]; then
            echo "No config file for CNV detection provided."; exit 1
        fi
        if [[ -z "${chr_lis}" ]]; then
            echo "No chromosome list for CNV detection provided."; exit 1
        fi
        sh ${sp}/bisample_compare_phase1.sh -c ${config} -l ${chr_lis}
    fi
    if ($PAIR) ;then
        if [[ -z "${config}" ]]; then
            echo "No config file for pairwise CNV detection provided."; exit 1
        fi
        if [[ -z "${chr_lis}" ]]; then
            echo "No chromosome list for CNV detection provided."; exit 1
        fi
        sh ${sp}/bisample_compare_phase1_2.sh -c ${config} -l ${chr_lis}
    fi
fi

if test $MODE = "SNP_extractor" ;then    # call CNV
    if [[ -z "${config}" ]]; then
        echo "No config file provided."; exit 1
    fi
    if ! [ -x "$(command -v bcftools)" ]; then
        echo 'Error: bcftools is not installed or put in PATH.' >&2
        exit 1
    fi
    while read lines;do
        VCF=`echo $lines|awk  '{print $1}'`
        VCF_ID=`echo $lines|awk '{print $2}'`
        output=`echo $lines|awk '{print $3}'`

        bcftools view -v snps -e 'MAF<0.01' --min-ac=1 -M2 -m2 ${VCF} -Ou|\
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n' -s ${VCF_ID} -o ${output}
    done < $config
    echo "SNP extraction Complete"
fi

if test $MODE = "DSR_counter" ;then
    if [[ -z "$config" ]]; then
        echo "No config file provided."; exit 1
    fi
    while read lines;do
        input=`echo $lines|awk -F '\t' '{print $1}'`
        output=`echo $lines|awk '{print $2}'`
        length=`echo $lines|awk '{print $3}'`
        python ${sp}/muti_func_snp_compare.py --DSR_count on -i ${input} --length ${length} \
         --DP_low ${DP_low} --DP_high ${DP_high} --GQ_sample ${GQ_sample} -b ${bin_size} -o ${output}
    done < $config
fi


if test $MODE = "SGR_PHR_definer" ;then
    if ($NC) ;then
        while read lines;do
           input=`echo $lines|awk '{print $1}'`
           output=`echo $lines|awk '{print $2}'`
           python ${sp}/muti_func_snp_compare.py --two_diff_level on -i ${input}  --LEVEL ${LEVEL} -o ${output}
        done < $config
    fi
    if ($PC) ;then
        while read lines;do
            input=`echo $lines|awk '{print $1}'`
            CNV_file=`echo $lines|awk '{print $2}'`
            output=`echo $lines|awk '{print $3}'`
            python ${sp}/muti_func_snp_compare.py --two_diff_level on -i ${input}  --LEVEL ${LEVEL} -o ${output}_tmp
            line=`cat $CNV_file|wc -l`
            if test $line = 0;then
                mv ${output}_tmp  ${output}
            else
                awk -F '\t' 'NR==FNR{a[$1$2]=$0;next}NR>FNR{if ($1$2 in a==1) {print a[$1$2]} else {print $0}}' $CNV_file ${output}_tmp > ${output}
                rm -f ${output}_tmp
            fi
        done < $config
    fi
fi

if test $MODE = "HMM_smoother" ;then
    if [ ! -r  $FolderNameList ]; then
        echo "Invalid --folder_lis provided."; exit 1
    fi
    if [ ! -d  $input ]; then
        echo "Invalid --input provided."; exit 1
    fi
    if [ "${input:0:1}" != "/" ]; then
        echo "--input accepts absolute path only."; exit 1
    fi
    if [ ! -d  $output ]; then
        echo "Invalid --output provided."; exit 1
    fi
    if [ "${output:0:1}" != "/" ]; then
        echo "--output accepts absolute path only."; exit 1
    fi

    (>&2 echo -n "Checking environment ... ")
    python3 0_check_lib.py
    (>&2 echo "Done")

    current=`date "+%Y-%m-%d-%H:%M:%S"`
    (>&2 echo "Smoothing project: "$current)

    (>&2 echo -n "Building project ... ")
    python3 1_copy_datafile.py $input $output $FolderNameList $current $Procs
    (>&2 echo "Done")

    (>&2 echo -n "Pre-processing data ... ")
    python3 2_sort_and_pre-process.py $output $current $Procs
    (>&2 echo "Done")

    if ($TrainFlag); then
        (>&2 echo -n "Generating training data ... ")
        python3 3_make_train_data.py $output $current
        (>&2 echo "Done")

        (>&2 echo "Training ... ")
        python3 4_train.py $output $current $Niter > $output/$current/5.3_remake_levelfile.3.py
        (>&2 echo "Training ... Done")

        (>&2 echo -n "Smoothing ... ")
        cat 5.1_remake_levelfile.1.py $output/$current/5.3_remake_levelfile.3.py 5.2_remake_levelfile.2.py > $output/$current/5_remake_levelfile.py
        python3 $output/$current/5_remake_levelfile.py $output $current $Procs
        (>&2 echo "Done")
    else
        (>&2 echo -n "Smoothing ... ")
        python3 5_remake_levelfile.py $output $current $Procs
        (>&2 echo "Done")
    fi
fi

exit 0
