#!/bin/sh
#

# Supporting functions, alphabetized...

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
chr_lis:,config:,config2:,\
output:,DP_low:,DP_high:,GQ_sample:,\
bin_size:,LEVEL:,version"

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
                    "CNV_detector"|"SNP_extractor"|"DSR_counter"|"SGR_PHR_definer")
                        [ "$2" = "CNV_detector" ] && printUsageCNV
                        [ "$2" = "SNP_extractor" ] && printUsageSNP
                        [ "$2" = "DSR_counter" ] && printUsageDSR
                        [ "$2" = "SGR_PHR_definer" ] && printUsageSGR
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
			--config) # -o outputFile or --output-file outputFile  Specify the output file
                config=$2
                shift 2  # Because output file is $2
                ;;
			--config2) # -o outputFile or --output-file outputFile  Specify the output file
                config2=$2
                shift 2  # Because output file is $2
                ;;
			--chr_lis) 
				chr_lis=$2
				shift 2
				;;
			--DP_low) # -i inputFile or --input-file inputFile  Specify the input file
				# Input file must be specified so $2 can be used
				DP_low=$2
				shift 2
				;;
			--DP_high) # -i inputFile or --input-file inputFile  Specify the input file
				# Input file must be specified so $2 can be used
				DP_high=$2
				shift 2
				;;
			--GQ_threshold) # -i inputFile or --input-file inputFile  Specify the input file
				# Input file must be specified so $2 can be used
				GQ_threshold=$2
				shift 2
				;;
			--bin_size) # -i inputFile or --input-file inputFile  Specify the input file
				# Input file must be specified so $2 can be used
				bin_size=$2
				shift 2
				;;
			--LEVEL) # -i inputFile or --input-file inputFile  Specify the input file
				# Input file must be specified so $2 can be used
				LEVEL=$2
				shift 2
				;;
			-v|--version) # -v or --version  Print the version
				printVersion
				exit 0
				;;
			--) # No more arguments , , , , 
				shift
                [ "$x" = "Pair_lis_generator" -o "$x" = "CNV_detector" -o "$x" = "SNP_extractor" -o "$x" = "DSR_counter" -o "$x" = "SGR_PHR_definer" ] && MODE=$1

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
    echo "Program: WheatComp (A pairwised comparison method to identify similar genetic regions and shared CNV regions.)
Version: 1.0"
    echo ""
	echo "Usage:   WheatComp [-v|--version] [-h|--help] <command> <argument>"
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

} 


printUsageCNV() {
	echo ""
    echo "About: detect CNV regions"
	echo "Usage: WheatComp CNV_detector <--chr_lis> [--single_CNV <--config1 <FILE>>] [--pair_CNV <--config2 <FILE>>]"
	echo ""
	echo "Options:"
    echo ""
    echo "   --chr_lis             chromosome list to detect CNV"
    echo ""
    echo "   --single_CNV          detect CNV region of samples listed"
    echo "                         produce three files: *.mask_CNV *.mask_CNV_deletion *.mask_CNV_duplication"
    echo "       --config <FILE>      file contains chromosome in BED file, path of BED file, path of BAM file and path of output directory"
    echo "                             e.g. chromosome BED BAM output_dev"
    echo ""
    echo "   --pair_CNV            identify sample specific CNV regions and shared CNV regions between sample pairs"
    echo "       --config2 <FILE>      file contains name and CNV file dircetory path of samples and output directory path"
    echo "                             e.g. sample1_name sample1_directory sample2_name sample2_directory pair_dircetory"
    echo ""
    echo "chromosome length in bed file shall not exceed the limit of bedtools (400Mb)"
    echo "see example/single_CNV.config and example/pair_CNV.config for more details (colunms separated by tab)"
}

printUsageSNP() {
	echo ""
    echo "About: extract SNP information from VCF for DSR calculation"
	echo "Usage: WheatComp SNP_extractor <--config <FILE>>"
	echo ""
	echo "Options:"
    echo ""
	echo "   --config <FILE>    file contains path of vcf files, sample ID list in vcf for extracting(separated by ','), path of output file"
    echo "                      e.g. VCF_file    vcf_ID   output_file"
    echo ""
    echo "see example/SNP_extractor.config for more details (colunms separated by tab)"
}

printUsageDSR() {
	echo ""
    echo "About: compute the DSR"
	echo "Usage: WheatComp DSR_counter <--config <file>>  [Options]"
    echo ""
    echo "   --config <FILE>      file contains path of files contain SNP inforamtion, path of output files and length of chromosomes (optional)"
    echo "                         e.g. SNP_file output_file length "
    echo ""
    echo "Options:"
    echo "   --DP_low <INT>            exclude SNP with DP <= <int>. default 3"
    echo ""
    echo "   --DP_high <INT>           exclude SNP with DP >= <int>. default 99"
    echo ""
    echo "   --GQ_sample <INT>         exclude SNP with GQ <= <int>. default 8"
    echo ""
    echo "   --bin_size <INT>          bln size. default 1000000"
    echo ""
    echo "see example/DSR_counter.config for more details (colunms separated by tab) "



}
printUsageSGR() {
	echo ""
    echo "About: detect SGR and PHR between sample pairs"
	echo "Usage: WheatComp DSR_counter <--config <file>> [--combine_CNV]"
	echo ""
	echo "Options:"
    echo ""
    echo "   --combine_CNV       combine CNV, SGR and PHR inforamtion between sample pairs together"
    echo ""
    echo "   --config <FILE>     file contains path of files contain DSR information and path of output files."
    echo "                       e.g. DSR_file  output"
    echo "                       if combines CNV information. config file shall contains path of CNV files and another path of output files."
    echo "                       e.g. DSR_file  output1 CNV_file output2"
}
# Print the program version
# - calling code needs to exit with the appropriate status
printVersion() {
	echo ""
	echo "WheatComp 1.0 2021-04"
    echo "License Expat: The MIT license
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law."
	echo ""
}

# Initialize variables
DP_low=3
DP_high=99
GQ_threshold=8
bin_size=1000000
LEVEL="10,1000"
MODE="ALL"
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
fi
echo $input $output



if test $MODE = "CNV_detector" ;then    # call CNV
    if [[ -z "${chr_lis}" ]]; then
        echo "No chromosome list provided."; exit 1
    fi
    if test $single_CNV = "on";then
        if [[ -z "${config}" ]]; then
            echo "No config file for CNV detection provided."; exit 1
        fi
        sh bisample_compare_phase1.sh   -c ${config}  -l ${chr_lis}
    fi
    if test $pair_CNV = "on";then
        if [[ -z "${config2}" ]]; then
            echo "No config file for pairwise sample CNV detection provided."; exit 1
        fi
        sh bisample_compare_phase1_2.sh -c ${config2}  -l ${chr_lis}  
fi

if test $MODE = "SNP_extractor" ;then    # call CNV
    if [[ -z "${config}" ]]; then
        echo "No config file provided."; exit 1
    fi
    while read lines;do
        VCF=`echo $lines|awk -F '\t' '{print $1}'`
        SAM_ID=`echo $lines|awk -F '\t' '{print $2}'`
        output=`echo $lines|awk '{print $4}'`

        bcftools view -v snps -e 'MAF<0.01' --min-ac=1 -M2 -m2 ${VCF} -Ou| \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n' -s ${SAM_ID} > ${output}
    done < $config
fi

if test $MODE = "DSR_counter" ;then
    if [[ -z "$config" ]]; then
        echo "No config file provided."; exit 1
    fi
    while read lines;do
        input=`echo $lines|awk -F '\t' '{print $1}'`
        output=`echo $lines|awk '{print $2}'`
        length=`echo $lines|awk '{print $3}'`
        ./muti_func_snp_compare.py --DSR_count on -i ${input} --length ${$length} \
         --DP_low ${DP_low} --DP_high ${DP_high} --GQ_sample ${GQ_sample} -b ${bin_size} -o ${output}
    done < $config
fi


if test $MODE = "SGR_PHR_definer" ;then
    if [[ -z "$config" ]]; then
        echo "No config file provided."; exit 1
    fi
    while read lines;do
        input=`echo $lines|awk '{print $1}'`
        output1=`echo $lines|awk '{print $2}'`
        ./muti_func_snp_compare.py --two_diff_level on -i ${input}  --LEVEL ${LEVEL} -o ${output}
        if test $combine_CNV = "on" ;then
            CNV=`echo $lines|awk '{print $3}'`
            output2=`echo $lines|awk '{print $4}'`
            line=`cat ${CNV}|wc -l`
            if test $line = 0;then
                awk '{print $0}' ${output1} > $output2
            else
                awk -F '\t' 'NR==FNR{a[$1]=$0;next}NR>FNR{if ($1 in a==1) {print a[$1]} else {print $0}}' ${CNV} ${output1} > ${output2}
            fi
        fi
    done < $config
fi
exit 0