#!/usr/bin/env bash
set -eo pipefail

######################################################
# Input Path 
IP_flag=false
IP=""
# Folder List File
FL_flag=false
FL=""
# Output Path
OP_flag=false
OP=""
# Multi-processing
MP_flag=false
MP="2"
# Re-train
RT_flag=false
# Re-train: n_iter
NI_flag=false
NI=60
######################################################

(>&2 echo -n "Checking environment ... ")
python3 0_check_lib.py
(>&2 echo "Done")

while getopts ":i:f:o:p:tn:" optname; do
    case "$optname" in
        "i")
            IP_flag=true
            IP=$OPTARG
        ;;
        "f")
            FL_flag=true
            FL=$OPTARG
        ;;
        "o")
            OP_flag=true
            OP=$OPTARG
        ;;
        "p")
            MP_flag=true
            MP=$OPTARG
        ;;
        "t")
            RT_flag=true
        ;;
        "n")
            NI_flag=true
            NI=$OPTARG
        ;;
        ":")
            (>&2 echo "No argument value for option $OPTARG")
            exit
        ;;
        "?")
            (>&2 echo "Unknown option $OPTARG")
            exit
        ;;
        *)
            (>&2 echo "Unknown error while processing options")
            exit
        ;;
    esac
done

if ($IP_flag && $FL_flag && $OP_flag) ;then
    current=`date "+%Y-%m-%d-%H:%M:%S"`
    (>&2 echo "Project: "$current)

    (>&2 echo -n "Building project ... ")
    python3 1_copy_datafile.py $IP $OP $FL $current $MP
    (>&2 echo "Done")

    (>&2 echo -n "Pre-processing data ... ")
    python3 2_sort_and_pre-process.py $OP $current $MP
    (>&2 echo "Done")

    if ($RT_flag); then
        (>&2 echo -n "Generating training data ... ")
        python3 3_make_train_data.py $OP $current
        (>&2 echo "Done")

        (>&2 echo "Training ... ")
        python3 4_train.py $OP $current $NI > $OP/$current/5.3_remake_levelfile.3.py
        (>&2 echo "Training ... Done")

        (>&2 echo -n "Smoothing ... ")
        cat 5.1_remake_levelfile.1.py $OP/$current/5.3_remake_levelfile.3.py 5.2_remake_levelfile.2.py > $OP/$current/5_remake_levelfile.py
        python3 $OP/$current/5_remake_levelfile.py $OP $current $MP
        (>&2 echo "Done")
    else
        (>&2 echo -n "Smoothing ... ")
        python3 5_remake_levelfile.py $OP $current $MP
        (>&2 echo "Done")
    fi
else
    (>&2 echo "Key parameter (-i, -f, -o) lacked.")
    exit
fi

