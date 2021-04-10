#!/bin/bash - 
#===============================================================================
#
# FILE: shell_pack.sh
# 
# USAGE: ./shell_pack.sh 
# 
# DESCRIPTION: 
# 
# OPTIONS: ---
# REQUIREMENTS: ---
# BUGS: ---
# NOTES: ---
# AUTHOR: lwq (28120), scue@vip.qq.com
# ORGANIZATION: 
# CREATED: 04/22/2015 02:38:01 PM CST
# REVISION: ---
#===============================================================================

#=== FUNCTION ================================================================
# NAME: usage
# DESCRIPTION: Display usage information.
#===============================================================================
function usage ()
{
 cat <<- EOT

Usage : $0 -p package -s script file1 file2 file3 ..

Options:
 -h|help Display this message
 -p|package The output package name
 -s|script The script will run when unpack package
 Other The all files what you want to pack

EOT
} # ---------- end of function usage ----------

#-----------------------------------------------------------------------
# Handle command line arguments
#-----------------------------------------------------------------------
package_name="a"
install_script="b"
while getopts ":hp:s:" opt
do
 case $opt in

 h|help ) usage; exit 0 ;;
 p|package ) package_name=$OPTARG ;;
 s|script ) install_script=$OPTARG ;;
 \? ) echo -e "\n Option does not exist : $OPTARG\n"
 usage; exit 1 ;;

esac # --- end of case ---
done
shift $(($OPTIND-1))

echo ${package_name} ${install_script}

