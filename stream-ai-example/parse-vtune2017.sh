#!/bin/bash 
# 
# Parses a VTune summary report for uncore memory access counts

usage () {
  echo "usage: $0 [-v | --verbose] [-s | --search_string name] file1 file2 ..."
  exit 0
}

# Initial parameter values
search="Uncore"
verbose=0

# Parse command line parameters
while [ "$#" -ge 1 ] ; do
  case "$1" in
    "-h" | "--help")           usage;;
    "-v" | "--verbose")        verbose=1; shift 1;;
    "-s" | "--search-string")  search=$2; shift 2;;
    *)                         break;;
  esac
done

# after parsing flags, make sure a file is specified
if [ $# -eq 0 ]; then
  usage
fi

# DDR
# UNC_M_CAS_COUNT.RD UNIT0 thru UNIT7
read_sum=([0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0)
# UNC_M_CAS_COUNT.WR UNIT0 thru UNIT7
write_sum=([0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0)

# MCDRAM
# UNC_E_RPQ_INSERTS UNIT0 thru UNIT7
mcd_read_sum=([0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0)
# UNC_E_WPQ_INSERTS UNIT0 thru UNIT7
mcd_write_sum=([0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0)

echo "Search stanza is \"$search\""
file_list=$@
files_parsed=0
mcdram="false";
for file in $file_list; do
  if [ -e $file ]; then
    exec < $file
    if [ $verbose -eq 1 ]; then echo "parsing $file"; fi
  else
    echo "***WARNING***: $file doesn't exist"
    continue
  fi
  # scan through file until the search string is found
  while read var1 var2; do
    if [ -n "$var2" ]; then 
      # search for match
      if [ "$var1" == "$search" ]; then
        if [ $verbose -eq 1 ]; then echo "$search found"; fi
        (( files_parsed++ ))
        # start reading lines looking for key values
        while read var1 var2; do
          # look for DDR counter
          unit=${var1#UNC_M_CAS_COUNT.}
          if [ $verbose -eq 1 ]; then echo "parsing $unit = $var2"; fi
          case "$unit" in
            "RD[UNIT0]")  (( read_sum[0] += $var2 ));;
            "RD[UNIT1]")  (( read_sum[1] += $var2 ));;
            "RD[UNIT2]")  (( read_sum[2] += $var2 ));;
            "RD[UNIT3]")  (( read_sum[3] += $var2 ));;
            "RD[UNIT4]")  (( read_sum[4] += $var2 ));;
            "RD[UNIT5]")  (( read_sum[5] += $var2 ));;
            "RD[UNIT6]")  (( read_sum[6] += $var2 ));;
            "RD[UNIT7]")  (( read_sum[7] += $var2 ));;
            "WR[UNIT0]")  (( write_sum[0] += $var2 ));;
            "WR[UNIT1]")  (( write_sum[1] += $var2 ));;
            "WR[UNIT2]")  (( write_sum[2] += $var2 ));;
            "WR[UNIT3]")  (( write_sum[3] += $var2 ));;
            "WR[UNIT4]")  (( write_sum[4] += $var2 ));;
            "WR[UNIT5]")  (( write_sum[5] += $var2 ));;
            "WR[UNIT6]")  (( write_sum[6] += $var2 ));;
            "WR[UNIT7]")  (( write_sum[7] += $var2 ));;
          esac
          # look for MCDRAM counter
          unit=${var1#UNC_E_RPQ_INSERTS}
          if [ $verbose -eq 1 ]; then echo "parsing $unit = $var2"; fi
          if [ $unit == "[UNIT0]" ]; then mcdram="true"; fi
          case "$unit" in
            "[UNIT0]")  (( mcd_read_sum[0] += $var2 ));;
            "[UNIT1]")  (( mcd_read_sum[1] += $var2 ));;
            "[UNIT2]")  (( mcd_read_sum[2] += $var2 ));;
            "[UNIT3]")  (( mcd_read_sum[3] += $var2 ));;
            "[UNIT4]")  (( mcd_read_sum[4] += $var2 ));;
            "[UNIT5]")  (( mcd_read_sum[5] += $var2 ));;
            "[UNIT6]")  (( mcd_read_sum[6] += $var2 ));;
            "[UNIT7]")  (( mcd_read_sum[7] += $var2 ));;
          esac
          unit=${var1#UNC_E_WPQ_INSERTS}
          if [ $verbose -eq 1 ]; then echo "parsing $unit = $var2"; fi
          case "$unit" in
            "[UNIT0]")  (( mcd_write_sum[0] += $var2 ));;
            "[UNIT1]")  (( mcd_write_sum[1] += $var2 ));;
            "[UNIT2]")  (( mcd_write_sum[2] += $var2 ));;
            "[UNIT3]")  (( mcd_write_sum[3] += $var2 ));;
            "[UNIT4]")  (( mcd_write_sum[4] += $var2 ));;
            "[UNIT5]")  (( mcd_write_sum[5] += $var2 ));;
            "[UNIT6]")  (( mcd_write_sum[6] += $var2 ));;
            "[UNIT7]")  (( mcd_write_sum[7] += $var2 ));;
          esac
        done
      fi 
    fi
  done
done

if [ $verbose -eq 1 ]; then
  echo "Total files parsed = $files_parsed"
fi

# print DDR out memory access summary
total_read=0
total_write=0
for i in "${!read_sum[@]}"; do
  echo "UNC_M_CAS_COUNT.RD[UNIT${i}] = ${read_sum[$i]}"
  (( total_read += (64 * read_sum[i]) ))
done
for i in "${!write_sum[@]}"; do
  echo "UNC_M_CAS_COUNT.WR[UNIT${i}] = ${write_sum[$i]}"
  (( total_write += (64 * write_sum[i]) ))
done
echo "--->DDR Report"
echo "--->Total Bytes read = $total_read"
echo "--->Total Bytes written = $total_write"
echo "--->Total Bytes = $(( total_read + total_write ))"
# print MCDRAM out memory access summary
if [ $mcdram == "true" ]; then
  total_read=0
  total_write=0
  for i in "${!mcd_read_sum[@]}"; do
    echo "UNC_E_RPQ_INSERTS[UNIT${i}] = ${mcd_read_sum[$i]}"
    (( total_read += (64 * mcd_read_sum[i]) ))
  done
  for i in "${!mcd_write_sum[@]}"; do
    echo "UNC_E_WPQ_INSERTS[UNIT${i}] = ${mcd_write_sum[$i]}"
    (( total_write += (64 * mcd_write_sum[i]) ))
  done
  echo "--->MCDRAM Report"
  echo "--->Total Bytes read = $total_read"
  echo "--->Total Bytes written = $total_write"
  echo "--->Total Bytes = $(( total_read + total_write ))"
fi
