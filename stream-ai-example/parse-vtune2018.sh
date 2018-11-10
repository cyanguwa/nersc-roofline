#!/bin/bash 
# 
# Parses a the output of a VTune version 2018 summary report for 
# uncore memory access counts using:
#
# amplxe-cl -report hw-events -r <result file> -group-by=package -format=csv -column=UNC_M_CAS_COUNT,UNC_E_RPQ_INSERTS,UNC_E_WPQ_INSERTS
#
# The above command outputs 2 lines of comma delimited column data, 
# the first line of a column is the name of the uncore counter field 
# and the second line of the column is the value

usage () {
  echo "usage: $0 [-v | --verbose] file1 file2 ..."
  exit 0
}

# Initial parameter values
verbose=0

# Parse command line parameters
while [ "$#" -ge 1 ] ; do
  case "$1" in
    "-h" | "--help")           usage;;
    "-v" | "--verbose")        verbose=1; shift 1;;
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

file_list=$@
files_parsed=0
for file in $file_list; do
  if [ -e $file ]; then
    exec < $file
    if [ $verbose -eq 1 ]; then echo "parsing $file"; fi
    (( ++files_parsed ))
  else
    echo "***WARNING***: $file doesn't exist"
    continue
  fi

  # build arrays with event name and counter values
  while read aline; do
    if [ $verbose -eq 1 ]; then echo "Parsing: $aline"; fi

    field0=${aline:0:7}
    case "$field0" in
      "Package" ) 
         IFS=',' read -ra events <<< "$aline"
         for ((i=1; i<${#events[@]}; ++i)); do
           # strip off useless text from front of string
           events[$i]=${events[$i]#Uncore Event Count:}
           if [ $verbose -eq 1 ]; then echo "Event $i is ${events[$i]}"; fi
         done ;;

      "package" ) 
        num=${aline:7:2}
        case "$num" in 
          "_0" )	# package_0
            IFS=',' read -ra counts0 <<< "$aline"
            for ((i=1; i<${#counts0[@]}; ++i)); do
              if [ $verbose -eq 1 ]; then echo "package_0 event $i is ${counts0[$i]}"; fi
            done ;;
          "_1" )	# package_1
            IFS=',' read -ra counts1 <<< "$aline"
            for ((i=1; i<${#counts1[@]}; ++i)); do
              if [ $verbose -eq 1 ]; then echo "package_1 event $i is ${counts1[$i]}"; fi
            done ;;
        esac ;;

      * )
        echo "WARNING: Unkown line quailfier -> $field0" ;;
    esac
  done # reading lines from a file

  # loop over all events
  for ((i=1; i<${#events[@]}; ++i)); do
    # check to see if more than a single socket with counter data
    if [ ${#counts1[@]} -eq 0 ]; then 
      value=${counts0[$i]}
    else
      value=$(( ${counts0[$i]} + ${counts1[$i]} ))
    fi

    if [ $verbose -eq 1 ]; then
      echo "processing ${events[$i]}, count = $value"
    fi

    # DDR uncore read counters
    if [ "${events[$i]:0:18}" == "UNC_M_CAS_COUNT.RD" ]; then
      unit=${events[$i]#UNC_M_CAS_COUNT.RD}
      case "$unit" in
        "[UNIT0]")  (( read_sum[0] += $value ));;
        "[UNIT1]")  (( read_sum[1] += $value ));;
        "[UNIT2]")  (( read_sum[2] += $value ));;
        "[UNIT3]")  (( read_sum[3] += $value ));;
        "[UNIT4]")  (( read_sum[4] += $value ));;
        "[UNIT5]")  (( read_sum[5] += $value ));;
        "[UNIT6]")  (( read_sum[6] += $value ));;
        "[UNIT7]")  (( read_sum[7] += $value ));;
      esac

    # DDR uncore write counters
    elif [ "${events[$i]:0:18}" == "UNC_M_CAS_COUNT.WR" ]; then
      unit=${events[$i]#UNC_M_CAS_COUNT.WR}
      case "$unit" in
        "[UNIT0]")  (( write_sum[0] += $value ));;
        "[UNIT1]")  (( write_sum[1] += $value ));;
        "[UNIT2]")  (( write_sum[2] += $value ));;
        "[UNIT3]")  (( write_sum[3] += $value ));;
        "[UNIT4]")  (( write_sum[4] += $value ));;
        "[UNIT5]")  (( write_sum[5] += $value ));;
        "[UNIT6]")  (( write_sum[6] += $value ));;
        "[UNIT7]")  (( write_sum[7] += $value ));;
      esac

    # MCDRAM uncore read counters
    elif [ "${events[$i]:0:17}" == "UNC_E_RPQ_INSERTS" ]; then
      mcdram=true;
      unit=${events[$i]#UNC_E_RPQ_INSERTS}
      case "$unit" in
        "[UNIT0]")  (( mcd_read_sum[0] += $value ));;
        "[UNIT1]")  (( mcd_read_sum[1] += $value ));;
        "[UNIT2]")  (( mcd_read_sum[2] += $value ));;
        "[UNIT3]")  (( mcd_read_sum[3] += $value ));;
        "[UNIT4]")  (( mcd_read_sum[4] += $value ));;
        "[UNIT5]")  (( mcd_read_sum[5] += $value ));;
        "[UNIT6]")  (( mcd_read_sum[6] += $value ));;
        "[UNIT7]")  (( mcd_read_sum[7] += $value ));;
      esac

    # MCDRAM uncore write counters
    elif [ "${events[$i]:0:17}" == "UNC_E_WPQ_INSERTS" ]; then
      mcdram=true;
      unit=${events[$i]#UNC_E_WPQ_INSERTS}
      case "$unit" in
        "[UNIT0]")  (( mcd_write_sum[0] += $value ));;
        "[UNIT1]")  (( mcd_write_sum[1] += $value ));;
        "[UNIT2]")  (( mcd_write_sum[2] += $value ));;
        "[UNIT3]")  (( mcd_write_sum[3] += $value ));;
        "[UNIT4]")  (( mcd_write_sum[4] += $value ));;
        "[UNIT5]")  (( mcd_write_sum[5] += $value ));;
        "[UNIT6]")  (( mcd_write_sum[6] += $value ));;
        "[UNIT7]")  (( mcd_write_sum[7] += $value ));;
      esac
    fi
  done # all events
done # all files

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
if [ $mcdram ]; then
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
