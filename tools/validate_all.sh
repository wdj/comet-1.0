#!/bin/bash
#==============================================================================
# Compare the validation files against the actual metrics files.
# Compare up a certain floating point tolerance.
# Can use parallel nodes with mpirun to make faster.
#==============================================================================

function process_files_simple
{
  local num_way="$1"
  shift
  local metricstxtfiles="$*"

  local metricstxtfile
  for metricstxtfile in $metricstxtfiles ; do

    if [ ! -e $metricstxtfile ] ; then
      echo "Checking $metricstxtfile ... Error: file does not exist."
      continue
    fi

    metricsvalidationfile=$(echo $metricstxtfile | sed -e 's/.txt$/_validate&/')

    if [ ! -e $metricsvalidationfile ] ; then
      echo "Checking $metricstxtfile ... Error: validation file does not exist."
      continue
    fi

    #echo "diff $F1 $F2"
    #diff <( sed -e 's/..$//' <$F1 ) <( sed -e 's/..$//' <$F2 )

    TOL=.00001
  
    F1=$metricstxtfile
    F2=$metricsvalidationfile
    C1=$(( 3 * $num_way + 1 ))
    C2=$(( 6 * $num_way + 2 ))

    echo -n "Checking $F1 ... " >/dev/null
    local numdiffs
    numdiffs=$(paste $F1 $F2 | awk '$'$C1' - $'$C2' > '$TOL' || $'$C2' - $'$C1' > '$TOL' {print $0 }' | wc -l)
    if [ $numdiffs = 0 ] ; then
      echo "PASSED." >/dev/null >/dev/null
    else
      echo "FAILED with $numdiffs diffs." >/dev/null
      echo "Checking $F1 ... FAILED with $numdiffs diffs."
    fi

  done
}

#------------------------------------------------------------------------------

function main
{
  if [ "$*" = "" ] ; then
    echo "Usage: ${0##*/} <num_way> <metrics_txt_file> ..."
    exit
  fi

  local num_way="$1"
  if [ "$num_way" -ne 2 -a "$num_way" -ne 3 ] ; then
    echo "Error: invalid value for num_way. $num_way" 1>&2
    exit 1
  fi

  shift
  local files_spec="$*"

  if [ -n "$PBS_NUM_NODES" ] ; then # If on rhea and in a batch job ...
    if [ -z "$OMPI_COMM_WORLD_SIZE" ] ; then # if not invoked by mpirun ...
      ppn=16
      tmpfile=tmp_$$
      # Store list of files in file because mpirun breaks command line breaks
      # if too many
      echo "$files_spec" > $tmpfile
      # Invoke mpirun, on this bash script
      mpirun -np $(( $PBS_NUM_NODES * $ppn )) --npernode $ppn $0 $num_way $tmpfile
      rm $tmpfile
    else # if invoked by mpirun ...
      files="$(cat $files_spec)" # retrieve the list of files
      # Let each mpi rank own a subset of the files
      local files_thisrank
      files_thisrank="$(echo $files \
        | tr ' ' '\12' \
        | awk 'NR%'${OMPI_COMM_WORLD_SIZE}'=='${OMPI_COMM_WORLD_RANK})"
      # Process the files serially
      process_files_simple $num_way $files_thisrank
    fi
  else # Not a rhea batch job ...
    # Process the files serially
    process_files_simple $num_way $files_spec
  fi
}

#------------------------------------------------------------------------------

main $@

#==============================================================================
