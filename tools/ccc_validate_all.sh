#!/bin/bash
#==============================================================================
# For each input metrics txt file, create a correponding file of metrics
# values "manually" calculated from the original file.
# Attempt to invoke parallel nodes with mpirun to make faster.
#==============================================================================

function process_files_simple
{
  local num_way="$1"
  shift
  local snpbinfile="$1"
  shift
  local lifile="$1"
  shift
  local metricstxtfiles="$*"

  local metricstxtfile
  for metricstxtfile in $metricstxtfiles ; do

    if [ ! -e $metricstxtfile ] ; then
      echo "Error: file $metricstxtfile does not exist." 1>&2
      exit 1
    fi

    local outfile
    outfile=$(echo $metricstxtfile | sed -e 's/\.txt.*$/_validate&/')

    cut -f1-$(( 2 * $num_way )) -d' ' <$metricstxtfile \
      | ccc_validate $num_way $snpbinfile $lifile > $outfile

    echo COMPLETED $metricstxtfile

  done
}

#------------------------------------------------------------------------------

function main
{
  if [ "$*" = "" ] ; then
    echo "Usage: ${0##*/} <num_way> <snpbinfile> <line_index_file> <metrics_txt_file> ..."
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
        | tail -n +3 \
        | awk 'NR%'${OMPI_COMM_WORLD_SIZE}'=='${OMPI_COMM_WORLD_RANK})"
      local snpbinfile
      snpbinfile=$(echo $files \
        | tr ' ' '\12' \
        | sed -n -e '1p')
      local line_index_file
      line_index_file=$(echo $files \
        | tr ' ' '\12' \
        | sed -n -e '2p')
      # Process the files serially
      process_files_simple $num_way $snpbinfile $line_index_file $files_thisrank
    fi
  else # Not a rhea batch job ...
    # Process the files serially
    process_files_simple $num_way $files_spec
  fi
}

#------------------------------------------------------------------------------

main $@

#==============================================================================
