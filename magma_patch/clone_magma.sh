#!/bin/bash
#==============================================================================
#
# Create a custom version of MAGMA that avoids namespace collisions
# with other such versions.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  # Parse command line.

  local tag="${1:-}"
  if [ "$tag" = "" -o "$tag" = "-h" -o "$tag" = "--help" ] ; then 
    echo "Usage: ${0##*/} <tag>"
    echo "  where <tag> is an alphanumeric string. magma_<tag> is the new library name."
    exit 0
  fi

  local match
  match="`echo \"$tag\" | sed -e 's/[a-zA-Z0-9]//g'`"

  if [ "$match" != "" ] ; then
    echo "clone_magma: invalid argument." 1>&2
    exit 1
  fi

  # Set up input and output directory names.

  local magma_name="magma-1.6.2"
  # NOTE: if the MAGMA librsry version is updated, the lists below may also need
  # to be updated and other changes made, e.g, patch files.

  local source_file=${magma_name}.tar.gz
  if [ ! -e $source_file ] ; then
    echo "Error: unable to locate source file. $source_file" 1>&2
    exit 1
  fi

  local source_dir=$magma_name
  if [ -e $source_dir ] ; then
    echo "Error: source directory already exists. $source_dir" 1>&2
    exit 1
  fi

  local target_dir=magma_${tag}.cloned

  if [ -e $target_dir ] ; then
    echo "Error: target directory already exists. $source_dir" 1>&2
    exit 1
  fi

  # Un-tar the (real) MAGMA library.

  echo "Cloning $source_dir to $target_dir ..."

  echo "Unpacking ..."

  gunzip <$source_file | tar xf -

  match=$(grep -ri $tag $source_dir <(echo $tag) | wc -l)

  if [ $match != 1 ] ; then
    echo "Error: tag string already occurs in Magma source." 1>&2
    exit 1
  fi

  # Remove a few unneeded files from the library.

  rm $(find $source_dir -name '._*' -print)
  rmdir $source_dir/exp/lib $source_dir/exp/quark/lib
  rmdir $source_dir/testing/checkdiag/lib

  # Modify possibly conflicting strings.

  echo "Modifying name-colliding strings in each file ..."

  local file
  echo -n "magma"
  for file in `grep -ril 'magma' $source_dir` ; do
    sed -i "s/[mM][aA][gG][mM][aA]/&_$tag/g" $file
  done

  local names="
  lapack_bool_const
  lapack_order_const
  lapack_trans_const
  lapack_uplo_const
  lapack_diag_const
  lapack_side_const
  lapack_norm_const
  lapack_dist_const
  lapack_sym_const
  lapack_pack_const
  lapack_vec_const
  lapack_range_const
  lapack_vect_const
  lapack_direct_const
  lapack_storev_const
  lapacke_bool_const
  lapacke_order_const
  lapacke_trans_const
  lapacke_uplo_const
  lapacke_diag_const
  lapacke_side_const
  lapacke_norm_const
  lapacke_dist_const
  lapacke_sym_const
  lapacke_pack_const
  lapacke_vec_const
  lapacke_range_const
  lapacke_vect_const
  lapacke_direct_const
  lapacke_storev_const
  zhetrf_nopiv
  zpanel_to_q
  zq_to_panel
  chetrf_nopiv
  cpanel_to_q
  cq_to_panel
  dsytrf_nopiv
  dpanel_to_q
  dq_to_panel
  ssytrf_nopiv
  spanel_to_q
  sq_to_panel
  swp2pswp

  lapacke_const
  zgehrd_data
  cgehrd_data
  sgehrd_data
  dgehrd_data

  lapack_const
  cublas_trans_const
  cublas_uplo_const
  cublas_diag_const
  cublas_side_const

  zlaset_lower_kernel
  zlaset_upper_kernel
  zlaset_full_kernel
  "

  local name
  for name in $names ; do
    echo -n " $name"
    for file in $(grep -ril "$name" $source_dir) ; do
      sed -i "s/$name/&_$tag/g" $file
    done
  done

  echo

  # Modify possibly conflicting filenames.
  # NOTE: this is somewhat brittle, may break if ordering is changed.

  echo "Changing names of name-colliding files and directories ..."

  names="
  zhetrf_nopiv
  zpanel_to_q
  chetrf_nopiv
  cpanel_to_q
  dsytrf_nopiv
  dpanel_to_q
  ssytrf_nopiv
  spanel_to_q
  "

  for name in $names ; do
    echo -n "${name}:"
    local item
    for item in $(find $source_dir -name "${name}*" -print | sort | tac) ; do
      echo -n " $item"
      local new_item
      new_item=$(echo $item | sed -e "s/\(.*\)\(${name}\)\(.*\)/\1\2_$tag\3/g")
      mv $item $new_item
    done
    echo
  done

  echo -n "magma:"
  for item in $(find $source_dir -iname '*magma*' -print | sort | tail -n +2 | tac) ; do
    echo -n " $item"
    new_item=$(echo $item | sed -e "s/\(.*\)\([mM][aA][gG][mM][aA]\)\(.*\)/\1\2_$tag\3/g")
    mv $item $new_item
  done

  echo

  # Complete.

  mv $source_dir $target_dir

  echo "Clone of $source_dir to $target_dir complete."

#cat <<EOF
#
#If no errors occurred above, the clone process is now complete.  To
#prepare the cloned Magma code for usage, please do the following:
#
#1. Copy the appropriate make.inc file into the new cloned Magma directory.
#
#2. Make the desired changes to Magma, such as redefinition of certain
#   mathematical operations. See in particular gemm_stencil.cuh,
#   gemm_stencil_defs.h, dgemm_tesla_T_N.cu in the blas subdirectory.
#
#3. If desired, add/commit the new Magma clone to the repository.
#
#4. Compile the source code to deploy for usage.
#
#EOF

}

#==============================================================================

main "$@"

#==============================================================================
