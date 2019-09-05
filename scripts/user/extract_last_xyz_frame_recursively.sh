#!/bin/bash
# this is a shell file that helps to run all the sbatch input files
#
#
# input.xyz:
# l 1   :
# l 2   :
# ...
# l x-2 : 
# l x-1 : $num_of_atoms
# l x   : Atom 
# l x+1 : type x y z ...
# ...
# l $num_of_lines : #last line of xyz file
# -----------------------------------------
#
# last_atom_line = x
# last_atom_line_minus_one = last_atom_line - 1
# We expect that:
# $num_of_atoms + (x-1) = $num_of_lines
# we check this.


MY_BATCH_FILE_INPUT='o_xyz.xyz'
MY_BATCH_FILE_OUTPUT='o_last_frame.xyz'



get_num_of_atoms_plus_two()
{
  last_id=$(($last_id-1))

  #the last line number containing 'Atom 
  last_atom_line=${tokenized_line_numbers[$last_id]}

  #the line before 'Atom' keyword containing number of atoms
  last_atom_line_minus_one=$(($last_atom_line-1))

  #get the number of xyz value line
  num_of_atoms=$(sed -n "$last_atom_line_minus_one"p $MY_BATCH_FILE_INPUT) 

  #add two lines because of number of atoms and 'Atom' keyword
  num_of_atoms_plus_two=$(($num_of_atoms+2))

  num_of_lines_using_data=$(($num_of_atoms+$last_atom_line))
}


RUN_FILES=$(find $PWD -type f -name $MY_BATCH_FILE_INPUT | sed -r 's|/[^/]+$||' |sort |uniq)
for i in $RUN_FILES
do
  #echo $i
  cd $i
  #ls

  num_of_lines=$(wc -l $MY_BATCH_FILE_INPUT)
  num_of_lines=( $num_of_lines )
  num_of_lines=${num_of_lines[0]}

  #line numbers containing 'Atom' keyword
  line_numbers_with_atom=$(grep -n "Atom" $MY_BATCH_FILE_INPUT | awk -F ":" '{print $1}' )


  tokenized_line_numbers=( $line_numbers_with_atom )

  tokenized_line_numbers_size=${#tokenized_line_numbers[@]}

  last_id=$tokenized_line_numbers_size

  get_num_of_atoms_plus_two

  # this is the case that the last frame is not complete
  if [ $num_of_lines -ne $num_of_lines_using_data ]
  then
    echo "Warning: incomplete last frame, i.e. num_of_lines:$num_of_lines -ne num_of_lines_using_data:$num_of_lines_using_data "
    echo "using one frame before the last"
    new_num_of_lines=$(($last_atom_line-2))
    get_num_of_atoms_plus_two

    if [ $new_num_of_lines -ne $num_of_lines_using_data ]
    then
      echo "ERROR: file $i/$MY_BATCH_FILE_INPUT is corrupted or the format is not supported"
      echo "in one frame before the last: new_num_of_lines:$new_num_of_lines -ne num_of_lines_using_data:$num_of_lines_using_data "
      continue
    fi
  fi

  #get the last frame and put it into another file
  tail -n $num_of_atoms_plus_two $MY_BATCH_FILE_INPUT > $MY_BATCH_FILE_OUTPUT
done



