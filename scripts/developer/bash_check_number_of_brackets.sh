#!/bin/bash

#=======================================================================
#
# Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
#
# This file is part of the CAVIAR package.
#
# The CAVIAR package is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 3.0 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CAVIAR distribution.
#
#=======================================================================

list=$(find . -type f)
for file in $list
do
  lbracket=$(grep -c '{' $file)
  rbracket=$(grep -c '}' $file)  
  echo $file  $lbracket $rbracket
  if [ $lbracket -ne $rbracket ]; then
    echo "WROOOOOOOOOOOOONG BRACKET"
  fi
done

    echo     
    echo "NAMESPACE"
    echo     
    
list=$(find . -type f)
for file in $list
do
  lbracket=$(grep -c '} // namespace caviar' $file)
  rbracket=$(grep -c 'namespace caviar {' $file)  
  echo $file  $lbracket $rbracket
  if [ $lbracket -ne $rbracket ]; then
    echo "WROOOOOOOOOOOOONG NAMESPACE"
  fi
done



