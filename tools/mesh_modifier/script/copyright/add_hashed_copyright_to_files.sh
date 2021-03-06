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

x=$(find . -type f -name '*.cmake')
y=$(find . -type f -name '*.caviar')
for i in $x $y
do
  if ! grep -q Copyright $i
  then
    cat copyright_hashed.txt $i >$i.new && mv $i.new $i
  fi
done

