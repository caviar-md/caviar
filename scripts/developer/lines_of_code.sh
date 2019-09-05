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

cd ../..
echo ""
echo "============ *.h ===================="
echo "============ *.h ===================="
echo "============ *.h ===================="
echo "Number of lines for h header files"
find . -name '*.h' | xargs wc -l

echo ""
echo "============ *.cpp ===================="
echo "============ *.cpp ===================="
echo "============ *.cpp ===================="
echo "Number of lines for cpp source files"
find . -name '*.cpp' | xargs wc -l

echo ""
echo "============ *.cmake ===================="
echo "============ *.cmake ===================="
echo "============ *.cmake ===================="
echo "Number of lines for cmake script files"
find . -name '*.cmake' | xargs wc -l

echo ""
echo "============ CMakeLists.txt ===================="
echo "============ CMakeLists.txt ===================="
echo "============ CMakeLists.txt ===================="
echo "Number of lines for CMakeLists.txt script files"
find . -name 'CMakeLists.txt' | xargs wc -l

echo ""
echo "=========== *.casl ===================="
echo "=========== *.casl ===================="
echo "=========== *.casl ===================="
echo "Number of lines for caviar script files"
find . -name '*.casl' | xargs wc -l

# The SLOCCount tool may help as well.
