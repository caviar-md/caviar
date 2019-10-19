
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

FOLDER_POV_MERGED='pov_merged/'
CAM_FILE='camera.pov'
MESH_FILE='o_mesh.pov'
SEARCH_DIR='../o_pov'

rm -r $FOLDER_POV_MERGED
mkdir $FOLDER_POV_MERGED

#============== merging files
echo "===== MERGING PART ====="

for ENTERY in "$SEARCH_DIR"/*
do
	FILE_NAME1=$(basename "$ENTERY")
	FILE_NAME2="m$FILE_NAME1"
#	echo "$FILE_NAME1"
#	echo "$FILE_NAME2"
#	echo "$FOLDER_POV_MERGED$FILE_NAME2" 
	cat $CAM_FILE $MESH_FILE $ENTERY > $FOLDER_POV_MERGED$FILE_NAME2
done

#============== png making
echo "===== PNG MAKING PART ====="

for ENTERY in "$FOLDER_POV_MERGED"/*
do
		povray $ENTERY +Q3 +W320 +H200
	#echo "$ENTERY"
done

#============== gif making
echo "===== GIF MAKING PART ====="

convert -delay 20 -loop 0 $FOLDER_POV_MERGED*.png myimage.gif

#==============
#eog myimage.gif

