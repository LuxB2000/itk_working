#!/bin/bash - 
#===============================================================================
#
#          FILE: test_apply_transform_to_file_of_index.sh
# 
#         USAGE: ./test_apply_transform_to_file_of_index.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 30/04/15 12:17
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

BIN=../../bin

INPUTFILE=/tmp/coordSeeds.txt
MODIFCOORD=/tmp/coordSeedModif2.txt
INPUTVOLUME=/home/plumat/Workspace/Data/guineapigs_MRI/GP3/20100628-GPCochleaLPSDay4/3D-data/vol-extracted/22-LC-only.mha
#TODO: make a real registration atlas - image
ATLAS=$INPUTVOLUME
DTRANSF=/tmp/dtransf.tfm
ITRANS=/tmp/itransf.tfm
ATLASR=/tmp/atlasR.mha

#./$BIN/MMRigid $INPUTVOLUME $ATLAS $ATLASR 3 20 $DTRANSF $ITRANS
./$BIN/applyTransform2PointFile $INPUTFILE $DTRANSF $MODIFCOORD


