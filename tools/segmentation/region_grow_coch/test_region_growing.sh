#!/bin/bash - 
#===============================================================================
#
#          FILE: test_level_set.sh
# 
#         USAGE: ./test_level_set.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 24/04/15 10:30
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

BINTOOLS=../../../bin

INPUTIM=/home/plumat/Workspace/Data/guineapigs_MRI/GP3/20100628-GPCochleaLPSDay4/3D-data/vol-extracted/22-LC-only.mha
COCHMESH=/tmp/coch-extracted.vtk
BM=/tmp/coch-segmented.mha
SEEDX=40
SEEDY=22
SEEDZ=12
SIGMA=0.05 # 0.05
ALPHA=-1000 #-300
BETA=2000 #2000
UPPERTHRESH=0.3
LOWERTHRESH=0.0
RADIUS=1
SEEDPATH=/tmp/coordSeeds.txt
./$BINTOOLS/regionGrowCoch $INPUTIM $BM $SIGMA $ALPHA $BETA $UPPERTHRESH $LOWERTHRESH $RADIUS $SEEDPATH
./$BINTOOLS/transformMask $BM $COCHMESH 2 1


