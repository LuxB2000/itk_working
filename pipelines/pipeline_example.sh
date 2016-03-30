#!/bin/bash - 
#===============================================================================
#
#          FILE: pipeline_example.sh
# 
#         USAGE: ./pipeline_example.sh 
# 
#   DESCRIPTION: 
#		WARNING: This pipeline is made for learning purpose. It doesn't compute
#		any actual concentration nor any 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ITK, Cmake
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Jerome Plumat, 
#  ORGANIZATION: University of Auckland
#       CREATED: 24/03/16 10:22
#      REVISION:  ---
#===============================================================================

set -o nounset        # Treat unset variables as an error

BINPATH=../bin
DATAPATH=../test_data
CONCENTRATIONPATH=$DATAPATH/concentrations
ATLASPATH=$DATAPATH/atlas

# inputs
# we assume that fixed and moving are normalized images. We put a phantom
# side to patient body, each pixel on a slice is divided by the average
# 
FIXEDIM=$DATAPATH/fixed_volume.mha
MOVINGIM1=$DATAPATH/moving_volume_t1.mha
MOVINGIM2=$DATAPATH/moving_volume_t2.mha
ATLASIM=$DATAPATH/atlas_volume.mha

# outputs
DEFORMEDIM1=$DATAPATH/deformed_volume1.mha
DEFORMEDIM2=$DATAPATH/deformed_volume2.mha
DEFORMEDATLAS=$ATLASPATH/deformed_atlas.mha
CONCVOLUME1=$CONCENTRATIONPATH/concentration_volume1.mha
CONCVOLUME2=$CONCENTRATIONPATH/concentration_volume2.mha
VOLUME4D=$CONCENTRATIONPATH/4D-volume.mha

# temporary data
TMPPATH=/tmp
DIRTRANSFORM=$TMPPATH/dir_transform.tfm
INVTRANSFORM=$TMPPATH/inv_transform.tfm
DEFFIELD=$TMPPATH/def_field.mha
NORMMOVING=$TMPPATH/norm_moving.mha
NORMFIXED=$TMPPATH/norm_fixed.mha

# Create the data
#----------------
mkdir $CONCENTRATIONPATH
mkdir $ATLASPATH

cd ../test_data
matlab -nojvm -nodisplay -nosplash -r create_data
cd ../pipelines


# Registrations
# -------------
# Do rigid registration between fixed and moving
./$BINPATH/MMRigid $FIXEDIM $MOVINGIM1 $DEFORMEDIM1 2 10 $DIRTRANSFORM $INVTRANSFORM
./$BINPATH/MMRigid $FIXEDIM $MOVINGIM2 $DEFORMEDIM2 2 10 $DIRTRANSFORM $INVTRANSFORM

# Do non rigid registration between fixed and atlas
./$BINPATH/MMRigid $FIXEDIM $ATLASIM $DEFORMEDATLAS 1 10 $DIRTRANSFORM $INVTRANSFORM 3 3 0
./$BINPATH/registerDemons $FIXEDIM $DEFORMEDATLAS $DEFORMEDATLAS $DEFFIELD

# Concentration
# -------------
# Do concentration estimation
# we assume that the fixed image is a pre-contrast agent injection image
# while the moving images are post-contrast agent injection
./$BINPATH/computeConcentrationNormIm $DEFORMEDIM1 $FIXEDIM $CONCVOLUME1
./$BINPATH/computeConcentrationNormIm $DEFORMEDIM2 $FIXEDIM $CONCVOLUME2

# Compute 4D volume
# -----------------
./$BINPATH/compute4dData $CONCENTRATIONPATH $VOLUME4D

# ROI measurements
# ----------------
# Extract the concentration in the different parts of the atlas
./$BINPATH/computeStat $VOLUME4D $ATLASPATH
