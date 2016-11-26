#!/bin/bash  
#--------------------------------------------------
# Script to make circos image.
# Please change the path for TABLEVIEWER_DIR
# and CIRCOS_DIR accordingly  
#
# Executed the following way:
#   $ make_circos_image.sh WORKING_DIR
#-------------------------------------------------- 
export LC_ALL=C

WORKING_DIR=$1
PERL=perl
CIRCOS_DIR=/u21/ewijaya/Tools/circos-0.64
BASENAME="image"
cd $WORKING_DIR
$PERL $CIRCOS_DIR/bin/circos -param random_string=$BASENAME -conf ./etc/circos-medium.conf
# $PERL $CIRCOS_DIR/bin/circos -param random_string=$BASENAME -conf ./etc/circos-clusterbased-medium.conf





