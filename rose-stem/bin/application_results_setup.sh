#!/bin/bash

# *****************************COPYRIGHT*******************************
# (C) Crown copyright 2024 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

# Prepare results location
mkdir -p $TASK_OUTPUT_DIR/results/
if [ $LUSTRE_FILESYSTEM ]; then
  # Set Lustre striping to maximum for results (performance)
  lfs setstripe -c -1 $TASK_OUTPUT_DIR/results/
fi
