#!/bin/bash
################################################
# Clean out old results and link input files.
################################################

prefix='paris2c'

# Save current directory
ORIG_DIR=$PWD

# Optional argument: path to scripts/ directory.
# Can be useful if this script was called from a different directory.
# If unset, assume we're already in it
SCRIPT_DIR=$1
if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR=$ORIG_DIR
fi

cd $SCRIPT_DIR

# Empty the run directory - but first make sure it exists!
if [ -d "../run" ]; then
  cd ../run
  rm -rf *
else
  echo 'Creating run directory'
  mkdir ../run
  cd ../run
fi

# Copy everything from the input directory
cp ../input/bathymetry.shice .
cp ../input/data.cal .
cp ../input/data.diagnostics .
cp ../input/data.exf .
cp ../input/data.kpp .
cp ../input/data.pkg .
cp ../input/eedata .
cp ../input/eedata.mth .
cp ../input/shelfice_topo.bin .
cp ../input/data .
cp ../input/data.obcs .
cp ../input/data.shelfice .
cp ../input/pload_$prefix.mdjwf .
cp ../input/*_ini_$prefix.bin .
cp ../input/OB*_$prefix.bin .
# Remove prefixes
for f in OB*.bin; do mv "$f" "$(echo "$f" | sed s/_$prefix//)"; done
for f in *ini*.bin; do mv "$f" "$(echo "$f" | sed s/_$prefix//)"; done
for f in *.mdjwf; do mv "$f" "$(echo "$f" | sed s/_$prefix//)"; done

# Link executable
ln -s ../build/mitgcmuv .

cd $ORIG_DIR
