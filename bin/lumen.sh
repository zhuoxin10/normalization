#!/bin/sh
#
# Shell script to run lumen measurements for a set of challenges
#===============================================================

# bindir is dir where executables and awk-script reside
bindir=/clusterdata/tvanwalsum/challengeSrcs/bin

# refdir is the base directory for the reference data
refdir=$1

# ptc dir is the base dir for the participants data
ptcdir=$2

# out dir is the base output dir for the processed participants data
outdir=$3

# challs are the challenges to process, e.g.:
challs="000 001 002 003 004 005 006 007 008 100 101 102 200 201 202"

# debugging output or not?
debug=""
#debug="--debug"

# for the cluster
LANG=C
LC_ALL=C
export LANG LC_ALL

for i in $challs; do
  if  [[ -d $refdir/challenge$i ]] ; then
    echo === processing using $refdir/challenge$i;
  else
    echo === error: no reference dir $refdir/challenge$i;
    continue;
  fi;
  if  [[ -d $ptcdir/challenge$i ]] ; then
    echo === processing using $ptcdir/challenge$i;
  else
    echo === error: no participant dir $ptcdir/challenge$i;
    continue;
  fi;
  # preprocessing ...
  echo === processing challenge $i ...

  if mkdir -p $outdir/challenge$i; then
    echo === successfully created output dir
  else
    echo === error in output dir creation, continuing with next challenge
    continue;
  fi;

  echo === clipping partial volume to reference roi ...
  if $bindir/clip-pv $debug \
    --pv-r $refdir/challenge$i/pv$i \
    --roi-r $refdir/challenge$i/roi$i \
    --pv-o $ptcdir/challenge$i/pv$i \
    --roi-o $ptcdir/challenge$i/roi$i \
    --pv-out $outdir/challenge$i/pv-c$i; then
    echo === done.
  else
    echo === error in clip-pv, continuing with next challenge;
    continue;
  fi;

  echo === calculating normals for partial volume ...
  if $bindir/normals $debug --input $outdir/challenge$i/pv-c$i; then
    echo === done.
  else
    echo === error in clip-pv, continuing with next challenge;
    continue;
  fi;

  echo === determining sdm from partial volume ...
  if $bindir/pv2sdm $debug --pv $outdir/challenge$i/pv-c$i \
    --sdm $outdir/challenge$i/sdm$i; then
    echo === done.
  else
    echo === error in pv2sdm, continuing with next challenge;
    continue;
  fi;

  echo === determining iso from sdm ...
  if $bindir/sdm2iso $debug --sdm $outdir/challenge$i/sdm$i \
    --iso $outdir/challenge$i/iso$i; then
    echo === done.
  else
    echo === error in sdm2iso, continuing with next challenge;
    continue;
  fi;

  # measures ...
  echo === determining measures for challenge $i ...

  echo === dice ...
  if $bindir/dice $debug --mask $refdir/challenge$i/ext$i \
    --pv-r $refdir/challenge$i/pv$i \
    --pv-o $outdir/challenge$i/pv-c$i \
     > $outdir/challenge$i/dice.txt; then
    echo === done.
  else
    echo === error in dice, continuing with next challenge;
    continue;
  fi;

  echo === surfdist to reference iso ...
  if $bindir/surfdist $debug --mask  $refdir/challenge$i/ext$i \
    --iso $refdir/challenge$i/iso$i \
    --sdm $outdir/challenge$i/sdm$i \
    > $outdir/challenge$i/surfa.txt; then
    echo === done.
  else
    echo === error in surfdist a, continuing with next challenge;
    continue;
  fi;
  
  echo === surfdist from reference sdm ... 
  if $bindir/surfdist $debug --mask  $refdir/challenge$i/ext$i \
    --iso $outdir/challenge$i/iso$i \
    --sdm $refdir/challenge$i/sdm$i \
    > $outdir/challenge$i/surfb.txt; then
    echo === done.
  else
    echo === error in surfdist b, continuing with next challenge;
    continue;
  fi;

  echo === combining results ...
  cat $outdir/challenge$i/dice.txt \
    $outdir/challenge$i/surfa.txt \
    $outdir/challenge$i/surfb.txt \
       | awk -f $bindir/combine-lumen-results.awk \
       > $outdir/challenge$i/lumen-results$i.txt
  echo === done, results are \(see also $outdir/challenge$i/lumen-results$i.txt\):
  echo =====================================
  echo -n challenge$i: ; cat $outdir/challenge$i/lumen-results$i.txt 
  echo =====================================
done
