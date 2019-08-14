#!/bin/bash -e

exe=../../cmake-build-release/minimize
outdir=output
mkdir -p $outdir


echo Minimize rec
name=$outdir/rec-min
$exe --prm parm.prm --rtf pdbamino.rtf --pdb BACE_4.rec.pdb --psf BACE_4.psf --out $name.pdb --nsteps 100 --json-log $name.json

echo Minimize rec 2 models
name=$outdir/rec2-min
$exe --prm parm.prm --rtf pdbamino.rtf --pdb BACE_4.rec.2models.pdb --psf BACE_4.psf --out $name.pdb --nsteps 100 --json-log $name.json

echo Rec+lig orig
name=$outdir/rec-lig-orig
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --out $name.pdb --nsteps 0

echo Score rec + lig
name=$outdir/rec-lig-score
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --out $name.pdb --score_only --json-log $name.json

echo Minimize rec + lig
name=$outdir/rec-lig-min
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --out $name.pdb --nsteps 1000 --json-log $name.json --fixed-pdb BACE_4.rec.pdb

echo Minimize rec + lig 2 models
name=$outdir/rec2-lig2-min
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.2models.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --lig-pdb BACE_4.lig.2modelsclose.pdb --out $name.pdb --nsteps 1000 --json-log $name.json --fixed-pdb BACE_4.rec.pdb

echo Minimize rec + lig 2 models and pointsprings
name=$outdir/rec2-lig2-pointspr
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.2models.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --lig-pdb BACE_4.lig.2modelsfar.pdb --out $name.pdb --nsteps 1000 --json-log $name.json --fixed-pdb BACE_4.rec.pdb --pointsprings pointsprings.txt

echo Minimize rec + lig 2 models and pairsprings
name=$outdir/rec2-lig2-pairspr
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.2models.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --lig-pdb  BACE_4.lig.2modelsfar.pdb --out $name.pdb --nsteps 1000 --json-log $name.json --fixed-pdb BACE_4.rec.pdb --pairsprings pairsprings.txt

echo Minimize rec + lig 2 models and density
name=$outdir/rec2-lig2-density
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.2models.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --lig-pdb  BACE_4.lig.2modelsclose.pdb --out $name.pdb --nsteps 1000 --json-log $name.json --fixed-pdb BACE_4.rec.pdb --fitting-pdblist <(echo density.pdb) --fitting-weight 10.0 

echo Minimize rec + lig 2 models and density and pointsprings-xl
name=$outdir/rec2-lig2-pointspr-density
$exe --prm parm.prm --rtf pdbamino.rtf --rec-pdb BACE_4.rec.2models.pdb --rec-psf BACE_4.psf --lig-json BACE_4.json --lig-pdb  BACE_4.lig.2modelsclose.pdb --out $name.pdb --nsteps 1000 --json-log $name.json --fixed-pdb BACE_4.rec.pdb --fitting-pdblist <(echo density.pdb) --fitting-weight 10.0 --pointsprings pointsprings-xl.txt



