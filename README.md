# NRGmin 
Utility for molecular mechanics energy minimization

## Installation (Linux only)

1. Install [libmol2](https://bitbucket.org/bu-structure/libmol2/src/master) (verified to work with commit [90b174e](https://bitbucket.org/bu-structure/libmol2/commits/90b174ebae0523bed7eda11b0c3fce23adb9d376), Check version: 0.14, Jansson version: 2.12)

2. Build minimization executable 

        mkdir build && cd build
        cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DNOE=OFF ..
        make && make test
       
3. The binary is in `build/nrgmin`

4. To build OpenMP enabled binary, add `-DOPENMP=YES` to `cmake` command, it will create additional executable
   `build/nrgmin.omp`. `--num-threads` controls the number of OpenMP threads.
   
CMake flags:  
   
   * `BUILD_TESTS` — set to `YES` to enable tests, and to `NO` if you don't need them.
   
   * `USE_LTO` – set to `YES` for Link Time Optimization.
   
   * `USE_SANITIZER` — enable AddressSanitizer. Not recommended for release builds, since it complicates linking to other projects.
   
   * `OPENMP` — set to `YES` to create a multithreaded binary `nrgmin.omp` in addition to a serial one.
   
   * `NOE` — set to `NO` if you want to drop NOE (NMR) calculation capability (e.g. if libmol2 was build without NOE).

## How to run

Full list of parameters can be viewed by running  `./nrgmin -h`

### Passing a molecule

A molecule can be passed using single file mode as a set of 4 files

```
./nrgmin --pdb mol.pdb --psf mol.psf --prm mol.prm --rtf mol.rtf
```

Or as a json file

```
./nrgmin --json mol.json
```

As well as in a rec/lig mode:

```
./nrgmin --rec-pdb rec.pdb --rec-psf rec.psf --rec-prm rec.prm --rec-rtf rec.rtf --lig-json lig.json
```

If you want to use parameters from json file and coordinates from pdb

```
./nrgmin --rec-pdb rec.pdb --rec-json rec.json --lig-pdb lig.pdb --lig-json lig.json
```

#### Multiple models, multithreading

PDB file can contain multiple models. `nrgmin.omp` parallelizes minimization for multiple models 
using the flag `--num-threads`. All available cores are used by default.

```
./nrgmin.omp --json mol.json --pdb mol_10models.pdb --num-threads 10
```

In `rec/lig` mode the receptor can contain a single model and the ligand and can have N models 
and vice versa. In this case the receptor model is copied N times and merged with the ligand models. 
Otherwise, the number of receptor and ligand models must match. Therefore, the following will work:

```
./nrgmin.omp --rec-pdb rec_1model.pdb --lig-pdb lig_10models.pdb --rec-psf rec.psf --rec-prm rec.prm --rec-rtf rec.rtf --lig-json lig.json

./nrgmin.omp --rec-pdb rec_10models.pdb --lig-pdb lig_10models.pdb --rec-psf rec.psf --rec-prm rec.prm --rec-rtf rec.rtf --lig-json lig.json

./nrgmin.omp --rec-pdb rec_10models.pdb --lig-pdb lig_1model.pdb --rec-psf rec.psf --rec-prm rec.prm --rec-rtf rec.rtf --lig-json lig.json
```

And the following will not:

```
./nrgmin.omp --rec-pdb rec_3models.pdb --lig-pdb lig_10models.pdb --rec-psf rec.psf --rec-prm rec.prm --rec-rtf rec.rtf --lig-json lig.json
```

### Energy function

Energy function setup is flexible, so you can switch terms on and off

```
./nrgmin --json mol.json --vdw-off --ace-on
```

And add restraining potentials

```
./nrgmin --json mol.json --pair-springs-txt springs.txt
```

List of restraints:      

* fixed atoms (`--fix-receptor`, `--fix-ligand`, `--fixed-pdb`)
* pairwise distance restraints (`--pair-springs-txt`)
* pointwise distance restraints (`--point-springs-txt`)
* density fitting (`--density-json`)
* NOE matrix fitting (`--noe-txt`, `--noe-json`) 

Instead of providing flags you can setup all minimization parameters and 
multiple minimization stages using a master json file

```
./nrgmin --setup-json setup.json
```

To see all flags:
```
./nrgmin -h
```

### Freezing atoms
Atoms can be frozen in several ways:   

1. `--fix-receptor` and `--fix-ligand` flags can be used in `rec/lig` mode
   

2. `--fixed-pdb` flag can be used to pass fixed atoms. Atoms in the minimized 
   atom group will be frozen, if they are closer than 0.01 A to any of the atoms in
   `--fixed-pdb`. `--fixed-pdb` can contain multiple models, in which case the
   number of models should match the number of models in the minimized atom group. 
   If it contains only one model and the minimized atom group has several, 
   then the same `--fixed-pdb` will be applied to all the models in the minimized 
   atom group.


3. Field `fixed` in `--setup-json` can specify IDs of frozen atoms (for example `fixed: [0, 1, 3, 4]`)


### Master json file format for --setup-json

Has two fields: `options`, which duplicates all existing flags, and `stages`. `Stages` is an array, 
where each entry describes a minimization stage and can use all existing flags. Fields specified in 
`options` will be used throughout all stages, unless specified differently in each stage. So if an option 
is not specified in stage parameters, its value is taken from `options` or is set to default, if it is not 
set there as well.

```
{
  "options": {
    "rec_pdb": "BACE_4_rec_2models.pdb",
    "rec_psf": "BACE_4_rec.psf",
    "rec_prm": "BACE_4_rec_parm.prm",
    "rec_rtf": "BACE_4_rec_pdbamino.rtf",
    "lig_json": "BACE_4_lig.json",
    "lig_pdb": "BACE_4_lig_2models_close.pdb"
  },
  "stages": [
    {
      "nsteps": 1000,
      "fix_ligand": true
    },
    {
      "nsteps": 1000,
      "fix_receptor": true
    },
    {
      "nsteps": 1000,
      "vdw03": false,
      "ace": true,
      "pointsprings": [
        {
          "weight": 10,
          "coords": [
            26.816,
            1.943,
            23.790
          ],
          "atoms": [
            3598
          ]
        }
      ],
      "pairsprings": [
        {
          "distance": 10,
          "weight": 10,
          "lerror": 1,
          "rerror": 1,
          "potential": 2,
          "average": 0,
          "group1": [
            3598
          ],
          "group2": [
            3630
          ]
        }
      ]
    }
  ]
}
```

### Pairsprings setup dictionary

As can be seen in the example of Master json above, to use the "pairsprings" potential energy in minimization, one dictionary needed to be created for every distance restraint and 8 different parameters need to be assigned in the dictionary: distance, weight, lerror, rerror, potential, average, group1, group2.  The "distance" denotes expected distance between the two sets of atoms. The "lerror" and "rerror" are the left and right tolerance errors from the expected distance respectively. The "weight" is the scale factor for the potential energy function. The "potential " is the type of penalty function for the distance restraints. There are three different types of penalty function: Square-Well (0), Biharmonic (1), Soft-Square (2). The "average" denotes the method to average the distances between the two selected sets of atoms, and two types of averaging can be selected: SUM (0) and R-6 (1). The "group1" and "group2" are two lists of indices that indicate the two sets of atoms in the restraint. If there are more than one atom in either set, the distances of all possible pairs between two groups will be computed. Then the average distance of all the computed distances will be calculated using "average" method. This averaged distance is the calculated distance for this restraint. More details of the description of distance restraints can be found here https://nmr.cit.nih.gov/xplor-nih/xplorMan/node376.html.

### Examples

More setup examples are [here](./examples/energy_setup) and usage examples [here](./examples/usage).

### All flags

#### Input/Output

Output control:
        
            --out-pdb Where to write the minimized molecule(s)
            --out-json Log file with energy terms
            --out-prefix Overrides the above two options
            --print-step Log energies for every step
            --print-stage Log energies for every stage
            --print-noe-matrix Log NOE matrix is NOE calculation is on
            --verbosity 0 - QUIET, 1 - ERROR, 2 - WARNING, 3 - INFO, >=4 - DEBUG
        
        
Single file mode:
        
            --psf Geometry
            --prm Forcefield parameters
            --rtf Topology file with residues description
            --pdb Atom coordinates (can contain multiple models)
            --json Json file with coordinates, geometry and force field parameters
        
Parameters for rec/lig mode
        
            --rec-psf Receptor geometry
            --rec-prm Receptor forcefield parameters
            --rec-rtf Receptor topology file
            --rec-pdb Receptor atom coordinates (can contain multiple models)
            --rec-json Receptor json file with coordinates, geometry and force field parameters
            --lig-psf Ligand geometry
            --lig-prm Ligand forcefield parameters
            --lig-rtf Ligand topology file with residues description
            --lig-pdb Ligand atom coordinates (can contain multiple models)
            --lig-json Ligand json file with coordinates, geometry and force field parameters
        
        
#### Minimization setup
        
            --nsteps Number of minimization steps
        
Energy terms switches. Everything is on by default except for GBSA
        
            --bonds-on/--bonds-off Bonds energy term
            --angles-on/--angles-off Angles
            --dihedrals-on/--dihedrals-off Dihedrals
            --impropers-on/--impropers-on Impropers
            --vdw-on/--vdw-off VDW
            --vdw03-on/--vdw03-off 1-4 VDW
            --eleng-on/--eleng-off Coulomb electrostatics
            --elengs03-on/--elengs03-off 1-4 Coulomb electrostatics
            --gbsa-on/gbsa-off GBSA
        
Fixed atoms
        
            --fixed-pdb PDB file with fixed atoms (uses coordinates, not atom IDs)
            --fix-receptor Fix receptor atoms (works only in rec/lig mode)
            --fix-ligand Fix ligand atoms (works only in rec/lig mode)
        
Distance restraints
        
            --pair-springs-txt Text file setup for pairwise restraints
            --point-springs-txt Text file setup for pointwise restraints
        
NOE setup
        
            --noe-txt Text file setup for NOE
            --noe-json Json file setup for NOE
        
Density setup
        
            --density-json Json file setup for density fitting
        
Global setup in json format
        
            --setup-json This file duplicates all the other options. It also allows for
                         multistage minimization protocol, where each stage can have
                         different minimization terms
        
Miscellaneous
            
            --num-threads Number of OpenMP threads to use. 0 (default) means using maximum number of threads
            --score-only Score only without doing minimization
            --help Show this message and exit

