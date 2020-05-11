# NRGmin 
Utility for molecular mechanics energy minimization

## Installation ###

1. Install [libmol2](https://bitbucket.org/bu-structure/libmol2/src/master)

2. Build minimization executable 

        mkdir build && cd build
        cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON ..
        make && make test
       
3. The binary is in `build/nrgmin`

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
./nrgmin ---rec-pdb rec.pdb --rec-psf rec.psf --rec-prm rec.prm --rec-rtf rec.rtf --lig-json lig.json
```

If you want to use parameters from json file and coordinates from pdb

```
./nrgmin ---rec-pdb rec.pdb --rec-json rec.json --lig-pdb lig.pdb --lig-json lig.json
```

Note: PDB file can contain multiple models (the number of models must match in rec/lig mode)

### Energy function

Energy function setup is flexible, so you can switch terms on and off

```
./nrgmin --json mol.json --vdw-off --gbsa-on
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
      "gbsa": true,
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
          "weight": 10,
          "length": 10,
          "error": 1,
          "atom1": 3598,
          "atom2": 3630
        }
      ]
    }
  ]
}
```

### Examples

More setup examples are [here](./examples/energy_setup) and usage examples [here](./examples/usage).

### All flags

#### Input/Output

Output control:
        
            --out-pdb Where to write the minimized molecule(s)
            --out-json Log file with energy terms
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
            --gbsa-on/gbsa-off GBSA
        
Fixed atoms (in rec/lig mode ligand atom IDs must be increased by the number of receptor atoms)
        
            --fixed-pdb PDB file with fixed atoms. Only atom ID is read (lines starting with ATOM)
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
        
            --score-only Score only without doing minimization
            --help Show this message and exit

