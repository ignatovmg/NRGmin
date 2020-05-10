# NRGmin 
Utility for molecular mechanics energy minimization

### Installation ###

1. Install [libmol2](https://bitbucket.org/bu-structure/libmol2/src/master)

2. Build minimization executable 

        mkdir build && cd build
        cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS ..
        make && make test
       
3. The binary is in `build/nrgmin`

### How to run

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
./nrgmin ---rec-pdb mol.pdb --rec-psf mol.psf --rec-prm mol.prm --rec-rtf mol.rtf --lig-json lig.json
```

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


### Master json format

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