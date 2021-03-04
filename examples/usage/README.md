# Usage example

Run
```
../../build/nrgmin --setup-json setup.json
```

The output will be in files `out.pdb` and `out.json`, which contains energy terms. `setup.json` has three
minimization stages and `out.json` will record final energy values for each stage and each model. If you
want to print every step, for example you can run

```
../../build/nrgmin --setup-json setup.json --print-step
```

Or add `"print_step": true` to `"options"` field in `setup.json` and rerun

```
../../build/nrgmin --setup-json setup.json
```

### Setup.json explained

The field `options` can contain all the flags listed in `--help`. Flags added to `options` will act globally, 
meaning that every minimization stage added to `stages` will use the value from `options` unless 
specified differently.

```
"options": {
    "rec_pdb": "BACE_4_rec_2models.pdb",
    "rec_psf": "BACE_4_rec.psf",
    "rec_prm": "BACE_4_rec_parm.prm",
    "rec_rtf": "BACE_4_rec_pdbamino.rtf",
    "lig_json": "BACE_4_lig.json",
    "lig_pdb": "BACE_4_lig_2models_close.pdb"
}
```

The above setup can be reproduced by using all the same flags via the command line:

```
../../build/nrgmin --rec-pdb BACE_4_rec_2models.pdb --rec-psf BACE_4_rec.psf --rec-prm BACE_4_rec_parm.prm --rec-rtf BACE_4_rec_pdbamino.rtf --lig-json BACE_4_lig.json --lig-pdb BACE_4_lig_2models_close.pdb 
```

The field `stages` is optional and can be added if you want to define 
separate minimization stages with different parameters. For example, if you want to run 
500 steps with the fixed ligand, and then 500 steps without torsions you can use the following:

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
            "nsteps": 500,
            "fix_ligand": true
        },
        {
            "nsteps": 500,
            "dihedrals": false
        }
    ]
}
```

This will be equivalent to specifying `nsteps` globally:

```
{
    "options": {
        "rec_pdb": "BACE_4_rec_2models.pdb",
        "rec_psf": "BACE_4_rec.psf",
        "rec_prm": "BACE_4_rec_parm.prm",
        "rec_rtf": "BACE_4_rec_pdbamino.rtf",
        "lig_json": "BACE_4_lig.json",
        "lig_pdb": "BACE_4_lig_2models_close.pdb",
        "nsteps": 500
    },
    "stages": [
        {
            "fix_ligand": true
        },
        {
            "dihedrals": false
        }
    ]
}
```