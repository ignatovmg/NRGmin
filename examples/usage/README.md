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