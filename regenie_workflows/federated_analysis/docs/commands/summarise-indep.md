# fedcat summarise-indep

Rewrite a finished run's [conditioning results](../guide/outputs.md#conditioning) from its per-step `summary_output` tables and the selected IDs in a step state file: to regenerate lost or changed results, or summarise runs from old fedcat, fedcond and looper versions. Nothing is uploaded and no state is changed.

```bash
fedcat summarise-indep --run-dir fedcat_runs/run1_BMI_singlevar_conditioning_run1 \
    --cond-type singles --pheno BMI
```

The last step, the state file holding `passing_ids` and the statistic column (`LOG10P`, or `P` in older runs) are detected from the run directory. Results go to `<run-dir>/<type>_cond_results`, unless it exists; choose another with `--output-dir`.

## Options

--8<-- "options/summarise-indep.md"
