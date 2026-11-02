# fedcat collinearity

Check how collinear candidate variants are with a selected (conditioned-on) set, from a PLINK2 pgen/pvar file set. Single-variant conditioning runs this in every cohort stratum to drop candidates too collinear with the selected set (`cond_check_coll`, `cond_rsq_threshold`).

```bash
fedcat collinearity --pfile chr1_genotypes --selected selected_ids.txt \
    --candidate candidate_ids.txt --out collinearity.tsv --threads 8
```

## Output

A TSV with a row per candidate:

| Column | Meaning |
|---|---|
| `ID` | Candidate variant ID |
| `cand_rsq` | R² of the candidate regressed on the selected set (its VIF R²) |
| `max_rsq` | Largest VIF R² in the selected set plus the candidate, each regressed on the others; also catches a candidate making a selected variant collinear |

A candidate whose pairwise r² with a selected variant exceeds `--rsq-threshold`, or that makes the set singular, gets `2.0` in both. Variants not in the pvar are skipped. The input files list one variant ID per line.

## Options

--8<-- "options/collinearity.md"
