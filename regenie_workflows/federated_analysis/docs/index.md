# ![FedCAT](fedcat_logo.png){ .fedcat-hero }

**Federated Conditional Association Testing** across All of Us (AoU) and UK Biobank (UKB).

`fedcat` runs a cross-biobank whole-genome sequencing association analysis from one command and one config file. It runs REGENIE in each biobank's own cloud, meta-analyses the results across biobanks and strata, then finds the independent signals by federated forward-selection conditioning. Individual-level data stays in each biobank; only summary statistics move.

```bash
fedcat discover --config my_config.json
```

## What a run does

| Phase | What it does |
|---|---|
| `step1` | REGENIE step 1 (whole-genome regression), per cohort and stratum |
| `step2_singles` | REGENIE step 2 single-variant tests, per cohort and stratum |
| `step2_aggregates` | REGENIE step 2 aggregate (burden) tests, per cohort and stratum |
| `meta_analysis` | Copies UKB summary statistics to the AoU bucket and meta-analyses everything, with ACAT tests of aggregates |
| `step0_subset` | Collects the significant meta-analysis results as the starting point for conditioning |
| `conditioning_singles` / `conditioning_aggregates` | Federated forward selection: test conditionally in every cohort, meta-analyse, select the top signal per chromosome, repeat |

See [How a run works](guide/how-it-works.md) for details.

## Key features

- **Up to 20 phenotypes at once**, using REGENIE's multi-phenotype support in the discovery steps (more with [phenotype grouping](guide/configuration.md#large-phenotype-sets)).
- **Resumable.** Progress is saved to state files, so a run picks up where it stopped after a transient failure or a config fix.
- **Runs from an AoU Jupyter terminal**, which also submits the UKB jobs to DNAnexus.
- **Standalone tools** for meta-analysis, ACAT, annotation conversion and collinearity checks (see [Commands](commands/index.md)).

## Where next

<div class="grid cards" markdown>

- **[Installation](getting-started/installation.md)**

    Install FedCAT on an AoU Jupyter instance and stage the workflows.

- **[Quick start](getting-started/quick-start.md)**

    Write a minimal config, dry-run it, and run your first analysis.

- **[Configuration](guide/configuration.md)**

    A guided tour of the config file, section by section.

- **[Commands](commands/index.md)**

    Every subcommand, with examples and all of its options.

</div>

## Contact

Questions and suggestions: Harry Wright (<h.wright4@exeter.ac.uk>). Logo by Theo Lupu (<t.lupu@exeter.ac.uk>).
