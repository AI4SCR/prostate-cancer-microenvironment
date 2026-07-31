# Agent instructions

## Never run long/heavy jobs on the login node

This repo is worked on from an HPC login node (`curnagl`), not a compute
node. The login node is shared by everyone on the cluster -- running
CPU/GPU-intensive or long-running work directly in the interactive shell
starves other users and can get the job killed or flagged by admins.

Before running any command expected to take more than a couple of minutes,
or that's CPU/GPU/memory-intensive (UMAP fits, k-means over millions of
cells, R survival/heatmap scripts over the full cohort, anything touching
the full `intensity.parquet`/`intensity_normalized.parquet` tables), submit
it through SLURM instead of running it directly:

```bash
sbatch --wrap="uv run python scripts/figures/figureX_something.py" \
  --job-name=figureX --time=04:00:00 --cpus-per-task=8 --mem=32G \
  --output=/path/to/logfile.log
```

Adjust `--time`/`--cpus-per-task`/`--mem` to the job. Use `squeue -u $USER`
to check status and tail the `--output` log for progress instead of running
the command inline and blocking on it.

Quick, cheap commands (reading a small parquet file, checking file
existence, git operations, editing code) are fine to run directly -- this
rule is about anything that would otherwise sit there consuming login-node
CPU/memory for real time.
