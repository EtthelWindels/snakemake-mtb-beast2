For first BEAST runs or overwriting existing runs:
- set `action: "overwrite"` in `analyses.yaml`

For resuming BEAST runs:
- set `action: "resume"` in `analyses.yaml`
- if slurm job is finished before number of iterations is reached: re-run snakemake
- if number of iterations is reached before convergence:
	- move .log and .trees files back to running folder
	- remove combined chains and trees
	- re-run snakemake
