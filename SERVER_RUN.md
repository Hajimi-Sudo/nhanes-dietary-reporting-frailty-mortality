# Server Run Notes

Run the formal analysis on a clean Linux environment with the required R and
Python versions. Do not commit downloaded participant-level files, processed
cohorts, or generated logs to the repository.

## Environment

- Ubuntu 22.04 or compatible Linux distribution
- R 4.1.2 or newer
- Python 3.10 or newer
- R packages listed in `requirements-r.txt`
- Python packages listed in `requirements-python.txt`

## Commands

From the repository root:

```bash
python -m pip install -r requirements-python.txt
bash run_pipeline.sh
```

Install the R packages before running the pipeline, for example:

```bash
Rscript -e 'install.packages(c("survey", "survival", "haven", "readr", "dplyr"), repos="https://cloud.r-project.org")'
```

The download stage writes
`logs/nhanes_download_manifest.csv`, including source URLs and SHA-256 hashes.
Keep that manifest with the run record, but review it before release because
it contains the download date and file metadata.

## Run record

Record the repository commit, R version, Python version, package versions, and
the date of execution alongside the generated results. Results are intentionally
ignored by Git; release only the aggregate tables or figures that are suitable
for public sharing.
