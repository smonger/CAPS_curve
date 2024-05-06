# Annotation and MAPS-based analysis of splicing prediction scores

1. Data
   1. Import the scores into Hail
   2. Add annotations for MAPS
   3. Group variants and calculate the number of singletons

2. Analysis
   1. Download grouped variants
   2. Calculate MAPS for each score-group combination
   3. Visualise the results

Tested with Hail version 0.2.107 and Snakemake 7.32

## The "Data" pipeline

`data/` contains [Hail](https://hail.is/) and Snakemake code that requires execution in Google Cloud and saves files to a Google Storage (GS) bucket.

### How to run

1. Create a new cluster: `hailctl dataproc start <cluster_name> --packages snakemake --requester-pays-allow-buckets gnomad-public-requester-pays --project <project_name> --bucket <bucket_name> --region <region> --num-workers <N> --image-version=2.0.27-debian10`
2. Connect to the cluster: `gcloud beta compute ssh <user_name>@<cluster_name>-m --project "<project_name>"`
3. `git clone` this repository and navigate to `data/`
4. Run the pipeline: `snakemake --cores all --configfile config.yaml --config gcp_rootdir="<bucket_name>/some_directory/"`

Alternatively, in Step 4 you can submit the pipeline as a job. Create `job.py` containing the following:
```
import snakemake
snakemake.main(
	[
		"--snakefile",
		"/path/to/Snakefile",
		"--cores",
		"all",
		"--configfile",
		"/path/to/config.yaml",
		"--config",
		'gcp_rootdir="<bucket_name>/some_directory/"',
	]
)
```
Submit the script with `hailctl dataproc submit <cluster_name> job.py`

## The "Analysis" pipeline

`analysis/` contains scripts that calculate and visualise MAPS scores using files created in `data/`.

### How to run

1. `git clone` the code for CAPS (https://github.com/VCCRI/CAPS) into the same root directory as `MAPS_for_splicing/`
2. Navigate to `MAPS_for_splicing/analysis/`
3. `snakemake --cores all --config gcp="True" gcp_rootdir="<bucket_name>/some_directory/"` (slower; will download the GS files each time) or `snakemake --cores all --config gcp="False"` (faster; will re-use the contents of `files/`; create the directory and download the required files from GS into it if necessary)
