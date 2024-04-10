# Annotation and MAPS-based analysis of splicing prediction scores

1. Data
   1. Import the scores into Hail
   2. Add annotations for MAPS
   3. Group variants and calculate the number of singletons

2. Analysis
   1. Download grouped variants
   2. Calculate MAPS for each score-group combination
   3. Visualise the results

## The "Data" pipeline

`data/` contains [Hail](https://hail.is/) and Snakemake code that requires execution in Google Cloud and saves files to a Google Storage (GS) bucket. Copies of the generated files are available in `files/`.

### How to run

1. Create a new cluster: `hailctl dataproc start <cluster_name> --packages snakemake==6.6.1 --requester-pays-allow-buckets gnomad-public-requester-pays --project <project_name> --bucket <bucket_name> --region <region> --num-workers <N> --image-version=2.0.27-debian10` (as of this writing, the command runs `gcloud dataproc clusters create` and installs the following packages: `--metadata=^|||^WHEEL=gs://hail-common/hailctl/dataproc/0.2.74/hail-0.2.74-py3-none-any.whl|||PKGS=aiohttp==3.7.4|aiohttp_session>=2.7,<2.8|asyncinit>=0.2.4,<0.3|bokeh>1.3,<2.0|boto3>=1.17,<2.0|botocore>=1.20,<2.0|decorator<5|Deprecated>=1.2.10,<1.3|dill>=0.3.1.1,<0.4|gcsfs==0.8.0|fsspec==0.9.0|humanize==1.0.0|hurry.filesize==0.9|janus>=0.6,<0.7|nest_asyncio|numpy<2|pandas>=1.1.0,<1.1.5|parsimonious<0.9|PyJWT|python-json-logger==0.1.11|requests==2.25.1|scipy>1.2,<1.7|tabulate==0.8.3|tqdm==4.42.1|google-cloud-storage==1.25.*`)
2. Connect to the cluster: `gcloud beta compute ssh <user_name>@<cluster_name>-m --zone "<zone>" --project "<project_name>"`
3. `git clone` this repository and navigate to `data/`
4. Run the pipeline: `snakemake --cores all --configfile config.yaml --config gcp_rootdir="<bucket_name>/some_directory/"`

## The "Analysis" pipeline

`analysis/` contains scripts that calculate and visualise MAPS scores using files created in `data/`.

### How to run

1. `git clone` the code for CAPS (https://github.com/VCCRI/CAPS) into the same root directory as `MAPS_for_splicing/`
2. Navigate to `MAPS_for_splicing/analysis/`
3. `snakemake --cores all --config gcp="False"` (faster: uses copies from `files/`) or `snakemake --cores all --config gcp="True" gcp_rootdir="<bucket_name>/some_directory/"` (slower: uses GS files)
