
import hail as hl
from preprocessing import preprocessing
from hail.expr.functions import int32, float64

#### Specify Exome or Genome here ####
exome_or_genome = "exome"

#### Specify whether t generate by_csq files. Only required when refitting CAPS
create_by_csq_file = False

#### Constants and file paths ####
xy_exome_total = 363624
xx_exome_total = 367323
xy_genome_total = 37237
xx_genome_total = 38942

gnomad_gs_path = "gs://gcp-public-data--gnomad/"
mikgud_gs_path = "gs://vccri-mikgud-uscentral1/"

exomes_ht_path = gnomad_gs_path + "release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht"
genomes_ht_path = gnomad_gs_path + "release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/"
context_ht_path = gnomad_gs_path + "resources/context/grch38_context_vep_annotated.v105.ht"
coverage_exomes_ht_path = gnomad_gs_path + "release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.ht/"
coverage_genomes_ht_path = gnomad_gs_path + "release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht/"
lof_metrics_by_gene = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

input_splicing_exomes = mikgud_gs_path + "hg38_gnomad_v4_exome_splicing_scores_27052024.bgz"
input_splicing_genomes = mikgud_gs_path + "hg38_gnomad_v4_genome_splicing_scores_27052024.bgz"
output_splicing_exomes = mikgud_gs_path + "hg38_gnomad_v4_exome_splicing_scores_27052024_data_output.bgz"
output_splicing_genomes = mikgud_gs_path + "hg38_gnomad_v4_genome_splicing_scores_27052024_data_output.bgz"
by_csq_file_genome = mikgud_gs_path + "gnomad_v4.1.0_genomes_by_csq.tsv"
by_csq_file_exome = mikgud_gs_path + "gnomad_v4.1.0_exomes_by_csq.tsv"

f1 = mikgud_gs_path + "f1.ht"
f2 = mikgud_gs_path + "f2.ht"
f3 = mikgud_gs_path + "f3.ht"
f4 = mikgud_gs_path + "f4.ht"

if exome_or_genome == "genome":
  ht_path = genomes_ht_path
  coverage_ht_path = coverage_genomes_ht_path
  input_file = input_splicing_genomes
  output_file = output_splicing_genomes
  xy_total = xy_genome_total
  xx_total = xx_genome_total
  by_csq_file = by_csq_file_genome
elif exome_or_genome == "exome":
  ht_path = exomes_ht_path
  coverage_ht_path = coverage_exomes_ht_path
  input_file = input_splicing_exomes
  output_file = output_splicing_exomes
  xy_total = xy_exome_total
  xx_total = xx_exome_total
  by_csq_file = by_csq_file_exome

#### RUN ####

t1 = hl.read_table(f1)
num_rows = t1.count()
num_cols = len(t1.row)
print(f'The table has {num_rows} rows and {num_cols} columns.')
t2 = hl.read_table(f2)
num_rows = t2.count()
num_cols = len(t2.row)
print(f'The table has {num_rows} rows and {num_cols} columns.')
t3 = hl.read_table(f3)
num_rows = t3.count()
num_cols = len(t3.row)
print(f'The table has {num_rows} rows and {num_cols} columns.')
t4 = hl.read_table(f4)
num_rows = t4.count()
num_cols = len(t4.row)
print(f'The table has {num_rows} rows and {num_cols} columns.')
