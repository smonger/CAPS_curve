from pathlib import Path
import hail as hl
from preprocessing import preprocessing
from hail.expr.functions import int32, float64

#### Specify Exome or Genome here ####
exome_or_genome = "genome"

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
mutation_rates_ht_path = gnomad_gs_path + "papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht"
context_ht_path = gnomad_gs_path + "resources/context/grch38_context_vep_annotated.v105.ht"
coverage_exomes_ht_path = gnomad_gs_path + "release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.ht"
coverage_genomes_ht_path = gnomad_gs_path + "release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht/"
lof_metrics_by_gene = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

input_splicing_exomes = mikgud_gs_path + "hg38_gnomad_v4_exome_splicing_scores_27052024.bgz"
input_splicing_genomes = mikgud_gs_path + "hg38_gnomad_v4_genome_splicing_scores_27052024.bgz"
output_splicing_exomes = mikgud_gs_path + "hg38_gnomad_v4_exome_splicing_scores_27052024_data_output.bgz"
output_splicing_genomes = mikgud_gs_path + "hg38_gnomad_v4_genome_splicing_scores_27052024_data_output.bgz"
by_csq_file_genome = mikgud_gs_path + "gnomad_v4.1.0_genomes_by_csq.tsv"
by_csq_file_exome = mikgud_gs_path + "gnomad_v4.1.0_exomes_by_csq.tsv"

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

gnomad = preprocessing(ht_path, context_ht_path, mutation_rates_ht_path, coverage_ht_path, {"female": xx_total, "male": xy_total})

if create_by_csq_file:
  gnomad.groupby("context", "ref", "alt", "mu", "methylation_level", "worst_csq", "coverage").aggregate(
    variant_count = hl.agg.count(), singleton_count=hl.agg.count_where(gnomad.freq[0].AC == 1)).export(by_csq_file)

data = hl.import_table(input_file, delimiter="\t", missing="")
data = data.filter(data.start != "start")

data = data.transmute(
            ds_ag=float64(data.DS_AG),
            ds_al=float64(data.DS_AL),
            ds_dg=float64(data.DS_DG),
            ds_dl=float64(data.DS_DL),
            sai_sum=float64(data.sai_sum),
            sai_max=float64(data.sai_max),
            sai_loss125=float64(data.sai_loss125),
            sai_loss150=float64(data.sai_loss150),
            sai_1xgb=float64(data.sai_1xgb),
            ssm_1e=float64(data.ssm_1e),
            ssm_1amne=float64(data.ssm_1amne),
            ssm_1amdi=float64(data.ssm_1amdi),
            ssm_1amnsu=float64(data.ssm_1amnsu),
            sai_2xgb=float64(data.sai_2xgb),
            ssm_2e=float64(data.ssm_2e),
            ssm_2amne=float64(data.ssm_2amne),
            ssm_2amdi=float64(data.ssm_2amdi),
            ssm_2amnsu=float64(data.ssm_2amnsu),
            absplice=float64(data.absplice),
            pangolin=float64(data.pangolin),
            start=int32(data.start))

# Index the variants for which scores are available to match the keys of Hail tables
data = data.key_by(locus=hl.locus(data.chr, data.start, reference_genome="GRCh38"), alleles=[data.ref, data.alt])

data = data[gnomad.key]
gnomad = gnomad.annotate(            
            ds_ag=data.ds_ag,
            ds_al=data.ds_al,
            ds_dg=data.ds_dg,
            ds_dl=data.ds_dl,
            sai_sum=data.sai_sum,
            sai_max=data.sai_max,
            sai_loss125=data.sai_loss125,
            sai_loss150=data.sai_loss150,
            sai_1xgb=data.sai_1xgb,
            ssm_1e=data.ssm_1e,
            ssm_1amne=data.ssm_1amne,
            ssm_1amdi=data.ssm_1amdi,
            ssm_1amnsu=data.ssm_1amnsu,
            sai_2xgb=data.sai_2xgb,
            ssm_2e=data.ssm_2e,
            ssm_2amne=data.ssm_2amne,
            ssm_2amdi=data.ssm_2amdi,
            ssm_2amnsu=data.ssm_2amnsu,
            absplice=data.absplice,
            pangolin=data.pangolin)

gnomad = gnomad.annotate(AC=gnomad.freq[0].AC)
gnomad = gnomad.filter(gnomad.worst_csq!="stop_gained")
gnomad.select(
            "context",
            "ref",
            "alt",
            "methylation_level",
            "mu",
            "worst_csq",
            "protein_coding",
            "coverage",
            "AC",
            "ds_ag",
            "ds_al",
            "ds_dg",
            "ds_dl",
            "sai_sum",
            "sai_max",
            "sai_loss125",
            "sai_loss150",
            "ssm_1xgb",
            "ssm_1e",
            "ssm_1amne",
            "ssm_1amdi",
            "ssm_1amnsu",
            "ssm_2xgb",
            "ssm_2e",
            "ssm_2amne",
            "ssm_2amdi",
            "ssm_2amnsu",
            "absplice",
            "pangolin").export(output_file)
