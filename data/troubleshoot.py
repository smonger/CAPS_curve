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

#gnomad = preprocessing(ht_path, context_ht_path, coverage_ht_path, {"female": xx_total, "male": xy_total})

from misc import (
    get_an_adj_criteria,
    prepare_ht,
    filter_vep_to_canonical_transcripts,
    get_worst_consequence_with_non_coding,
)


data_ht = ht_path
context_ht = context_ht_path
coverage_ht = coverage_ht_path
sex_split = {"female": xx_total, "male": xy_total}

context = hl.read_table(context_ht)
coverage_ht = hl.read_table(coverage_ht)
ht = hl.read_table(data_ht)
ht = ht.annotate(coverage=coverage_ht[ht.locus].median_approx)



# Define the column of interest
column_of_interest = ht.freq[0].AN

# Aggregate to get the 25th and 75th percentiles
percentiles = ht.aggregate(hl.agg.approx_quantiles(column_of_interest, [0.25, 0.5, 0.80, 0.90]))

# Compute the mode
mode = ht.aggregate(hl.agg.mode(column_of_interest))

# Extract percentiles
p10, p25, p75, p90 = percentiles

print(f"Mode: {mode}")
print(f"25th Percentile: {p25}")
print(f"50th Percentile: {p50}")
print(f"80th Percentile: {p80}")
print(f"90th Percentile: {p90}")

    # Allele number (AN) adjustment. This retains only those variants,
    # in which the call was made in at least 80% (see "an_cutoff") of
    # potential carriers.
#    ht = ht.filter(get_an_adj_criteria(ht, sex_split))

    # Filter the table so that only those variants that have AF>0 and
    # filter PASS are retained. The first condition is necessary
    # because in gnomAD variants that were excluded from the analysis
    # through QC have AF=0.
#    ht = ht.filter(
#        (ht.freq[0].AF > 0)
#        & (ht.filters.length() == 0)
#        & (ht.vep.variant_class == "SNV")
#    )

##    ht = filter_vep_to_canonical_transcripts(ht)
 #   ht = get_worst_consequence_with_non_coding(ht)

#    context = context[ht.key]
    # The 2020 version of MAPS uses methylation.
    # Function "prepare_ht" annotates the input table with methylation level,
    # coverage (optional), CpG/Non-CpG info, context for mutability
    # (ref allele in the middle plus two bases to the left and to the right)
    # and other information.
#    ht = prepare_ht(
#        ht.annotate(context=context.context, methylation=context.methyl_mean),
#        trimer=True,
#        annotate_coverage=False,
#    )

    #ht = ht.annotate(
    #    mu=mutation_rates[
    #        hl.struct(
    #            context=ht.context,
    #            ref=ht.ref,
    #            alt=ht.alt,
    #            methylation_level=ht.methylation_level,
    #        )
    #    ].mu_snp
    #)

#    return ht
