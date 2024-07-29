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

#filter out X and Y
filtered_ht = ht.filter((ht.locus.contig != 'X') & (ht.locus.contig != 'Y'))

column_of_interest = ht.freq[0].AN
percentiles = ht.aggregate(hl.agg.approx_quantiles(column_of_interest, [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]))
p5, p10, p20, p30, p40, p50, p60, p70, p80, p90 = percentiles
total_samplesX1point6 = 0.8 * 2 * (xx_total + xy_total)

print(f"Dataset: {exome_or_genome}")
print(f"Total samples X 0.8 X 2: {total_samplesX1point6}")
print(f"5th Percentile: {p5}")
print(f"10th Percentile: {p10}")
print(f"20th Percentile: {p20}")
print(f"30th Percentile: {p30}")
print(f"40th Percentile: {p40}")
print(f"50th Percentile: {p50}")
print(f"60th Percentile: {p60}")
print(f"70th Percentile: {p70}")
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
