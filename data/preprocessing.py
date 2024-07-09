import hail as hl

from misc import (
    get_an_adj_criteria,
    prepare_ht,
    filter_vep_to_canonical_transcripts,
    get_worst_consequence_with_non_coding,
)


def preprocessing(
    data_ht,
    context_ht,
    mutation_rates_ht,
    coverage_ht,
    sex_split,
):
    """Preprocessing steps (QC and annotations).

    data_ht -- WES or WGS variants (Hail table)
    context_ht -- context (Hail table)
    mutation_rates_ht -- mutability/mutation rates (Hail table)
    coverage_ht -- coverage (Hail table)
    sex_split -- how many males, how many females
    """

    context = hl.read_table(context_ht)
    mutation_rates = hl.read_table(mutation_rates_ht)
    coverage_ht = hl.read_table(coverage_ht)

    ht = hl.read_table(data_ht)

    ht = ht.annotate(coverage=coverage_ht[ht.locus].median_approx)

    # Allele number (AN) adjustment. This retains only those variants,
    # in which the call was made in at least 80% (see "an_cutoff") of
    # potential carriers.
    ht = ht.filter(get_an_adj_criteria(ht, sex_split))

    # Filter the table so that only those variants that have AF>0 and
    # filter PASS are retained. The first condition is necessary
    # because in gnomAD variants that were excluded from the analysis
    # through QC have AF=0.
    ht = ht.filter(
        (ht.freq[0].AF > 0)
        & (ht.filters.length() == 0)
        & (ht.vep.variant_class == "SNV")
    )

    ht = filter_vep_to_canonical_transcripts(ht)
    ht = get_worst_consequence_with_non_coding(ht)

    context = context[ht.key]
    # The 2020 version of MAPS uses methylation.
    # Function "prepare_ht" annotates the input table with methylation level,
    # coverage (optional), CpG/Non-CpG info, context for mutability
    # (ref allele in the middle plus two bases to the left and to the right)
    # and other information.
    ht = prepare_ht(
        ht.annotate(context=context.context, methylation=context.methyl_mean),
        trimer=True,
        annotate_coverage=False,
    )

    ht = ht.annotate(
        mu=mutation_rates[
            hl.struct(
                context=ht.context,
                ref=ht.ref,
                alt=ht.alt,
                methylation_level=ht.methylation_level,
            )
        ].mu_snp
    )

    return ht
