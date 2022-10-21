# This is simply a rewrite of the entire snakemake pipeline "Data" --
# with less flexibility and more hard-coded variables

import hail as hl
from typing import Union


def trimer_from_heptamer(
    t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    trimer_expr = hl.cond(hl.len(t.context) == 7, t.context[2:5], t.context)
    return (
        t.annotate_rows(context=trimer_expr)
        if isinstance(t, hl.MatrixTable)
        else t.annotate(context=trimer_expr)
    )


def annotate_variant_types(
    t: Union[hl.MatrixTable, hl.Table], heptamers: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (
        ((t.ref == "A") & (t.alt == "G"))
        | ((t.ref == "G") & (t.alt == "A"))
        | ((t.ref == "T") & (t.alt == "C"))
        | ((t.ref == "C") & (t.alt == "T"))
    )
    cpg_expr = (
        (t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1 : mid_index] == "C")
    ) | (
        (t.ref == "C")
        & (t.alt == "T")
        & (t.context[mid_index + 1 : mid_index + 2] == "G")
    )
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (
        hl.case()
        .when(t.cpg, "CpG")
        .when(t.transition, "non-CpG transition")
        .default("transversion")
    )
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )
    else:
        return t.annotate(
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )


def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return (
        hl.switch(base)
        .when("A", "T")
        .when("T", "A")
        .when("G", "C")
        .when("C", "G")
        .default(base)
    )


def reverse_complement_bases(
    bases: hl.expr.StringExpression,
) -> hl.expr.StringExpression:
    return hl.delimit(
        hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), ""
    )


def collapse_strand(
    ht: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    collapse_expr = {
        "ref": hl.cond(
            ((ht.ref == "G") | (ht.ref == "T")),
            reverse_complement_bases(ht.ref),
            ht.ref,
        ),
        "alt": hl.cond(
            ((ht.ref == "G") | (ht.ref == "T")),
            reverse_complement_bases(ht.alt),
            ht.alt,
        ),
        "context": hl.cond(
            ((ht.ref == "G") | (ht.ref == "T")),
            reverse_complement_bases(ht.context),
            ht.context,
        ),
        "was_flipped": (ht.ref == "G") | (ht.ref == "T"),
    }
    return (
        ht.annotate(**collapse_expr)
        if isinstance(ht, hl.Table)
        else ht.annotate_rows(**collapse_expr)
    )


def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter(
            (hl.len(ht.ref) == 1)
            & (hl.len(ht.alt) == 1)
            & ht.context.matches(f"[ATCG]{{{str_len}}}")
        )
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows(
            (hl.len(ht.ref) == 1)
            & (hl.len(ht.alt) == 1)
            & ht.context.matches(f"[ATCG]{{{str_len}}}")
        )
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    if annotate_coverage:
        annotation["exome_coverage"] = ht.coverage.exomes.median
    return (
        ht.annotate(**annotation)
        if isinstance(ht, hl.Table)
        else ht.annotate_rows(**annotation)
    )


def get_an_adj_criteria(
    hail_table,
    sex_split,
    an_cutoff: float = 0.8,
):
    """Get lower bound allele number (AN) thresholds.

    hail_table -- variants in ht format
    sex_split -- how many males, how many females
    an_cutoff -- percent of genotype calls
    """
    return (
        hl.case()
        .when(
            hail_table.locus.in_autosome_or_par(),
            hail_table.freq[0].AN >= an_cutoff * 2 * sum(sex_split.values()),
        )
        .when(
            hail_table.locus.in_x_nonpar(),
            hail_table.freq[0].AN
            >= an_cutoff * (sex_split["male"] + sex_split["female"] * 2),
        )
        .when(
            hail_table.locus.in_y_nonpar(),
            hail_table.freq[0].AN >= an_cutoff * sex_split["male"],
        )
        .or_missing()
    )


def annotate_quartiles(ht, variable):
    quartiles = ht.aggregate(
        hl.agg.approx_quantiles(ht[variable], [0, 0.25, 0.5, 0.75, 1])
    )
    ht = ht.annotate(
        quartile=hl.case()
        .when(ht[variable] <= quartiles[1], 1)
        .when(
            (ht[variable] > quartiles[1]) & (ht[variable] <= quartiles[2]),
            2,
        )
        .when(
            (ht[variable] > quartiles[2]) & (ht[variable] <= quartiles[3]),
            3,
        )
        .when(ht[variable] > quartiles[3], 4)
        .or_missing()
    )
    return ht


splicesm = (
    hl.import_table(
        "gs://vccri-mikgud/prediction_scores_v1.txt",
        delimiter="\t",
        missing="",
        # TODO: this is ad-hoc, make this work for a more
        # general case of [1,N] scores
        types={
            "spliceAI_max": hl.tfloat64,
            "spliceAI_xgb": hl.tfloat64,
            "spliceSM": hl.tfloat64,
            "start": hl.tint32,
        },
    )
    # TODO: this is ad-hoc, make this work for a more general
    # case of [1,N] scores
    .select(
        "chr", "start", "ref", "alt", "spliceAI_max", "spliceAI_xgb", "spliceSM"
    ).key_by("chr", "start", "ref", "alt")
)

# Load all exonic variants from gnomAD v2 and perform QC
exomes = hl.read_table("gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht/")
exomes = exomes.filter(
    # Filter the table so that only those variants that have
    # AF>0 and filter PASS are retained. The first condition
    # is necessary because in gnomAD variants that were
    # excluded from the analysis through QC have AF=0.
    (exomes.freq[0].AF > 0)
    & (exomes.filters.length() == 0)
    # SpliceSM can be calculated for SNVs as well as indels,
    # but MAPS only works on SNVs
    & (exomes.vep.variant_class == "SNV")
)
# Allele number (AN) adjustment
# TODO: make this a CONST variable
exomes = exomes.filter(get_an_adj_criteria(exomes, {"female": 57787, "male": 67961}))

# Annotate the variants with SpliceSM scores
# TODO: this is ad-hoc, make this work for a more general case
# of [1,N] scores
exomes = exomes.annotate(
    spliceAI_max=splicesm[
        hl.struct(
            chr=exomes.locus.contig,
            start=exomes.locus.position,
            ref=exomes.alleles[0],
            alt=exomes.alleles[1],
        )
    ].spliceAI_max,
    spliceAI_xgb=splicesm[
        hl.struct(
            chr=exomes.locus.contig,
            start=exomes.locus.position,
            ref=exomes.alleles[0],
            alt=exomes.alleles[1],
        )
    ].spliceAI_xgb,
    spliceSM=splicesm[
        hl.struct(
            chr=exomes.locus.contig,
            start=exomes.locus.position,
            ref=exomes.alleles[0],
            alt=exomes.alleles[1],
        )
    ].spliceSM,
)
exomes = annotate_quartiles(exomes, "spliceAI_max").rename(
    {"quartile": "spliceAI_max_quartile"}
)
exomes = annotate_quartiles(exomes, "spliceAI_xgb").rename(
    {"quartile": "spliceAI_xgb_quartile"}
)
exomes = annotate_quartiles(exomes, "spliceSM").rename(
    {"quartile": "spliceSM_quartile"}
)

# Methylation, mutability and other context data
context = hl.read_table("gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht")
context = context[exomes.key]
exomes = prepare_ht(
    exomes.annotate(context=context.context, methylation=context.methylation),
    trimer=True,
    annotate_coverage=False,
)
mutation_rates = hl.read_table("gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht")
exomes = exomes.annotate(
    mu=mutation_rates[
        hl.struct(
            context=exomes.context,
            ref=exomes.ref,
            alt=exomes.alt,
            methylation_level=exomes.methylation_level,
        )
    ].mu_snp
)

exomes.group_by(
    "context",
    "ref",
    "alt",
    "methylation_level",
    "mu",
    "spliceSM_quartile",
    "spliceAI_xgb_quartile",
    "spliceAI_max_quartile",
).aggregate(
    variant_count=hl.agg.count(),
    singleton_count=hl.agg.count_where(exomes.freq[0].AC == 1),
).export(
    "gs://vccri-mikgud/SpliceSM_quartiles.tsv"
)
