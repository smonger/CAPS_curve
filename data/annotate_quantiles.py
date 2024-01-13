import hail as hl
def annotate_quantiles(ht, variable):
    """Group into bins by exponential quantiles.

    ht -- variants (Hail table)
    variable -- variable of interest
    """

    quantiles = ht.aggregate(
        # 1 - 1 / 2 ^ seq(0, 9)
        hl.agg.approx_quantiles(
            ht[variable],
            [
                0.0000,
                0.5000,
                0.7500,
                0.8750,
                0.9375,
                0.9688,
                0.9844,
                0.9922,
                0.9961,
                0.9980,
            ],
        )
    )
    ht = ht.annotate(
        variable_bin=hl.case()
        .when(ht[variable] >= quantiles[9], 9)
        .when(ht[variable] >= quantiles[8], 8)
        .when(ht[variable] >= quantiles[7], 7)
        .when(ht[variable] >= quantiles[6], 6)
        .when(ht[variable] >= quantiles[5], 5)
        .when(ht[variable] >= quantiles[4], 4)
        .when(ht[variable] >= quantiles[3], 3)
        .when(ht[variable] >= quantiles[2], 2)
        .when(ht[variable] >= quantiles[1], 1)
        .when(ht[variable] >= quantiles[0], 0)
        .or_missing()
    )
    return ht
