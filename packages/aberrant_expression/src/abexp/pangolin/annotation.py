import polars as pl

# the 9 columns of a GFF3 file
GFF3_SCHEMA = {
    'seqid': pl.String,
    'source': pl.String,
    'type': pl.String,
    'start': pl.Int64,
    'end': pl.Int64,
    'score': pl.String,
    'strand': pl.String,
    'phase': pl.String,
    'attributes': pl.String,
}


def _attribute(name):
    """The value of the GFF3 attribute `name`, null if the feature has none."""
    return pl.col('attributes').str.extract(f'(?:^|;){name}=([^;]*)')


def read_gff3_genes(gff3_file, transcript_tags):
    """The genes of a GENCODE GFF3 file and their annotated splice sites, as Pangolin uses them.

    The selection is that of AbExp's Pangolin annotation database: all genes on the + or - strand, and the exons with
    one of `transcript_tags`. The splice sites of a gene are the first and last bases of these exons. An exon belongs
    to the gene with its `gene_id` on its chromosome. A gene without such exons has no sites.

    Args:
      gff3_file: GENCODE GFF3 file, may be gzipped
      transcript_tags: the exons with one of these tags mark the splice sites, e.g. ['Ensembl_canonical']

    Returns:
      polars DataFrame with one row per gene, in the order of the file, and the columns chrom, start (0-based), end,
      strand, gene_id and sites (the sorted 1-based positions of the splice sites).
    """
    gff = pl.read_csv(gff3_file, separator='\t', has_header=False, comment_prefix='#', quote_char=None,
                      schema=GFF3_SCHEMA)
    genes = gff.filter((pl.col('type') == 'gene') & pl.col('strand').is_in(['+', '-'])).select(
        pl.col('seqid').alias('chrom'),
        (pl.col('start') - 1).alias('start'),
        'end',
        'strand',
        _attribute('gene_id').alias('gene_id'),
    )
    tagged = _attribute('tag').str.split(',').list.eval(pl.element().is_in(list(transcript_tags))).list.any()
    sites = (
        gff
        .filter((pl.col('type') == 'exon') & tagged)
        .select(
            pl.col('seqid').alias('chrom'),
            _attribute('gene_id').alias('gene_id'),
            pl.concat_list('start', 'end').alias('sites'),
        )
        .explode('sites', empty_as_null=False)
        .unique()
        .group_by('chrom', 'gene_id')
        .agg(pl.col('sites').sort())
    )
    return genes.join(sites, on=['chrom', 'gene_id'], how='left', maintain_order='left').with_columns(
        pl.col('sites').fill_null([])
    )
