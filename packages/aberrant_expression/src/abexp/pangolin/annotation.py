import polars as pl

from abexp.utils.gff3 import read_gff3


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
      strand, gene_id and sites (the sorted 0-based positions of the splice sites, like start).
    """
    gff = read_gff3(gff3_file, ('gene_id', 'tag'))
    genes = gff.filter((pl.col('type') == 'gene') & pl.col('strand').is_in(['+', '-'])).select(
        'chrom', 'start', 'end', 'strand', 'gene_id',
    )
    tagged = pl.col('tag').str.split(',').list.eval(pl.element().is_in(list(transcript_tags))).list.any()
    sites = (
        gff
        .filter((pl.col('type') == 'exon') & tagged)
        .select(
            'chrom',
            'gene_id',
            # the first and last bases of the exon, 0-based
            pl.concat_list('start', pl.col('end') - 1).alias('sites'),
        )
        .explode('sites', empty_as_null=False)
        .unique()
        .group_by('chrom', 'gene_id')
        .agg(pl.col('sites').sort())
    )
    return genes.join(sites, on=['chrom', 'gene_id'], how='left', maintain_order='left').with_columns(
        pl.col('sites').fill_null([])
    )
