import numpy as np
import polars as pl
import polars.testing
import pytest
import torch
import torch.nn.functional as F

import synthetic as s
from abexp.pangolin import Pangolin, PangolinModels, PangolinNet, read_gff3_genes

# The expected rows come from upstream Pangolin with the test models, see make_expected.py. PyTorch sums in another
# order for the unpadded convolutions and the batches, which changes the scores in the last float32 digits.
ABS_TOL = 1e-6


def make_net(model_idx, n_channels=s.N_CHANNELS):
    """A network with the kernel sizes and dilations of Pangolin and the fixed random weights of synthetic.py."""
    net = PangolinNet(n_channels=n_channels)
    shapes = {name: tuple(value.shape) for name, value in net.state_dict().items()}
    net.load_state_dict({name: torch.from_numpy(value) for name, value in
                         s.random_state_dict(shapes, model_idx).items()})
    return net.eval()


@pytest.fixture(scope='module')
def models():
    # the order of make_expected.py: 3 replicates for each head
    return PangolinModels([[make_net(3 * head + replicate) for replicate in range(3)] for head in range(4)])


@pytest.fixture(scope='module')
def genome(tmp_path_factory):
    path = tmp_path_factory.mktemp('genome')
    s.write_fasta(path / 'genome.fa')
    s.write_gff3(path / 'genes.gff3')
    return path


@pytest.fixture
def predict(genome, models, tmp_path):
    def predict(*records, mask=True):
        s.write_vcf(tmp_path / 'variants.vcf', records)
        genes = read_gff3_genes(genome / 'genes.gff3', s.TRANSCRIPT_TAGS)
        pangolin = Pangolin(str(genome / 'genome.fa'), genes, models, distance=s.DISTANCE, mask=mask)
        return pangolin.predict_df(str(tmp_path / 'variants.vcf'))
    return predict


def assert_rows(df, *rows):
    expected = pl.DataFrame(rows, schema=Pangolin.SCHEMA, orient='row')
    pl.testing.assert_frame_equal(df, expected, check_exact=False, rel_tol=0, abs_tol=ABS_TOL)


def test_snv_plus_strand(predict):
    assert_rows(
        predict(s.SNV_PLUS),
        ('chr1:6000:T>G', 'PLUS1', 0.03144197538495064, 0, -0.007421374320983887, 20, []),
    )


def test_snv_minus_strand(predict):
    assert_rows(
        predict(s.SNV_MINUS),
        ('chr1:7150:A>C', 'MINUS1', 0.03211745619773865, 0, -0.004649688955396414, -50, []),
    )


def test_insertion(predict):
    # the outputs of the first base and the 4 inserted bases collapse into their maximum
    assert_rows(
        predict(s.INSERTION),
        ('chr1:5950:C>CTTAG', 'PLUS1', 0.12284623831510544, 0, -0.018105357885360718, -50, []),
    )


def test_deletion(predict):
    # the 3 deleted bases have no alt output and count as 0
    assert_rows(
        predict(s.DELETION),
        ('chr1:5980:GCTA>G', 'PLUS1', 0.04620047410329183, 23, -0.0019849836826324463, 40, []),
    )


def test_multiallelic_record(predict):
    # one row per ALT allele
    assert_rows(
        predict(s.MULTIALLELIC),
        ('chr1:5960:T>A', 'PLUS1', 0.012907475233078003, 20, 0.0, -50, []),
        ('chr1:5960:T>C', 'PLUS1', 0.025438398122787476, 0, 0.0, -50, []),
    )


def test_overlapping_genes_on_one_strand(predict):
    # one prediction, masked with the splice sites of each gene
    assert_rows(
        predict(s.OVERLAPPING_GENES),
        ('chr1:6100:T>C', 'PLUS1', 0.03092496655881405, 0, -0.006370723247528076, 10, []),
        ('chr1:6100:T>C', 'PLUS2', 0.03092496655881405, 0, -0.000657667696941644, 50, []),
    )


def test_genes_on_both_strands(predict):
    assert_rows(
        predict(s.BOTH_STRANDS),
        ('chr1:6580:G>A', 'PLUS1', 0.012087802402675152, -1, -0.002651393413543701, 20, []),
        ('chr1:6580:G>A', 'PLUS2', 0.012087802402675152, -1, 0.0, -50, []),
        ('chr1:6580:G>A', 'MINUS1', 0.06992657482624054, 0, -0.003699402092024684, 40, []),
    )


def test_gene_without_splice_sites(predict):
    # without splice sites, the masking sets every loss to 0
    assert_rows(
        predict(s.NO_SITES),
        ('chr1:7400:G>T', 'MINUS2', 0.02329273521900177, 0, 0.0, -50, ['NoAnnotatedSitesToMaskForThisGene']),
    )


def test_masking_off(predict):
    # the same variant as in test_snv_plus_strand; the loss at an unannotated site counts
    assert_rows(
        predict(s.SNV_PLUS, mask=False),
        ('chr1:6000:T>G', 'PLUS1', 0.03144197538495064, 0, -0.044891681522130966, 0, []),
    )


def test_n_in_window(predict):
    # N is all zeros in the one-hot encoding
    assert_rows(
        predict(s.N_IN_WINDOW),
        ('chr2:6020:T>C', 'PLUS3', 0.029540440067648888, 0, -0.004769896622747183, 20, []),
    )


def test_chromosome_start(predict):
    # the rows of chr3n:6120:A>G, where 6000 N precede the same sequence
    assert_rows(
        predict(s.CHROM_START),
        ('chr3s:120:A>G', 'PLUS4', 0.04733463004231453, 0, -0.0016608238220214844, 30, []),
    )


def test_chromosome_end(predict):
    # the rows of chr4n:5900:A>G, where 6000 N follow the same sequence
    assert_rows(
        predict(s.CHROM_END),
        ('chr4s:5900:A>G', 'MINUS3', 0.03115805983543396, 0, -0.001647988916374743, -50, []),
    )


def test_variant_at_gene_start(predict):
    # upstream Pangolin misses PLUS1; the row is that of upstream Pangolin with PLUS1 starting at 5797
    assert_rows(
        predict(s.GENE_START),
        ('chr1:5800:A>C', 'PLUS1', 0.018580934032797813, -5, -0.018892541527748108, 0, []),
    )


def test_deletion_into_gene(predict):
    # the ref allele overlaps PLUS1; upstream Pangolin only checks the position of the variant, see
    # test_variant_at_gene_start
    assert_rows(
        predict(s.DELETION_INTO_GENE),
        ('chr1:5798:ACAT>A', 'PLUS1', 0.05373251438140869, 0, -0.5394827524820963, 2, []),
    )


def test_vcf_chromosome_without_chr_prefix(predict):
    # the variant of test_snv_plus_strand; the FASTA file names the chromosome chr1
    assert_rows(
        predict(('1', 6000, 'T', 'G')),
        ('1:6000:T>G', 'PLUS1', 0.03144197538495064, 0, -0.007421374320983887, 20, []),
    )


def test_skipped_variants(predict, caplog):
    df = predict(
        # outside of the genes
        ('chr1', 3000, 'A', 'C'),
        # ref differs from the FASTA file
        ('chr1', 6000, 'C', 'G'),
        # neither an SNV, an MNV, an insertion nor a deletion
        ('chr1', 6200, 'CCT', 'GA'),
    )
    assert_rows(df)
    assert [r.getMessage() for r in caplog.records] == [
        'Skipping variant chr1:6000:C>G: ref differs from the FASTA file (T)',
        'Skipping variant chr1:6200:CCT>GA: neither an SNV, an MNV, an insertion nor a deletion',
    ]


def test_read_gff3_genes(genome):
    genes = read_gff3_genes(genome / 'genes.gff3', s.TRANSCRIPT_TAGS).filter(pl.col('chrom') == 'chr1')
    expected = pl.DataFrame({
        'chrom': ['chr1'] * 4,
        'start': [5799, 6079, 6549, 7299],
        'end': [6600, 6900, 7200, 7600],
        'strand': ['+', '+', '-', '-'],
        'gene_id': ['PLUS1', 'PLUS2', 'MINUS1', 'MINUS2'],
        'sites': [[5800, 5900, 6020, 6110, 6500, 6600], [6080, 6150, 6300, 6350, 6800, 6900],
                  [6550, 6620, 6700, 6760, 7100, 7200], []],
    })
    pl.testing.assert_frame_equal(genes, expected)


def padded_forward(net, x, head):
    """Pangolin's forward pass: padded convolutions, then `flank` bases cropped on each side."""
    conv = net.conv1(x)
    skip = net.skip(conv)
    dense_convs = iter(net.convs)
    for i, block in enumerate(net.resblocks):
        dilation = block.conv1.dilation[0]
        padding = dilation * (block.conv1.kernel_size[0] - 1) // 2
        out = F.conv1d(torch.relu(block.bn1(conv)), block.conv1.weight, block.conv1.bias, padding=padding,
                       dilation=dilation)
        out = F.conv1d(torch.relu(block.bn2(out)), block.conv2.weight, block.conv2.bias, padding=padding,
                       dilation=dilation)
        conv = out + conv
        if i in net.skip_blocks:
            skip = skip + next(dense_convs)(conv)
    skip = skip[:, :, net.flank:skip.shape[2] - net.flank]
    return torch.softmax(net.heads[head](skip), dim=1)[:, 1]


def test_unpadded_net_matches_padded_net():
    # the size of Pangolin's network; 2 sequences with 30 outputs each
    net = make_net(0, n_channels=32)
    bases = np.random.default_rng(0).integers(0, 4, (2, 2 * net.flank + 30))
    x = torch.from_numpy(np.eye(4, dtype=np.float32)[bases].transpose(0, 2, 1).copy())
    with torch.inference_mode():
        for head in range(4):
            out = net(x, head)
            assert out.shape == (2, 30)
            np.testing.assert_allclose(out.numpy(), padded_forward(net, x, head).numpy(), rtol=0, atol=1e-6)
