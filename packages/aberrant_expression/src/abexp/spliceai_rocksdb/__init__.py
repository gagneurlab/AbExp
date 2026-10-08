"""SpliceAI scores from SpliceAI-RocksDB, and from SpliceAI for the variants that are not in the database.

Ported from spliceai_rocksdb (https://github.com/gagneurlab/spliceai_rocksdb, commit 3c40d6e), on kipoiseq2
instead of kipoiseq and cyvcf2. Needs the extra `rocksdb`. spliceai_rocksdb has no license file; its setup.py
declares license="MIT license" and author="M. Hasan Çelik". See LICENSE.
"""
from abexp.spliceai_rocksdb.spliceai import SpliceAI

__all__ = [
    'SpliceAI',
]
