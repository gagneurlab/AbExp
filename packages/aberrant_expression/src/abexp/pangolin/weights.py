# The published weights of Pangolin (Zeng and Li, Genome Biology 2022, https://github.com/tkzeng/Pangolin): their file
# names, their SHA-256 sums and their download. The weights are GPL-3 and not part of this package, see README.md.
# This module needs only the standard library.
import hashlib
import logging
import shutil
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

# pinned commit of the Pangolin repository
PANGOLIN_COMMIT = '5cf94b8db938c658391b4305cd7ce33297d44ff7'
MODEL_URL = f'https://raw.githubusercontent.com/tkzeng/Pangolin/{PANGOLIN_COMMIT}/pangolin/models/{{}}'
# The published models of the ensemble, in pangolin/models of the Pangolin repository: 3 replicates for each of the
# 4 tissues. MODEL_FILES[head][replicate] is the file name; final.<replicate>.<2 * head>.3.v2 predicts with the
# output head 2 * head + 1 (conv_last1, 3, 5 or 7).
MODEL_FILES = tuple(
    tuple(f'final.{replicate}.{2 * head}.3.v2' for replicate in (1, 2, 3))
    for head in range(4)
)
# the SHA-256 sum of each model file at PANGOLIN_COMMIT
MODEL_SHA256 = {
    'final.1.0.3.v2': 'f0478fab173b75f7f7e9fe96688bad6c50fa4a46d70557f423b110caaf565501',
    'final.2.0.3.v2': 'c4c6bb4880fa6fb28b14182ae3ea0600edb07056158f55325b5e6e6e48fc9f26',
    'final.3.0.3.v2': 'ec685a6e7105a4486c1f89a005458a13deb3fe7171f13d434f4877e386d10676',
    'final.1.2.3.v2': '559c05de3e1ce65c2515ca3e92ef85edb0ec2e47686ca58060e25891ce06eb3a',
    'final.2.2.3.v2': '48758ba8b95eee9aa9feea52672ef06ca1b34111299c27f8a710f734d8b9aae5',
    'final.3.2.3.v2': '7cb576c2b24db4fdd6970c4ca4fb7c20ae1b1d8ae80645ebbe689848b5743129',
    'final.1.4.3.v2': 'c50b12e0c0af776d5674ca5e346493f8265783494d4df383364de9c1136657f6',
    'final.2.4.3.v2': 'e03303bed4fd6f135ec0f6c1b192cce954ea42d0646f44d17b4a6fbb2b1f610e',
    'final.3.4.3.v2': '9476d2e25520d7ff15bece0cd5d3b657e3b1dd3cc5fcab1d9c3b62bea7a0c5b6',
    'final.1.6.3.v2': '2aae563fa18a8a9b6699c6c96e0d32b8ec7543f8f805fb3bc9de77302cc9f66e',
    'final.2.6.3.v2': '7d3c0b1b2a60067b940dec315567874fbc8bcd322f1b7c76bf969f51f0f53f7f',
    'final.3.6.3.v2': '756e7721a382cace24e9bfea5b543af5623f2487d9a3efe7385e9c76367005fd',
}


def check_sha256(path, sha256):
    """Raise a ValueError if the SHA-256 sum of the file `path` differs from `sha256`."""
    with open(path, 'rb') as f:
        actual = hashlib.file_digest(f, 'sha256').hexdigest()
    if actual != sha256:
        raise ValueError(f'The SHA-256 sum of {path} is {actual}, not {sha256}')


def download_model(file, models_dir):
    """Download one published model of Pangolin, the file `file` of `MODEL_FILES`, into `models_dir`.

    The file comes from pangolin/models of the Pangolin repository at `PANGOLIN_COMMIT`, 2.9 MB. If it exists with the
    SHA-256 sum of `MODEL_SHA256`, it is kept. Otherwise it is downloaded to `<file>.part`, checked against its SHA-256
    sum and then renamed.

    Args:
      file: name of the model file, e.g. 'final.1.0.3.v2'
      models_dir: folder for the model file; created if missing

    Raises:
      ValueError: if the existing or the downloaded file has another SHA-256 sum. A downloaded file is deleted then.
      urllib.error.URLError: if the download fails
    """
    models_dir = Path(models_dir)
    models_dir.mkdir(parents=True, exist_ok=True)
    path = models_dir / file
    if path.exists():
        check_sha256(path, MODEL_SHA256[file])
        return
    url = MODEL_URL.format(file)
    logger.info('Downloading %s', url)
    part = models_dir / f'{file}.part'
    with urllib.request.urlopen(url, timeout=60) as response, open(part, 'wb') as f:
        shutil.copyfileobj(response, f)
    try:
        check_sha256(part, MODEL_SHA256[file])
    except ValueError:
        part.unlink()
        raise
    part.replace(path)


def download_models(models_dir):
    """Download the published models of Pangolin, the files of `MODEL_FILES`, into `models_dir`.

    Calls `download_model` for each file, so the files that exist with the right SHA-256 sum are kept.

    Args:
      models_dir: folder for the model files; created if missing

    Raises:
      ValueError: if an existing or a downloaded file has another SHA-256 sum. A downloaded file is deleted then.
      urllib.error.URLError: if a download fails
    """
    for files in MODEL_FILES:
        for file in files:
            download_model(file, models_dir)
