import pyranges as pr
import tensorflow as tf
import pathlib


def gtf_to_pandas(gtf: str | pathlib.Path):
    """
    Read GTF file to pandas DataFrame
    :param gtf: Path to GTF file
    :return:
    """
    return pr.read_gtf(gtf, as_df=True, duplicate_attr=True)


class RandomModel(tf.keras.Model):
    """
    A random model for testing purposes.
    """

    def __init__(self, lamda=10):
        super().__init__()
        self.lamda = lamda

    def predict_on_batch(self, input_tensor):
        # tf.random.set_seed(42)
        return {
            'human': tf.abs(tf.random.poisson((input_tensor.shape[0], 896, 5313,), lam=self.lamda)),
            'mouse': tf.abs(tf.random.poisson((input_tensor.shape[0], 896, 1643), lam=self.lamda)),
        }
