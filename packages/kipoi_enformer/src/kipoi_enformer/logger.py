import logging

__all__ = ['logger', 'setup_logger']

logger = logging.getLogger('kipoi_enformer')


def setup_logger(level=logging.INFO):
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger.setLevel(level)
    return logger
