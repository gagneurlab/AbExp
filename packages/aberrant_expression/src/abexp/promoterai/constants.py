# The input length of the published PromoterAI models, the `--input_length` of their scoring script. The models take
# any sequence length, but they were trained and published with this one.
DEFAULT_INPUT_LENGTH = 20_480

# PromoterAI rounds the difference of its outputs to this many decimals before the tanh, see PromoterAI.predict.
SCORE_DECIMALS = 4

# The bound on the positions (sequences x bases) per call of the model. The MLP of each block of the published models has
# 4 * 1024 = 4096 units, and a TensorFlow GPU kernel fails on a tensor of 2**31 elements or more ("Check failed:
# work_element_count >= 0"), so a call must hold fewer than 2**31 / 4096 = 2**19 positions, e.g. at most 25 sequences
# of 20480 bases. PromoterAI.predict splits larger batches into several calls.
MAX_POSITIONS_PER_CALL = 2 ** 19
