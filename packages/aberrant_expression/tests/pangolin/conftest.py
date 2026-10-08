import torch

# pytest-xdist starts one worker per core, and PyTorch starts one thread per core in each worker. The workers then
# compete for the cores: on 16 cores, the tests take 340 s instead of 6 s. One thread per worker avoids that.
torch.set_num_threads(1)
