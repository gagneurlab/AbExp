# The network of Pangolin (Zeng and Li, Genome Biology 2022, https://github.com/tkzeng/Pangolin), written for the
# published weights of Pangolin. The weights are GPL-3, see README.md.
from pathlib import Path

import numpy as np
import torch
from torch import nn

from abexp.pangolin.weights import MODEL_FILES, MODEL_SHA256, PANGOLIN_COMMIT, check_sha256

# the 16 residual blocks of Pangolin: kernel size and dilation of their two convolutions
KERNEL_SIZES = (11,) * 8 + (21,) * 4 + (41,) * 4
DILATIONS = (1,) * 4 + (4,) * 4 + (10,) * 4 + (25,) * 4
N_CHANNELS = 32
# after every 4th residual block, and after the last one, a 1x1 convolution adds the features to the skip connection
SKIP_INTERVAL = 4


class ResidualBlock(nn.Module):
    """Pre-activation residual block: batch norm, ReLU and convolution, twice, plus the input.

    The convolutions are unpadded, so the output is `crop` bases shorter than the input on each side.
    """

    def __init__(self, n_channels, kernel_size, dilation):
        super().__init__()
        self.bn1 = nn.BatchNorm1d(n_channels)
        self.conv1 = nn.Conv1d(n_channels, n_channels, kernel_size, dilation=dilation)
        self.bn2 = nn.BatchNorm1d(n_channels)
        self.conv2 = nn.Conv1d(n_channels, n_channels, kernel_size, dilation=dilation)
        self.crop = dilation * (kernel_size - 1)

    def forward(self, x):
        out = self.conv1(torch.relu(self.bn1(x)))
        out = self.conv2(torch.relu(self.bn2(out)))
        return out + x[:, :, self.crop:x.shape[2] - self.crop]


class PangolinNet(nn.Module):
    """Pangolin's network with unpadded convolutions.

    Pangolin pads its convolutions and then crops `flank` bases on each side of the output. The outputs that it keeps
    never see the padding. So unpadded convolutions, whose features shrink from block to block, give the same outputs
    up to the order of the float sums, with less computation. The module and parameter names are those of Pangolin's
    state dicts.

    Args:
      n_channels: number of channels of all hidden layers
      kernel_sizes: kernel size of each residual block
      dilations: dilation of each residual block
    """

    def __init__(self, n_channels=N_CHANNELS, kernel_sizes=KERNEL_SIZES, dilations=DILATIONS):
        super().__init__()
        self.conv1 = nn.Conv1d(4, n_channels, 1)
        self.skip = nn.Conv1d(n_channels, n_channels, 1)
        self.resblocks = nn.ModuleList(
            ResidualBlock(n_channels, kernel_size, dilation) for kernel_size, dilation in zip(kernel_sizes, dilations)
        )
        n_blocks = len(self.resblocks)
        self.skip_blocks = {i for i in range(n_blocks) if (i + 1) % SKIP_INTERVAL == 0 or i + 1 == n_blocks}
        self.convs = nn.ModuleList(nn.Conv1d(n_channels, n_channels, 1) for _ in self.skip_blocks)
        # Four output heads, each a softmax over 2 channels (conv_last1, 3, 5 and 7), and four sigmoid heads
        # (conv_last2, 4, 6 and 8). Only channel 1 of the softmax heads is used; the sigmoid heads are in the state
        # dicts only.
        self.conv_last1 = nn.Conv1d(n_channels, 2, 1)
        self.conv_last2 = nn.Conv1d(n_channels, 1, 1)
        self.conv_last3 = nn.Conv1d(n_channels, 2, 1)
        self.conv_last4 = nn.Conv1d(n_channels, 1, 1)
        self.conv_last5 = nn.Conv1d(n_channels, 2, 1)
        self.conv_last6 = nn.Conv1d(n_channels, 1, 1)
        self.conv_last7 = nn.Conv1d(n_channels, 2, 1)
        self.conv_last8 = nn.Conv1d(n_channels, 1, 1)
        self.heads = (self.conv_last1, self.conv_last3, self.conv_last5, self.conv_last7)
        # bases of context on each side of an output
        self.flank = sum(block.crop for block in self.resblocks)

    def forward(self, x, head):
        """The splice site usage of output head `head` (0 to 3).

        Args:
          x: one-hot encoded sequences, float tensor of shape (batch, 4, length)
          head: index of the output head

        Returns:
          tensor of shape (batch, length - 2 * flank)
        """
        conv = self.conv1(x)
        skip = self.skip(conv)
        # bases cropped from each side of `conv` and of `skip`
        conv_crop = 0
        skip_crop = 0
        dense_convs = iter(self.convs)
        for i, block in enumerate(self.resblocks):
            conv = block(conv)
            conv_crop += block.crop
            if i in self.skip_blocks:
                crop = conv_crop - skip_crop
                skip = skip[:, :, crop:skip.shape[2] - crop] + next(dense_convs)(conv)
                skip_crop = conv_crop
        return torch.softmax(self.heads[head](skip), dim=1)[:, 1]


class PangolinModels:
    """The ensemble of Pangolin: several replicates of a network for each output head.

    Args:
      nets: `nets[head][replicate]` predicts with output head `head`. All nets need the same flank.
    """

    def __init__(self, nets):
        self.nets = tuple(tuple(replicates) for replicates in nets)
        flanks = {net.flank for replicates in self.nets for net in replicates}
        if len(flanks) != 1:
            raise ValueError(f'The nets have different flanks: {sorted(flanks)}')
        self.flank = flanks.pop()
        self.device = next(self.nets[0][0].parameters()).device

    @classmethod
    def from_dir(cls, models_dir, device=None):
        """Load the published models of Pangolin, the files of `MODEL_FILES`, from `models_dir`.

        `download_models` downloads the files. Each file must have the SHA-256 sum of `MODEL_SHA256`, so that the files
        of another Pangolin commit, e.g. from an older version of this package, do not load.

        Args:
          models_dir: folder with the model files
          device: torch device. Default is the GPU if torch finds one, as in Pangolin, else the CPU.

        Raises:
          ValueError: if a file has another SHA-256 sum
        """
        if device is None:
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
        nets = []
        for files in MODEL_FILES:
            replicates = []
            for file in files:
                path = Path(models_dir) / file
                try:
                    check_sha256(path, MODEL_SHA256[file])
                except ValueError as e:
                    raise ValueError(f'{e}, the sum at Pangolin commit {PANGOLIN_COMMIT}. Delete the folder '
                                     f'{models_dir} and download the files again.') from e
                net = PangolinNet()
                net.load_state_dict(torch.load(path, map_location=device, weights_only=True))
                replicates.append(net.to(device).eval())
            nets.append(replicates)
        return cls(nets)

    def predict(self, x):
        """The splice site usage of each head and replicate.

        Args:
          x: one-hot encoded sequences, float32 array of shape (batch, 4, length)

        Returns:
          float32 array of shape (heads, replicates, batch, length - 2 * flank)
        """
        x = torch.from_numpy(np.ascontiguousarray(x)).to(self.device)
        with torch.inference_mode():
            return np.stack([
                np.stack([net(x, head).cpu().numpy() for net in replicates])
                for head, replicates in enumerate(self.nets)
            ])
