# `inference`

This directory contains the code and outputs for runnning gene program inference on the objects contained in the [`datasets` directory](../datasets/). Subdirectories are first broken up by dataset, with each containing the outputs of the 3 inference methods: factor_analysis, cNMF and Topyfic.

The main outputs of inference are again `h5mu` files mirroring the input format, but with the addition of the inferred gene programs. These files are then used as input to the evaluation and dashboard steps. See the [main README.md](../../README.md) for more information on the input and output formats for inference.

# Note:
Currently, the inference method procedure is not completely standardized across each method. You will likely need to dig into the code for each method to get it running on a new dataset. This is a known issue and will be addressed in future releases.