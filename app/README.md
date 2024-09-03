# Gene Program Evaluation Dashboard v0.0.1

## Overview
1. [Installation](#installation)
2. [Expected input](#expected-input)


## Installation
1. Install dash
```bash
pip install dash dash-bootstrap-components dash-dangerously-set-inner-html
```

2. Run the app
```bash
python app.py
```

3. Open the app in your browser: 
If you are running the app locally, you can open it in your browser by visiting [http://](http://0.0.0.0:8050/)
If you are running the app on a server, run
```bash
ssh -L 8050:localhost:8050 <username>@<server>
```

## Expected input
The dashboard expects that you have run the snakemake evaluation pipeline and that the output is structured as follows:
    
```bash
path_pipeline_outs/
├── cNMF_60
├── cNMF_59
├── cNMF_58
...
├── cNMF_1
└── output.h5mu
```

where `output.h5mu` contains the MuData object with keys `input` that contains the original input matrix and `cNMF_60`, `cNMF_59`, ..., `cNMF_1` containing the cell and gene loadings for each of the inference runs.

The keys of the MuData and the output subdirectories in `path_pipeline_outs` should be named accordingly to the desired type of analysis:
1. A **single run of a method** with a fixed number of components. The MuData should have 1 key other than `input` that can be named anything. In this scenario, a cross run analysis (see below) will be omitted.
2. A **cross k analysis** should share the same base name and be suffixed by an integer that indicates the number of components used in inference. e.g. `cNMF_60`, `cNMF_59`, ..., `cNMF_1`. Here the assumption is that the number of components is varied across the same method.
3. A **cross method analysis** should have a key for each method. e.g. `cNMF`, `Topyfic`, etc. Here the assumption is that the number of components is fixed and the method is varied.
4. A **cross k and method** analysis is a combination of 2) and 3). The keys should have the same base name for each method and be suffixed by an integer that indicates the number of components used in inference. e.g. `cNMF_60`, `cNMF_59`, ..., `cNMF_1`, `Topyfic_60`, `Topyfic_59`, ..., `Topyfic_1`.

Note that if you run different methods with different numbers of components and don't suffix the keys with the number of components, the dashboard will not perform any comparisons across values of k.

## Dashboard

### Overview page
The overview page will describe both the inputs to the method and...