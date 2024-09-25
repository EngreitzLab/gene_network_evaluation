# Instructions for running dashboard for the 2024 IGVF Gene Program Jamboree
This page details how to run the gene program interpretation and annotation dashboard for the 2024 IGVF Gene Program Jamboree. For more information on gene program inference and evaluation, please see the [main repository README](../README.md).

# Overview
1. [Installation](#installation)
2. [Datasets](#datasets)
3. [Running the dashboard](#running-the-dashboard)

## Installation
Make sure you have a Python environment manager like [`pyenv`](https://github.com/pyenv/pyenv) or `conda` installed (if you don't have one installed, we recommend [installing miniconda](#install-miniconda)).

Once you are in your Python environment, install the required packages by running the following command:

```bash
cd gene_network_evaluation/app/
pip install -r requirements.txt
```

Let an organizer know if you run into any issues with the installation.

## Datasets
You will you will need a [Synapse account](https://www.synapse.org/Home:x) to access the large files considered for this jamboree.
First create an account, then post your username to the IGVF Slack channel **#jamboree-gene-programs-2024** to be granted access to these files. If you are interested in more general IGVF Synapse access, message Ben Hitz on slack with your Synapse username.

The dashboard expects that you have run the both gene inference and evaluation pipeline in advance (see main [`README`](../README.md) for more details). For the purposes of the jamboree, we have done this for you for 4 datasets (in alphabetical order):

1. `Bridge_Samples` -- IGVF left cerebral cortex brain samples mice. UCI performed snRNA-seq or "Split-seq" using the Parse Biosciences platform
2. `CharacterizationMcGinnis_Dataset6` -- H1 pancreatic diff, 7 timepoints (n=7; D3/7/11/18/32)
3. `Endothelial` -- TeloHAEC (endothelial model cell line) Perturb-seq dataset from [1]
4. `iPSC_EC` -- induced pluripotent stem cell (iPSC) derived endothelial cells (ECs) over a 3 day differentiation protocol

For each dataset we have run 3 program inference methods (`factor_analysis`, `cNMF`, and `Topyfic`, with a few exceptions) and several evaluations. The results are organized [on Synapse](https://www.synapse.org/Synapse:syn63392535) with the following structure:

- `datasets` - Contains the raw data for the datasets that can be input into the gene program inference methods
- `evaluation` - Contains the scripts and results for running the evaluation pipeline on the gene program inference results
- `get_data` - Contains scripts for downloading the data from Synapse
- `inference` - Contains the scripts and results for running the gene program inference methods
- `report` - Contains the scripts and results for generating the report

Each of these subdirectories contains a `README.md` and is subsequently organized by dataset.

Many of these files are quite large. Single files can be easily downloaded via the browser, but you can also download files programmatically using the [synapse command line client or Python API](https://pypi.org/project/synapseclient/). For example, you can run the following command to download the `Endothelial` evaluation `h5mu` file:

```bash
# Synapse may require you to generate a personal access token to authenticate. You can do so in your profile settings on Synapse.
synapse get syn63392881 --downloadLocation ../examples/evaluation/Endothelial/cNMF
```

Or you can follow the the [`examples/get_data/download_from_synapse.ipynb`](../examples/get_data/download_from_synapse.ipynb) notebook to download files using the Python API. Note that running the full notebook will download files for several datasets, so pick and choose the data you want to download when using this option.

## Running the dashboard
Running the dashboard requires that you have already run the evaluation pipeline (see [datasets](#datasets)) and that you pass in a `yaml` format config with the following information:

```yaml
path_evaluation_outs:  # List of paths to the evaluation output directories
  - "/path/to/evaluation/output/directory_1"
  - "/path/to/evaluation/output/directory_2"
path_evaluation_config: "/path/to/evaluation/config/evaluation.yaml"  # Path to the evaluation configuration file
path_mdata: "/path/to/eval.h5mu"  # Path to the MuData object containing the data_key and program_key
path_report_out: "/path/to/report/output/directory"  # Path to the directory where the report will be saved
data_key: "rna"  # Key used to store the original single-cell data in the MuData object
categorical_keys:   # List of categorical keys to use for dashboard
  - "sample"
continuous_keys:  # List of continuous keys to use for dashboard
  - "n_counts"
annotations_loc: "annotations.csv"  # Name of file in the report directory to dump annotations
```

We have provided an example config file for the [`Endothelial`](Endothelial.yaml) dataset.

```yaml
path_evaluation_outs: ["../examples/evaluation/Endothelial/cNMF"]
path_evaluation_config: "../examples/evaluation/Endothelial/cNMF/evaluation_pipeline.yml"
path_mdata: "../examples/evaluation/Endothelial/cNMF/cNMF.h5mu"
path_report_out: "../examples/report/Endothelial/cNMF"
prog_keys: ["cNMF"]
data_key: "rna"
categorical_keys: ["batch", "sample"]
continuous_keys: ["n_counts"]
annotations_loc: "annotations.csv"  # if none defaults to annotations.csv
```

All the files in the paths specified above were downloaded when you cloned the GitHub repository *except* for 1 file that is too large to store on GitHub. That file is available for [download from Synapse](https://www.synapse.org/Synapse:syn63392881). You can download it using the following command:

```bash
# Synapse may require you to generate a personal access token to authenticate. You can do so in your profile settings on Synapse.
synapse get syn63392881 --downloadLocation ../examples/evaluation/Endothelial/cNMF
```

We can now pass in `Endothelial.yaml` in as the `--config` argument of the following command:

```bash
python app.py --config Endothelial.yaml
```

Once the dashboard is running, you can access it via your default web browser by visiting [http://0.0.0.0:8050/](http://0.0.0.0:8050/).

If you are running the app on a server, you will need to open an SSH tunnel to access the dashboard. You can do this by running the following command in a terminal on your local machine:

```bash
ssh -L 8050:localhost:8050 <username>@<server>
```

Once the tunnel is open, you can access the dashboard by visiting [http://localhost:8050/](http://localhost:8050/) in your web browser.

# Dashboard features
We have hosted the Endothelial cell Perturb-seq example at the following link: http://34.169.124.229:8080/

# Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
Start by downloading the appropriate Miniconda installer for your operating system from https://docs.anaconda.com/miniconda/miniconda-other-installer-links/.

Run the installer using the command line. For example, if you downloaded the installer for Linux, you would run the following command:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Press enter when prompted

After hitting enter, you will then be presented with the license/terms and conditions. If you want to skip to the end, hit `q` and then accept the license terms by typing _yes_. 

Miniconda will then present you with an installation location. Press enter to confirm the location (or specify your own if you want to do a bit of organization).

Miniconda is now installing! This may take a bit so don't get frustrated. Leave your terminal open and let this run.

![image](https://github.com/user-attachments/assets/4808b7c3-5310-4f4d-b709-b926196e9329)

Finally, Miniconda will ask if you want to run `conda init` to configure your account to automatically use conda on each login. Type __yes__ and hit enter.

Time to check if this worked. Type:
    
```bash
source ~/.bashrc
```

Followed by:

```bash
conda --version
```

You should see the following output:

```bash
conda 24.5.0
```

You should also see your prompt change to something like:
```bash
(base) bash-5.1$
```

Now create an environment for jamboree by running the following command:
```bash
conda create -n jamboree-gene-programs python=3.10.8 -y
```

Activate (enter) the environment we just created:
```bash
conda activate jamboree-gene-programs
```

You can now go back to the [Installation](#installation) section to install the required packages.

# References
[1] Schnitzler, G. R., Kang, H., Fang, S., Angom, R. S., Lee-Kim, V. S., Ma, X. R., ... & Engreitz, J. M. (2024). Convergence of coronary artery disease genes onto endothelial cell programs. Nature, 626(8000), 799-807.

```bibtex
@ARTICLE{Schnitzler2024-xi,
  title     = "Convergence of coronary artery disease genes onto endothelial
               cell programs",
  author    = "Schnitzler, Gavin R and Kang, Helen and Fang, Shi and Angom,
               Ramcharan S and Lee-Kim, Vivian S and Ma, X Rosa and Zhou,
               Ronghao and Zeng, Tony and Guo, Katherine and Taylor, Martin S
               and Vellarikkal, Shamsudheen K and Barry, Aurelie E and
               Sias-Garcia, Oscar and Bloemendal, Alex and Munson, Glen and
               Guckelberger, Philine and Nguyen, Tung H and Bergman, Drew T and
               Hinshaw, Stephen and Cheng, Nathan and Cleary, Brian and Aragam,
               Krishna and Lander, Eric S and Finucane, Hilary K and
               Mukhopadhyay, Debabrata and Gupta, Rajat M and Engreitz, Jesse M",
  journal   = "Nature",
  publisher = "Springer Science and Business Media LLC",
  volume    =  626,
  number    =  8000,
  pages     = "799--807",
  abstract  = "Linking variants from genome-wide association studies (GWAS) to
               underlying mechanisms of disease remains a challenge1-3. For some
               diseases, a successful strategy has been to look for cases in
               which multiple GWAS loci contain genes that act in the same
               biological pathway1-6. However, our knowledge of which genes act
               in which pathways is incomplete, particularly for
               cell-type-specific pathways or understudied genes. Here we
               introduce a method to connect GWAS variants to functions. This
               method links variants to genes using epigenomics data, links
               genes to pathways de novo using Perturb-seq and integrates these
               data to identify convergence of GWAS loci onto pathways. We apply
               this approach to study the role of endothelial cells in genetic
               risk for coronary artery disease (CAD), and discover 43 CAD GWAS
               signals that converge on the cerebral cavernous malformation
               (CCM) signalling pathway. Two regulators of this pathway, CCM2
               and TLNRD1, are each linked to a CAD risk variant, regulate other
               CAD risk genes and affect atheroprotective processes in
               endothelial cells. These results suggest a model whereby CAD risk
               is driven in part by the convergence of causal genes onto a
               particular transcriptional pathway in endothelial cells. They
               highlight shared genes between common and rare vascular diseases
               (CAD and CCM), and identify TLNRD1 as a new, previously
               uncharacterized member of the CCM signalling pathway. This
               approach will be widely useful for linking variants to functions
               for other common polygenic diseases.",
  month     =  feb,
  year      =  2024,
  language  = "en"
}
```
