# Gene Program Evaluation Dashboard v0.0.1

## Overview
1. [Installation](#installation)
2. [Expected input](#expected-input)


## Installation
1. Install dash
```bash
pip install dash dash-bootstrap-components
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

## Dashboard

### Overview page
The overview page will describe both the inputs to the method and...



# 3) Install Miniconda

Installing multiple packages that are all compatible with each other can be a painful process for even seasoned bioinformaticians. Luckily for us, there exist package installation aides called ["package managers"](https://en.wikipedia.org/wiki/Package_manager) that make our lives a whole lot easier. [Many package managers exist](https://en.wikipedia.org/wiki/List_of_software_package_management_systems), but we will be using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for this bootcamp as its very flexible and lightweight.

Miniconda is pretty easy to install itself. Start by copying the the following file to your home directory:

```bash
scp /tscc/nfs/home/hkcarter/Miniconda3-latest-Linux-x86_64.sh ~/.
```

This file is a bash script (set of code instructions) that will install Miniconda on your account. To run the script, type the following command:
```bash
cd ~
bash Miniconda3-latest-Linux-x86_64.sh
```

Press enter when prompted

After hitting enter, you will then be presented with the license/terms and conditions. If you want to skip to the end, hit `q` and then accept the license terms by typing _yes_. 

Miniconda will then present you with an installation location as `/tscc/nfs/home/etrain##/miniconda3`. Press enter to confirm the location (or specify your own if you want to do a bit of organization -- see the pro tip at the end of this doc).

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

# 4. Setting up your "base" environment

Miniconda works by putting downloaded software into containers known as [environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments). This allows you to create different containers/environments that have different purposes.

When you first install Miniconda, you are put in a default environment called `base`. This is the environment that you are in now and will be whenever you log in to TSCC.

In base, we will need something called Jupyter notebooks for later in the bootcamp (don't worry if you don't know what those are for now). To install Jupyter, run the following command:

```bash
conda install -c conda-forge jupyter jupyterlab -y
```

<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    Up until recently, conda performance was painfully slow. This has since been remedied to some extent, but I'm not sure if it's been fixed entirely. A much faster alternative to conda is a package manager called <a href="https://mamba.readthedocs.io/">mamba</a> that is a drop-in replacement for conda.
    <br><br>
    <strong>To install mamba on an existing installation of conda, use the following:</strong>
    <br>
    <code>conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'</code>
    Once installed you can replace <code>conda</code> with <code>mamba</code> in any future command.
  </p>
</div>

# 5. Creating an environment for running an `rna-seq` analysis

During bootcamp, we will be honing our bioinformatic skills using an RNA-seq analysis. We will get more into the details in the first couple days, but for now, we need to install some software.

As a rule of thumb I install very little in my base environment (Jupyter being the exception). This is to avoid a bloated `base` environment that can cause performance issues.

Instead, we will create a new environment specifically for this bootcamp. Run the following command
```bash
conda create -n 2024-mstp-bootcamp python=3.11 r-base=4.3.1 -y
```

Let's break this down
- `conda create -n 2024-mstp-bootcamp` - create a new environment called `2024-mstp-bootcamp`
- `python=3.11` - install python version 3.11 in this environment
- `r-base=4.3.1` - install R version 4.3.1 in this environment
- `-y` - automatically say yes to any prompts

We can now 'activate' (enter) the environment we just created:
```bash
conda activate 2024-mstp-bootcamp
```

You should see your prompt change to something like:
```bash
(2024-mstp-bootcamp) bash-5.1$
```

This indicates that we are in the bootcamp environment, we can now install *most* of the necessary packages for RNA-seq analysis:

```bash
conda install -c conda-forge -c bioconda numpy pandas matplotlib seaborn STAR fastqc samtools bzip2 subread -y
```

Let's break this down
- `conda install -c conda-forge -c bioconda` - install packages from the conda-forge and bioconda channels
- `numpy pandas matplotlib seaborn` - install the python packages numpy, pandas, matplotlib, and seaborn
- `STAR fastqc samtools bzip2 subread` - install the programs STAR, fastqc, samtools, bzip2, and subread

Some packages are not available via conda and instead can be installed via the Python package manager [`pip`](https://pip.pypa.io/en/stable/).

Lucky for us, `pip` comes default when a new Python environment is created in conda, and conda and pip are very compatible. To install the packages we want, all we have to do is:
```bash
pip install decoupler pydeseq2 scanpy sanbomics gseapy PyWGCNA
```

Great! Hopefully these ran successfully for you. We will talk more about the packages and what they are used for in the actual bootcamp.

There is one last thing we need to do. Jupyter notebooks have no way of knowing where these programs are unless we tell them. We need to install something called ipykernel:
```bash
conda install -c anaconda ipykernel -y
```

and then create a "kernel" (Jupyter jargon) that knows where the software we just installed lives:
    
```bash
python -m ipykernel install --user --name 2024-mstp-bootcamp --display-name "Python 3.11 R 4.3.1 2024-mstp-bootcamp"
```

One last time for this notebook, let's break this down:
- `python -m ipykernel install` - run the command to install a new kernel
- `--user` - install the kernel for the current user only, as opposed to system-wide
- `--name 2024-mstp-bootcamp` - name the kernel `2024-mstp-bootcamp`, this should match the conda environment name
- `--display-name "Python 3.11 R 4.3.1 2024-mstp-bootcamp"` - display the kernel as "Python 3.11 R 4.3.1 2024-mstp-bootcamp" in Jupyter

You can now exit the interactive session by typing `exit` or `CTRL-D`.





# Install pyenv

Pyenv
brew install pyenv
brew install pyenv-virtualenv
pyenv install 3.10.8
pyenv virtualenv 3.10.8 jamboree-gene-programs




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