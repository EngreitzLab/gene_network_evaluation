# `report`

This directory contains the code and outputs for creating dashboards on the objects contained in the [`evaluation` directory](../evaluation/). Subdirectories are first broken up by dataset, and then by `prog_key`. 

The main outputs are `annotations.csv` files that contain program level annotations.

## Configuration files
Running the dashboard also makes use of a configuration file. `dashboard_pipeline.yaml` contains a template dashboard pipeline configuration file. This file can be copied and modified to build a dashboard for a given dataset and set of gene programs **that have already been evaluated with this pipeline**. It is theoretically possible to build a dashboard for gene programs that have not been evaluated with this pipeline, but this is not recommended. Running the dashboard is done as follows:

```bash
python app/app.py --config /path/to/config.yaml
```

Minimally, the configuration file should contain the following information:
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

An example dashboard configuration file can be found here [`examples/report/iPSC_EC/cNMF/cNMF_30/report.yaml`](/examples/report/iPSC_EC/cNMF/cNMF_30/report.yaml).