# $$\color{red}\textbf{WORK} \space \textbf{IN} \space \textbf{PROGRESS}$$

# Snakemake pipeline to perform gene program evaluation

## Roadmap
1. Integrate outputs into a single mudata.
3. Implement plotting & interactive dashboard.
4. Enforce compatibility with standards
    * https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html
    

## OpenTargets Credential Generation Instructions
The disease-based benchmarking relies on the `resources/OpenTargets_L2G_Filtered.csv.gz` file, which is generated first by querying hundreds of GWAS from OpenTargets Genetics (output into `resources/OpenTargets_L2G_noQC.csv.gz`), and then by filtering the query down to the set of GWAS and Locus2Gene (L2G) scores used by the enrichment script. If you wish to re-generate the `resources/OpenTargets_L2G_Filtered.csv.gz` file, then you will need to provide Google Cloud credentials that will enable you to run an OpenTargets query through Google BigQuery.

If you have not done this before, the steps below explain how to setup a Google Cloud Service Account to use Google BigQuery.

Querying the first 1TB/month of public BigQuery data (including OpenTargets) is free, but there are fees beyond that, so you will need to link payment information to your credentials to run queries. Unless you are doing many very large queries, you should not reach the 1TB limit.

1. Go to https://console.cloud.google.com/
2. If you do not have a project already, create one.
3. Go to https://console.cloud.google.com/iam-admin/serviceaccounts & click `+ CREATE SERVICE ACCOUNT`
4. Enter a descriptive Service Account name like "query-opentargets" & set the role to "owner"
5. Once back on the Service Accounts main page, click the "Actions" menu button > "Manage Keys"
6. `ADD KEY` > `Create New Key` > Select`JSON`.
7. A `.json` file will automatically download. **Be careful with your `.json`. You will need to call this file, but DO NOT upload it to GitHub as these are your private credentials and malicious actors could use your credentials to rack up expensive queries**
8. Add the absolute path to the credentials .json to `smk/config/config.yaml` where it says `credentials:`

Once your service account is setup, you should be able to access data in the BigQuery database
`bigquery-public-data.open_targets_genetics`.  You can view the table schema in the OpenTargets Genetics database by going to 
https://console.cloud.google.com/bigquery?p=bigquery-public-data&d=open_targets_genetics&page=dataset .

Note that there is a separate database for the OpenTargets Platform `bigquery-public-data.open_targets_platform` that contains additional information that you may wish to map to the OpenTargets Genetics data you are pulling in, for example known drugs that target GWAS genes, although this database is not used in this codebase.
