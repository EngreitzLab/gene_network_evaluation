import os
import pandas as pd

# Output directory from config
outdir: config['outdir']
currdir = os.getcwd()

rule prep_metacells:
    input:
        data=config['input_loc']
    output:
        path_metacell_stats="{outdir}/metacell_stats.csv",
        path_metacell_rds="{outdir}/metacells.rds",
        path_umap="{outdir}/umap.png",
        path_out="{outdir}/prepMetacells.rds"
    params:
        groupby=config['metacell_groupby'],
        genes=config['genes'],
        nearest_neighbors=config['nearest_neighbors'],
        reduction=config['reduction'],
        min_cells=config['min_cells'],
        target_metacells=config['target_metacells'],
        max_shared=config['max_shared'],
        name=config['name'],
        random_state=config['random_state']
    log:
        "{outdir}/logs/prep_metacells.log"
    benchmark:
        "{outdir}/benchmarks/prep_metacells.txt"
    shell:
        """
        Rscript workflow/scripts/prepMetacells.R \
        --path_data {input.data} \
        --path_output {output.path_out} \
        --path_metacell_stats {output.path_metacell_stats} \
        --path_metacell_rds {output.path_metacell_rds} \
        --path_umap {output.path_umap} \
        --group_by {params.groupby} \
        --genes {params.genes} \
        --reduction {params.reduction} \
        --min_cells {params.min_cells} \
        --target_metacells {params.target_metacells} \
        --max_shared {params.max_shared} \
        --nearest_neighbors {params.nearest_neighbors} \
        --name {params.name} \
        --seed {params.random_state}
        """

rule test_soft_powers:
    input:
        path_rds="{outdir}/prepMetacells.rds"
    output:
        path_power_table="{outdir}/powerTable.tsv",
        path_png="{outdir}/softPowers.png",
        path_out="{outdir}/testSoftPowers.rds"
    params:
        groupby=config['softpower_groupby'],
        groups=config['softpower_groups'],
        name=config['name'],
        random_state=config['random_state']
    log:
        "{outdir}/logs/test_soft_powers.log"
    benchmark:
        "{outdir}/benchmarks/test_soft_powers.txt"
    shell:
        """
        Rscript workflow/scripts/testSoftPowers.R \
        --path_rds {input.path_rds} \
        --path_output {output.path_out} \
        --path_table {output.path_power_table} \
        --path_png {output.path_png} \
        --group_by {params.groupby} \
        --groups {params.groups} \
        --name {params.name} \
        --seed {params.random_state}
        """

rule extract_soft_power:
    input:
        path_power_table="{outdir}/powerTable.tsv"
    output:
        path_soft_power="{outdir}/selected_softPower.txt"
    run:
        power_tbl = pd.read_csv(f"{config['outdir']}/powerTable.tsv", sep="\t")
        soft_power = power_tbl[power_tbl["SFT.R.sq"] > 0.8].iloc[0, 0]
        with open(output.path_soft_power, "w") as f:
            f.write(str(soft_power))
        print(f"Selected soft power: {soft_power}")

rule construct_network:
    input:
        path_rds="{outdir}/testSoftPowers.rds",
        soft_power_file="{outdir}/selected_softPower.txt"
    output:
        path_dendrogram="{outdir}/dendrogram.png",
        path_module_sizes="{outdir}/moduleSizes.tsv",
        path_out="{outdir}/constructNetwork.rds"
    params:
        name=config['name'],
        random_state=config['random_state']
    log:
        "{outdir}/logs/construct_network.log"
    benchmark:
        "{outdir}/benchmarks/construct_network.txt"
    shell:
        """
        Rscript workflow/scripts/constructNetwork.R \
        --path_rds {input.path_rds} \
        --path_output {output.path_out} \
        --path_dendrogram {output.path_dendrogram} \
        --path_module_sizes {output.path_module_sizes} \
        --power $(cat {input.soft_power_file}) \
        --name {params.name} \
        --seed {params.random_state}
        """
    
rule analyze_modules:
    input:
        path_rds="{outdir}/constructNetwork.rds"
    output:
        path_MEs="{outdir}/MEs.tsv",
        path_modules="{outdir}/modules.tsv",
        path_out="{outdir}/analyzeModules.rds"
    params:
        harmonize_by=config['modules_harmonize_by'],
        groupby=config['modules_groupby'],
        groups=config['modules_groups'],
        name=config['name'],
        random_state=config['random_state']
    log:
        "{outdir}/logs/analyze_modules.log"
    benchmark:
        "{outdir}/benchmarks/analyze_modules.txt"
    shell:
        """
        Rscript workflow/scripts/analyzeModules.R \
        --path_rds {input.path_rds} \
        --path_output {output.path_out} \
        --path_MEs {output.path_MEs} \
        --path_modules {output.path_modules} \
        --harmonize_by {params.harmonize_by} \
        --group_by {params.groupby} \
        --groups {params.groups} \
        --name {params.name} \
        --seed {params.random_state}
        """

rule prep_outputs:
    input:
        path_input=config['input_loc'],
        path_modules="{outdir}/modules.tsv",
        path_MEs="{outdir}/MEs.tsv",
    output:
        path_mdata="{outdir}/wgcna.h5mu",
    log:
        "{outdir}/logs/prep_outputs.log"
    benchmark:
        "{outdir}/benchmarks/prep_outputs.txt"
    shell:
        """
        python workflow/scripts/prep_outputs.py \
        --path_input {input.path_input} \
        --path_modules {input.path_modules} \
        --path_MEs {input.path_MEs} \
        --path_mdata {output.path_mdata}
        """


ruleorder: prep_metacells > test_soft_powers > extract_soft_power > construct_network > analyze_modules > prep_outputs

