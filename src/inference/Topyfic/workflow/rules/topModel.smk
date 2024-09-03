import pandas as pd
import Topyfic
from scripts import make_topmodel

configfile: 'config/config.yaml'

# Rule to run topModel for one input adata
rule run_top_model:
    input:
        expand(f"{config['workdir']}/{{name}}/{{n_topic}}/train/train_{{name}}_{{n_topic}}.p",
            name=config["names"],
            n_topic=config["n_topics"])
    output:
        f"{config['workdir']}/{{name}}/{{n_topic}}/topmodel/topModel_{{name}}_{{n_topic}}.p",
    params:
        name=lambda wildcards: wildcards.name,
        n_topic=lambda wildcards: wildcards.n_topic,
    run:
        # df, res = make_topmodel.find_best_n_topic(n_topics=config["n_topics"],
        #     names=params.name,
        #     topmodel_outputs=config['top_model']['workdir'])
        #
        # df.to_csv(f"{config['top_model']['workdir']}/k_N_{params.name}.csv",index=False)
        train = Topyfic.read_train(f"{config['workdir']}/{params.name}/{params.n_topic}/train/train_{params.name}_{params.n_topic}.p")

        make_topmodel.make_top_model(trains=[train],
            adata_paths=config['count_adata'][params.name],
            n_top_genes=int(config['top_model']['n_top_genes']),
            resolution=float(config['top_model']['resolution']),
            max_iter_harmony=int(config['top_model']['max_iter_harmony']),
            min_cell_participation=None if config['top_model']['min_cell_participation'] == "None" else float(
                config['top_model']['min_cell_participation']),
            topmodel_output=f"{config['workdir']}/{params.name}/{params.n_topic}/topmodel/")

# Rule to run topModel for merging multiple input adatas
if 'merge' in config:
    rule run_multiple_top_model:
        input:
            expand(f"{config['workdir']}/{{name}}/{{n_topic}}/train/train_{{name}}_{{n_topic}}.p",
                name=config["names"],
                n_topic=config["n_topics"]),
            expand(f"{config['workdir']}/{{name}}/{{n_topic}}/topmodel/topModel_{{name}}_{{n_topic}}.p",
               name=config["names"],
               n_topic=config["n_topics"]),
        output:
            f"{config['workdir']}/topModel_{'_'.join(config['names'])}.p",
            f"{config['workdir']}/analysis_{'_'.join(config['names'])}.p",
        params:
            name=lambda wildcards: wildcards.name,
            n_topic=lambda wildcards: wildcards.n_topic,
        run:
            df, res = make_topmodel.find_best_n_topic(n_topics=config["n_topics"],
                names=config["names"],
                topmodel_outputs=config['merging']['workdir'])

            df.to_csv(f"{config['workdir']}/k_N.csv",index=False)

            pd.DataFrame.from_dict(res, orient='index').to_csv(f"{config['workdir']}/best_k.csv")

            adata_paths = []
            for name in config['names']:
                adata_paths.append(config['count_adata'][name])

            trains = []
            for name in res.keys():
                n_topic = res[name]
                train = Topyfic.read_train(f"{config['workdir']}/{name}/{n_topic}/train/train_{name}_{n_topic}.p")
                trains.append(train)

            make_topmodel.make_top_model(trains=trains,
                adata_paths=adata_paths,
                n_top_genes=int(config['top_model']['n_top_genes']),
                resolution=float(config['top_model']['resolution']),
                max_iter_harmony=int(config['top_model']['max_iter_harmony']),
                min_cell_participation=None if config['top_model']['min_cell_participation'] == "None" else float(
                    config['top_model']['min_cell_participation']),
                topmodel_output=f"{config['workdir']}")
