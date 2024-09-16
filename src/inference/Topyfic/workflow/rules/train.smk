from scripts import make_train

configfile: 'config/config.yaml'

# Rule to run make_single_train_model
rule run_single_train:
    output:
        f"{config['workdir']}/{{name}}/{{n_topic}}/train/train_{{name}}_{{n_topic}}_{{random_state}}.p",
    params:
        name=lambda wildcards: wildcards.name,
        n_topic=lambda wildcards: wildcards.n_topic,
        random_state=lambda wildcards: wildcards.random_state
    run:
        make_train.make_single_train_model(name=params.name,
            adata_path=config['count_adata'][params.name],
            k=int(params.n_topic),
            random_state=int(params.random_state),
            train_output=f"{config['workdir']}/{params.name}/{params.n_topic}/train/")

# Rule to run make_train_model
rule run_train_model:
    input:
        expand(f"{config['workdir']}/{{name}}/{{n_topic}}/train/train_{{name}}_{{n_topic}}_{{random_state}}.p",
            name=config["names"],
            n_topic=config["n_topics"],
            random_state=config["train"]["random_states"])
    output:
        f"{config['workdir']}/{{name}}/{{n_topic}}/train/train_{{name}}_{{n_topic}}.p",
    params:
        name=lambda wildcards: wildcards.name,
        n_topic=lambda wildcards: wildcards.n_topic,
    run:
        make_train.make_train_model(name=params.name,
            adata_path=config['count_adata'][params.name],
            k=int(params.n_topic),
            n_runs=config['train']['n_runs'],
            random_state=config['train']['random_states'],
            train_output=f"{config['workdir']}/{params.name}/{params.n_topic}/train/")