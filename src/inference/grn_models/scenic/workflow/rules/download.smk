import os

# Get working dir
currdir = os.getcwd()

rule download_tf_annotations:
    output: 
        os.path.join(currdir, "resources", "motif_annotations", os.path.basename(config['motif_annotations']))
    params: 
        link = config['motif_annotations']
    run:
        shell("mkdir -p resources/motif_annotations")
        shell("wget -c {params.link} -O {output}")

rule download_tf_list:
    output: 
        os.path.join(currdir, "resources", "tf_lists", os.path.basename(config['tf_list']))
    params: 
        link = config['tf_list']
    run:
        shell("mkdir -p resources/tf_lists")
        shell("wget -c {params.link} -O {output}")

rule download_rankings_db:
    output: 
        os.path.join(currdir, "resources", "rankings_db", os.path.basename(config['rankings_db']))
    params:
        rankings = config['rankings_db'],
    run:
        shell("mkdir -p resources/rankings_db")
        shell("wget -c {params.rankings} -O {output}")
