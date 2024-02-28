# organism from config
organism = config['organism']

rule download_blacklist:
    output: 
        "resources/blacklists/{organism}.bed"
    params:
        blacklist = config['blacklist']
    run:
        shell("mkdir -p resources/blacklists")
        shell("wget -c {params.blacklist} -O resources/blacklists/{organism}.bed.gz")
        shell("gunzip resources/blacklists/{organism}.bed.gz")
        
rule download_tf_annotations:
    output: 
        "resources/annotations/{organism}.tbl"
    params: 
        link = config['tf_annotations']
    run:
        shell("mkdir -p resources/annotations")
        shell("wget -c {params.link} -P resources/annotations")

rule download_tf_list:
    output: 
        "resources/tf_lists/{organism}.txt"
    params: 
        link = config['tf_list']
    run:
        shell("mkdir -p resources/tf_lists")
        shell("wget -c {params.link} -O resources/tf_lists/{organism}.txt")

rule download_databases:
    output: 
        "resources/rankings/{organism}.feather",
        "resources/scores/{organism}.feather"
    params:
        rankings = config['rankings_db'],
        scores = config['scores_db']
    run:
        shell("mkdir -p resources/rankings")
        shell("mkdir -p resources/scores")
        shell("wget -c {params.rankings} -P resources/rankings")
        shell("wget -c {params.scores} -P resources/scores")
