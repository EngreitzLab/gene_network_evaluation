# organism from config
organism = config['organism']

rule download_genome_sizes:
    output: 
        "resources/genome_sizes/{organism}.txt"
    params:
        genome_sizes = config['genome_sizes']
    run:
        shell("mkdir -p resources/genome_sizes")
        shell("wget -c {params.genome_sizes} -O resources/genome_sizes/{organism}.txt")
