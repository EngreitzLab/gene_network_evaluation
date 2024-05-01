# organism from config
organism = config['organism']
if organism == "human":
    genome = "hg38"
elif organism == "mouse":
    genome = "mm10"
else:
    raise ValueError("Unknown organism: {organism}")
print(f"Organism: {organism}, Genome: {genome}")

rule download_genome_sizes:
    output: 
        "resources/genome_sizes/{organism}.txt"
    params:
        genome_sizes = config['genome_sizes']
    run:
        shell("mkdir -p resources/genome_sizes")
        shell("wget -c {params.genome_sizes} -O resources/genome_sizes/{organism}.txt")

rule download_genomes:
    output: 
        d=directory('resources/genomes/{organism}')
    shell:
        "python worfklow/scripts/download_genomes.py \
        -d {params.gdir} \
        -o {organism}"

