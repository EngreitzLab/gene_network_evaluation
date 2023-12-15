# Download HOCOMOCO
rule download_HOCOMOCO:
    output: 'resources/hocomoco.meme'
    log: 'logs/3_0_downloaded_hocomoco.log'
    params: link = config['hocomoco']
    shell: "wget -c {params.link} -O resources/hocomoco.meme"

# Input seqs

# Run FIMO

# Run enrichment