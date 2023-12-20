# Download HOCOMOCO
rule download_HOCOMOCO:
    output: 'resources/hocomoco_meme.meme'
    log: 'logs/3_0_downloaded_hocomoco.log'
    params: link = config['hocomoco']
    shell: "wget -c {params.link} -O resources/hocomoco_meme.meme"

# Convert coordinate to seq

# Run FIMO

# Run enrichment