# Download HOCOMOCO
rule download_HOCOMOCO:
    output: 'resources/hocomoco_meme.meme'
    log: 'logs/3_0_downloaded_hocomoco.log'
    params: link = config['hocomoco']
    shell: "wget -c {params.link} -O resources/hocomoco_meme.meme"

# Perform enhancer motif matching and enrichment
rule run_enhancer_motif_enrichment:
    input: 
        'evaluation_mdata.h5mu',
        'resources/hocomoco_meme.meme'
    output: 
        touch('touch/3_1_run_enhancer_motif_enrichment.done')
    params:
        seq_class = 'enhancer',
        sig = 0.05
    log: 'logs/3_1_run_enhancer_motif_enrichment.log'
    script: '../scripts/3_run_motif_enrichment.py'

# Perform promoter motif matching and enrichment
rule run_promoter_motif_enrichment:
    input: 
        'evaluation_mdata.h5mu',
        'resources/hocomoco_meme.meme'
    output: 
        touch('touch/3_2_run_promoter_motif_enrichment.done')
    params:
        seq_class = 'promoter',
        sig = 0.05
    log: 'logs/3_2_run_promoter_motif_enrichment.log'
    script: '../scripts/3_run_motif_enrichment.py'
