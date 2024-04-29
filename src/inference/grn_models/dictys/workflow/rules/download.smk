import os
print(os.getcwd())

rule download:
	params:
		dir = config['working_dir']
	run:
		shell("mkdir -p {params.dir}")
		shell("wget -q -o /dev/null -O {params.dir}/motifs.motif 'https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif'")
		shell("dictys_helper genome_homer.sh hg38 {params.dir}/genome")
		shell("wget -q -o /dev/null -O {params.dir}/gene.gtf.gz http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz")
		shell("gunzip {params.dir}/gene.gtf.gz ")
		shell("dictys_helper gene_gtf.sh {params.dir}/gene.gtf {params.dir}/gene.bed")
