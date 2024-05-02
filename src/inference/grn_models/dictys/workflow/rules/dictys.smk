import os

print(os.getcwd())

rule setup_directory:
	input:
		data=f"{config['working_dir']}/{config['input_data']}"
	output:
		peak_file=expand(f"{config['working_dir']}/static/Subset{{cluster}}/all_peak.bed", cluster=config['clusters']),
		working_dir=directory(f"{config['working_dir']}/static"),
		cluster_dir=directory(expand(f"{config['working_dir']}/static/Subset{{cluster}}", cluster=config['clusters']))
	shell:
		"python workflow/scripts/setup_cell_cluster.py --mudata {input.data} --data_dir {output.working_dir}"

rule fragment_to_bam:
	input:
		all_peak=f"{config['working_dir']}/static/Subset{{cluster}}/all_peak.bed",
		frag=f"{config['working_dir']}/{config['fragment']}",
		working_dir=f"{config['working_dir']}/static"
	output:
		bam=f"{config['working_dir']}/static/Subset{{cluster}}/reads.bam",
		bai=f"{config['working_dir']}/static/Subset{{cluster}}/reads.bai"
	shell:
		"zcat {input.frag} | "
		"python workflow/scripts/fragment_to_bam_nofrag.py --cluster_number {wildcards.cluster} --data_dir {input.working_dir} | "
		"samtools view -b | samtools sort -o {output.bam} && samtools index {output.bam} {output.bai} " 
            
rule bam_to_overlapping_peak:
	input:
		bam=f"{config['working_dir']}/static/Subset{{cluster}}/reads.bam",
		bai=f"{config['working_dir']}/static/Subset{{cluster}}/reads.bai",
		peak_file=f"{config['working_dir']}/static/Subset{{cluster}}/all_peak.bed",
		working_dir=f"{config['working_dir']}/static/Subset{{cluster}}"
	output:
		macs_peaks=f"{config['working_dir']}/static/Subset{{cluster}}/reads_peaks.narrowPeak",
		overlap_peaks=f"{config['working_dir']}/static/Subset{{cluster}}/peaks.bed"
	shell:
		"macs2 callpeak -t {input.bam} -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --verbose 4 "
		"--call-summits -q 0.05 --name reads --outdir {input.working_dir} && " 
		"bedtools intersect -a {input.peak_file} -b {output.macs_peaks} -wa | bedtools sort -i - | uniq > {output.overlap_peaks}"

rule peak_to_footprints:
	input:
		peak_file=f"{config['working_dir']}/static/Subset{{cluster}}/peaks.bed",
		working_dir=f"{config['working_dir']}/static/Subset{{cluster}}",
		motif=f"{config['working_dir']}/{config['motif']}",
		genome=f"{config['working_dir']}/{config['genome']}",
		annotation=f"{config['working_dir']}/{config['gene_annotation']}"
	params:
		n_threads=config['n_threads']
	output:
		tf_to_gene_mask=f"{config['working_dir']}/static/Subset{{cluster}}/binlinking.tsv.gz"
	shell:
		"python workflow/scripts/chromatin_proc.py --subdir {input.working_dir} --motif {input.motif} "
		"--genome {input.genome} --gene_annotation {input.annotation} --threads {params.n_threads}"
        
rule infer_grn:
	input:
		tf_to_gene_mask=f"{config['working_dir']}/static/Subset{{cluster}}/binlinking.tsv.gz",
		directory=f"{config['working_dir']}/static/Subset{{cluster}}"
	params:
		n_threads=config['n_threads'],
		device=config['device']
	output:
		grn_weights=f"{config['working_dir']}/static/Subset{{cluster}}/net_iweight.tsv.gz"
	shell:
		"python workflow/scripts/infer_grn.py --subdir {input.directory} --device {params.device} --threads {params.n_threads}"
        
rule aggregate_results:
	input:
		data=f"{config['working_dir']}/{config['input_data']}",
		weights_files=expand(f"{config['working_dir']}/static/Subset{{cluster}}/net_iweight.tsv.gz", cluster=config['clusters']),
		input_dir=expand(f"{config['working_dir']}/static/Subset{{cluster}}", cluster=config['clusters'])
	output:
		output_dir=directory(config['output_dir']),
		output_mudata=f"{config['output_dir']}/dictys.h5mu"
	shell:
		"python workflow/scripts/aggregate_output.py --mudata {input.data} --input_dir {input.input_dir} --output_dir {output.output_dir}"
