import os
import re
import snakemake.io
from glob import glob

#Run BLASTN for sxt genes on all GLAMR (project set) samples (assemblies) of interest
rule make_sxtAGHI_db:
    input: 
        "blast/database/nucl/{blast_name}.fasta"
    output: 
        "blast/database/nucl/{blast_name}.fasta.nin"
    log: 
        "blast/logs/database/{blast_name}.txt"
    conda: "config/blast.yml"
    resources: 
        cpus=1, mem_mb=5000, time_min=10000
    shell: 
        """
        makeblastdb -in {input} -dbtype nucl -logfile {log}
        """

rule blastn_assemblies_vs_sxt:
    input: 
        blast_db = "blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        assembly = "blast/query/assemblies_vs_sxt/{sample}_final.contigs.renamed.fa"
    output:
         "blast/output/sxt/blastn/assemblies/{sample}__{blast_name}_contigs_blastn.tsv"
    log:
        "blast/logs/output/sxt/blastn/assemblies/{sample}__{blast_name}_contigs_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000
    shell:
        """
        blastn -db {input.blast_db} -query {input.assembly} -out {output} -outfmt '6 std qcovs stitle' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_assemblies_vs_sxt:
    input: 
        expand("blast/output/sxt/blastn/assemblies/{sample}__{blast_name}_contigs_blastn.tsv", sample = glob_wildcards("blast/query/assemblies_vs_sxt/{sample}_final.contigs.renamed.fa",followlinks=True).sample, blast_name = "blastn_db_sxt_AGHI_FINAL")

rule blastn_raw_bins_vs_sxt:
    input: 
        blast_db = "ada_bins/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        raw_bin = "blast/query/raw_bins_vs_sxt/{sample_binner_number}.fa"
    output:
         "blast/output/sxt/blastn/raw_bins/{sample_binner_number}__{blast_name}_raw_bin_blastn.tsv"
    log:
        "blast/logs/output/sxt/blastn/raw_bins/{sample_binner_number}__{blast_name}_raw_bin_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000
    shell:
        """
        blastn -db {input.blast_db} -query {input.raw_bin} -out {output} -outfmt '6 std qcovs stitle' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_raw_bins_vs_sxt:
    input: 
        expand("blast/output/sxt/blastn/raw_bins/{sample_binner_number}__{blast_name}_raw_bin_blastn.tsv", sample_binner_number = glob_wildcards("blast/query/raw_bins_vs_sxt/{raw_bin}.fa",followlinks=True).raw_bin, blast_name = "blastn_db_sxt_AGHI_FINAL")

rule drep_ada: 
    input: 
        bins = expand("ada_bins/raw_bins/{bin}.fa", bin = glob_wildcards("ada_bins/raw_bins/{bin}.fa").bin)
    output:
        main_dir = directory("ada_bins/drep_{ani}")
    conda: "config/drep.yml"
    resources: cpus=16, mem_mb=250000, time_min=2880,
    shell:
        """
  	    dRep dereplicate --S_ani 0.{wildcards.ani} {output.main_dir} -g {input.bins}
        """

rule run_ada_drep:
    input:
         expand("ada_bins/drep_{ani}", ani = ["95", "97", "98", "99"])

rule drep_ada_redo: 
    input: 
        bins = expand("ada_bins/raw_bins_ada_redo/{bin}.fa", bin = glob_wildcards("ada_bins/raw_bins_ada_redo/{bin}.fa").bin)
    output:
        main_dir = directory("ada_bins_ada_redo/drep_{ani}")
    conda: "config/drep.yml"
    resources: cpus=30, mem_mb=250000, time_min=2880,
    shell:
        """
  	    dRep dereplicate --S_ani 0.{wildcards.ani} {output.main_dir} -g {input.bins}
        """

rule run_ada_drep_ada_redo:
    input:
         expand("ada_bins_ada_redo/drep_{ani}", ani = ["98", "99"])

rule run_gtdbtk_GToTree:
    input:
        genomes = "ada_bins/gtdbtk/genomes_GToTree_99/",
        refs = "/geomicro/data2/kiledal/references/gtdbtk/release207_v2"
    output: directory("ada_bins/gtdbtk/output")
    conda: "config/gtdbtk.yml"
    resources: cpus=50
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fa --genome_dir {input.genomes} --out_dir {output} --cpus {resources.cpus}
        """

rule checkm_FINAL_MAGs:
    input: "FINAL_MAGs/"
    output:
        dir = directory("FINAL_MAGs/checkm"),
        results = "FINAL_MAGs/checkm.txt"
    conda: "config/checkm.yml"
    resources: cpus=16, mem_mb=80000, time_min=2880
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {input} {output.dir}
        """

rule drep_FINAL_MAGs: 
    output:
        main_dir = directory("FINAL_MAGs/drep_99")
    conda: "config/drep.yml"
    resources: cpus=16, mem_mb=250000, time_min=2880,
    shell:
        """
  	    dRep dereplicate --S_ani 0.99 {output.main_dir} -g FINAL_MAGs/*.fa
        """
        
rule gtdbtk_FINAL_MAGs:
    input:
        genomes = "FINAL_MAGs/",
        refs = "/geomicro/data2/kiledal/references/gtdbtk/release207_v2"
    output: directory("FINAL_MAGs/gtdbtk/output")
    conda: "config/gtdbtk.yml"
    resources: cpus=50
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fa --genome_dir {input.genomes} --out_dir {output} --cpus {resources.cpus}
        """

#Notebook3: sxtWXD and sxtA (positive control) presence/absence in sample 471 (primary bin: samp_471_concoct_280 - Dolichospermum sp000312705)
rule read_mapping_471_sxtWXDA:
    input: 
        f_reads = "mapping/reads/{sample}_sxtWXDA/reads/decon_fwd_reads_fastp.fastq.gz",
        r_reads = "mapping/reads/{sample}_sxtWXDA/reads/decon_rev_reads_fastp.fastq.gz",
        ref = "mapping/reads/{sample}_sxtWXDA/database/sxtWXDA.fasta"
    output: 
        sam = temp("mapping/reads/{sample}_sxtWXDA/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/{sample}_sxtWXDA/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/{sample}_sxtWXDA/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/{sample}_sxtWXDA/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/{sample}_sxtWXDA/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=40
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80

            samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
            samtools index -@ {resources.cpus} {output.bam}
        """

rule ref_read_mapping_pileup_471_sxtWXDA:
    input:
        bam = "mapping/reads/{sample}_sxtWXDA/output/bam/{sample}_{ref_seqs}_mapped.bam",
        ref = "mapping/reads/{sample}_sxtWXDA/database/{ref_seqs}.fasta"
    output:
        pileup = "mapping/reads/{sample}_sxtWXDA/output/pileup/{sample}_{ref_seqs}_pileup.txt"
    conda: "config/minimap2.yml"
    log: "mapping/reads/{sample}_sxtWXDA/output/log/{sample}_{ref_seqs}_pileup.log"
    resources: cpus=40
    shell:
        """
        samtools mpileup -f {input.ref} -o {output.pileup} {input.bam}
        """

rule run_read_mapping_471_sxtWXDA:
    input:
        expand("mapping/reads/{sample}_sxtWXDA/output/pileup/{sample}_{ref_seqs}_pileup.txt", sample = "samp_471", ref_seqs = "sxtWXDA")

rule blastn_GToTree_bins_vs_sxta: 
    input: 
        blast_db = "blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        tree_bin = "ada_bins_ada_redo/drep_99/GToTree/{sample}.fna"
    output:
         "blast/output/sxt/blastn/GToTree_bins/{sample}__{blast_name}_bin_blastn.tsv"
    log:
        "blast/logs/output/sxt/blastn/GToTree_bins/{sample}__{blast_name}_bin_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=24, mem_mb=16000, time_min=10000
    shell:
        """
        blastn -db {input.blast_db} -query {input.tree_bin} -out {output} -outfmt '6 std qcovs stitle' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_GToTree_bins_vs_sxta:
    input: 
        expand("blast/output/sxt/blastn/GToTree_bins/{sample}__{blast_name}_bin_blastn.tsv", sample = glob_wildcards("ada_bins_ada_redo/drep_99/GToTree/{sample}.fna",followlinks=True).sample, blast_name = "blastn_db_sxtA_FINAL")

rule read_mapping_GLAMRraw_sxtAll_mm2:
    input: 
        f_reads = "mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__fwd.fastq.gz",
        r_reads = "mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__rev.fastq.gz",
        ref = "mapping/reads/GLAMR_sxtAll_mm2/database/{ref_seqs}.fa"
    output: 
        sam = temp("mapping/reads/GLAMR_sxtAll_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/GLAMR_sxtAll_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/GLAMR_sxtAll_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/GLAMR_sxtAll_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/GLAMR_sxtAll_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_GLAMRraw_sxtAll_mm2:
    input:
        expand("mapping/reads/GLAMR_sxtAll_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__fwd.fastq.gz").sample, ref_seqs = "GLAMRsxtAll")

rule run_gtdbtk_GToTree_FINAL:
    input:
        genomes = "ada_bins_ada_redo/July23_2_Plus2073",
        refs = "/geomicro/data2/kiledal/references/gtdbtk/release207_v2"
    output: directory("ada_bins_ada_redo/July23_2_Plus2073/gtdbtk/output")
    conda: "config/gtdbtk.yml"
    resources: cpus=50
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fa --genome_dir {input.genomes} --out_dir {output} --cpus {resources.cpus}
        """

rule read_mapping_GLAMRraw_WLERef_mm2: 
    input: 
        f_reads = "mapping/reads/GLAMR_WLERef_mm2/reads/{sample}__fwd.fastq.gz",
        r_reads = "mapping/reads/GLAMR_WLERef_mm2/reads/{sample}__rev.fastq.gz",
        ref = "mapping/reads/GLAMR_WLERef_mm2/database/{ref_seqs}.fa"
    output: 
        sam = temp("mapping/reads/GLAMR_WLERef_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/GLAMR_WLERef_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/GLAMR_WLERef_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/GLAMR_WLERef_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/GLAMR_WLERef_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 80 \
            --minId 90
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_GLAMRraw_WLERef_mm2:
    input:
        expand("mapping/reads/GLAMR_WLERef_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/GLAMR_WLERef_mm2/reads/{sample}__fwd.fastq.gz").sample, ref_seqs = "WLERefCat3")

rule read_mapping_SpecificityCheck_mm2:
    input: 
        f_reads = "mapping/reads/MAG_specificity_mm2/reads/{sample}__fwd.fastq.gz",
        r_reads = "mapping/reads/MAG_specificity_mm2/reads/{sample}__rev.fastq.gz",
        ref = "mapping/reads/MAG_specificity_mm2/database/{ref_seqs}.fa"
    output: 
        sam = temp("mapping/reads/MAG_specificity_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/MAG_specificity_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/MAG_specificity_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/MAG_specificity_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/MAG_specificity_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 80 \
            --minId 90
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_SpecificityCheck_mm2_LE20_WE8:
    input:
        expand("mapping/reads/MAG_specificity_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/MAG_specificity_mm2/reads/{sample}__fwd.fastq.gz").sample, ref_seqs = "LE20WE8")

rule run_read_mapping_SpecificityCheck_mm2_LE16_WE8_Aug:
    input:
        expand("mapping/reads/MAG_specificity_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/MAG_specificity_mm2/reads/{sample}__fwd.fastq.gz").sample, ref_seqs = "LE16WE8Aug")

rule checkm_FINAL_MAGs_qc_rm:
    input: "FINAL_MAGs/qc_remove/"
    output:
        dir = directory("FINAL_MAGs/qc_remove/checkm"),
        results = "FINAL_MAGs/qc_remove/checkm.txt"
    conda: "config/checkm.yml"
    resources: cpus=16, mem_mb=80000, time_min=2880
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {input} {output.dir}
        """

rule prodigal_mags:
    input:
        bin = "prodigal/FINAL_MAGs_cp_Dec05_23/{bin}.fa"
    output:
        genes = "prodigal/{bin}.fasta",
        gbk = "prodigal/{bin}.gbk",
        proteins = "prodigal/{bin}.faa"
    conda: "config/main.yml"
    log: "prodigal/logs/{bin}.log"
    shell:
        """
        prodigal \
            -i {input.bin} \
            -a {output.proteins} \
            -d {output.genes} \
            -o {output.gbk} \
            1>{log} 2>&1
        """

rule run_prodigal_mags:
    input:
        expand("prodigal/{bin}.fasta", bin = glob_wildcards("prodigal/FINAL_MAGs_cp_Dec05_23/{mag}.fa").mag)

rule kofam_scan:
    input:
        proteins = rules.prodigal_mags.output.proteins,
        profile = "/geomicro/data2/kiledal/GLAMR/data/reference/kegg/kofamscan/profiles",
        ko_list = "/geomicro/data2/kiledal/GLAMR/data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "kofam_scan/{bin}_kofam_results.txt"
    conda: "config/kofamscan.yml"
    log: "kofam_scan/logs/{bin}.log"
    resources: cpus=64, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --format=detail-tsv \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.bin}_kofamscan \
            --ko-list {input.ko_list} {input.proteins}
        """
        
rule run_kofam_scan:
    input:
        expand("kofam_scan/{bin}_kofam_results.txt", bin = glob_wildcards("prodigal/FINAL_MAGs_cp_Dec05_23/{mag}.fa").mag)

rule KEGGdecoder_bins:
    input:
        "keggdec/input/{bin}_keggdec_input.tsv"
    output:
        "keggdec/output/{bin}_kegg_decoder_list.txt"
    conda: "keggdecoder"
    shell:
        """
        KEGG-decoder --input {input} --output {output} --vizoption static
        """

rule run_KEGGdecoder_bins:
    input:
        expand("keggdec/output/{bin}_kegg_decoder_list.txt", bin = glob_wildcards("prodigal/FINAL_MAGs_cp_Dec05_23/{mag}.fa").mag)

rule drep_ada_redo_2153: 
    input: 
        bins = expand("2153/raw_bins_ada_redo_2153/{bin}.fa", bin = glob_wildcards("2153/raw_bins_ada_redo_2153/{bin}.fa").bin)
    output:
        main_dir = directory("2153/drep_{ani}")
    conda: "config/drep.yml"
    resources: cpus=30, mem_mb=250000, time_min=2880,
    shell:
        """
  	    dRep dereplicate --S_ani 0.{wildcards.ani} {output.main_dir} -g {input.bins}
        """

rule run_drep_ada_redo_2153:
    input:
         expand("2153/drep_{ani}", ani = ["98", "99"])

#Notebook3: sxtWXD and sxtA (positive control) presence/absence in sample 471 (primary bin: samp_471_concoct_280 - Dolichospermum sp000312705)
rule read_mapping_471_sxtWXDA_REDO:
    input: 
        f_reads = "mapping/reads/{sample}_sxtWXDA/reads/decon_fwd_reads_fastp.fastq.gz",
        r_reads = "mapping/reads/{sample}_sxtWXDA/reads/decon_rev_reads_fastp.fastq.gz",
        ref = "mapping/reads/{sample}_sxtWXDA/database/sxtWXDA.fasta"
    output: 
        sam = temp("mapping/reads/{sample}_sxtWXDA_REDO/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/{sample}_sxtWXDA_REDO/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/{sample}_sxtWXDA_REDO/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/{sample}_sxtWXDA_REDO/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/{sample}_sxtWXDA_REDO/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=30
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 40 \
            --minId 80

            samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
            samtools index -@ {resources.cpus} {output.bam}
        """

rule ref_read_mapping_pileup_471_sxtWXDA_REDO:
    input:
        bam = "mapping/reads/{sample}_sxtWXDA_REDO/output/bam/{sample}_{ref_seqs}_mapped.bam",
        ref = "mapping/reads/{sample}_sxtWXDA/database/{ref_seqs}.fasta"
    output:
        pileup = "mapping/reads/{sample}_sxtWXDA_REDO/output/pileup/{sample}_{ref_seqs}_pileup.txt"
    conda: "config/minimap2.yml"
    log: "mapping/reads/{sample}_sxtWXDA_REDO/output/log/{sample}_{ref_seqs}_pileup.log"
    resources: cpus=30
    shell:
        """
        samtools mpileup -f {input.ref} -o {output.pileup} {input.bam}
        """

rule run_read_mapping_471_sxtWXDA_REDO:
    input:
        expand("mapping/reads/{sample}_sxtWXDA_REDO/output/pileup/{sample}_{ref_seqs}_pileup.txt", sample = "samp_471", ref_seqs = "sxtWXDA")

#Notebook3: sxtWXD and sxtA (positive control) presence/absence in sample_2153 
rule read_mapping_2153_sxtWXDA_REDO:
    input: 
        f_reads = "mapping/reads/{sample}_sxtWXDA_REDO/reads/decon_fwd_reads_fastp.fastq.gz",
        r_reads = "mapping/reads/{sample}_sxtWXDA_REDO/reads/decon_rev_reads_fastp.fastq.gz",
        ref = "mapping/reads/{sample}_sxtWXDA_REDO/database/sxtWXDA.fasta"
    output: 
        sam = temp("mapping/reads/{sample}_sxtWXDA_REDO/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/{sample}_sxtWXDA_REDO/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/{sample}_sxtWXDA_REDO/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/{sample}_sxtWXDA_REDO/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/{sample}_sxtWXDA_REDO/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=30
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 40 \
            --minId 80

            samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
            samtools index -@ {resources.cpus} {output.bam}
        """

rule ref_read_mapping_pileup_2153_sxtWXDA_REDO:
    input:
        bam = "mapping/reads/{sample}_sxtWXDA_REDO/output/bam/{sample}_{ref_seqs}_mapped.bam",
        ref = "mapping/reads/{sample}_sxtWXDA_REDO/database/{ref_seqs}.fasta"
    output:
        pileup = "mapping/reads/{sample}_sxtWXDA_REDO/output/pileup/{sample}_{ref_seqs}_pileup.txt"
    conda: "config/minimap2.yml"
    log: "mapping/reads/{sample}_sxtWXDA_REDO/output/log/{sample}_{ref_seqs}_pileup.log"
    resources: cpus=30
    shell:
        """
        samtools mpileup -f {input.ref} -o {output.pileup} {input.bam}
        """

rule run_read_mapping_2153_sxtWXDA_REDO:
    input:
        expand("mapping/reads/{sample}_sxtWXDA_REDO/output/pileup/{sample}_{ref_seqs}_pileup.txt", sample = "samp_2153", ref_seqs = "sxtWXDA")

rule drep_ALL: 
    input: 
        bins = expand("blast/query/raw_bins_vs_sxt/{bin}.fa", bin = glob_wildcards("blast/query/raw_bins_vs_sxt/{bin}.fa").bin)
    output:
        main_dir = directory("drep_all/drep_{ani}")
    conda: "config/drep.yml"
    resources: cpus=30, mem_mb=250000, time_min=2880,
    shell:
        """
  	    dRep dereplicate --S_ani 0.{wildcards.ani} {output.main_dir} -g {input.bins}
        """

rule run_drep_ALL:
    input:
         expand("drep_all/drep_{ani}", ani = ["98", "99"])

rule bakta_FINAL_MAGs:
    input:
        "prodigal/FINAL_MAGs_cp_Dec05_23/{bin}.fa"
    output:
        directory("bakta/output/{bin}")
    conda: "bakta"
    resources: cpus=48
    shell:
        """
        bakta --db /geomicro/data2/pdenuyl2/databases/bakta/db/ {input} --output {output} --keep-contig-headers --force --threads {resources.cpus}
        """

rule run_bakta_FINAL_MAGs:
    input:
        expand("bakta/output/{bin}", bin = glob_wildcards("prodigal/FINAL_MAGs_cp_Dec05_23/{mag}.fa").mag)

rule blastn_coassembly_bins_vs_sxt: 
    input: 
        blast_db = "blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        raw_bin = "blast/query/coassembly_bins_vs_sxt/{sample_binner_number}.fa"
    output:
        "blast/output/sxt/blastn/coassemblies/{sample_binner_number}__{blast_name}_coassembly_bin_blastn.tsv"
    log:
        "blast/logs/output/sxt/blastn/coassemblies/{sample_binner_number}__{blast_name}_coasssembly_bin_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000 
    shell: 
        """
        blastn -db {input.blast_db} -query {input.raw_bin} -out {output} -outfmt '6 std qcovs stitle' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_coassembly_bins_vs_sxt:
    input: 
        expand("blast/output/sxt/blastn/coassemblies/{sample_binner_number}__{blast_name}_coassembly_bin_blastn.tsv", sample_binner_number = glob_wildcards("blast/query/coassembly_bins_vs_sxt/{raw_bin}.fa",followlinks=True).raw_bin, blast_name = "blastn_db_sxt_AGHI_FINAL")

rule checkm_coassembly_initial:
    input: "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/FINAL_MAGs/sxt2/"
    output:
        dir = directory("/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/FINAL_MAGs/sxt2/checkm"),
        results = "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/FINAL_MAGs/sxt2/checkm.txt"
    conda: "config/checkm.yml"
    resources: cpus=16, mem_mb=80000, time_min=2880
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {input} {output.dir}
        """

rule run_gtdbtk_coassembly:
    input:
        genomes = "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/FINAL_MAGs/sxt2/",
        refs = "/geomicro/data2/pdenuyl2/databases/gtdbtk/release207_v2"
    output: directory("/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/FINAL_MAGs/sxt2/gtdbtk")
    conda: "config/gtdbtk.yml"
    resources: cpus=16
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fa --genome_dir {input.genomes} --out_dir {output} --cpus {resources.cpus}
        """

rule run_gtdbtk_sph_refs:
    input:
        genomes = "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/GToTree/GToTree_drep98_sxt2/Sphaerospermopsis_NCBI/selected_sequences_GTDBTK/",
        refs = "/geomicro/data2/pdenuyl2/databases/gtdbtk/release207_v2"
    output: directory("/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/GToTree/GToTree_drep98_sxt2/Sphaerospermopsis_NCBI/selected_sequences_GTDBTK/gtdbtk")
    conda: "config/gtdbtk.yml"
    resources: cpus=16
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fa --genome_dir {input.genomes} --out_dir {output} --cpus {resources.cpus}
        """

rule read_mapping_GLAMRraw_sxtWVX_mm2:
    input: 
        f_reads = "mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__fwd.fastq.gz",
        r_reads = "mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__rev.fastq.gz",
        ref = "mapping/reads/GLAMR_sxtWVX_mm2/database/{ref_seqs}.fa"
    output: 
        sam = temp("mapping/reads/GLAMR_sxtWVX_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/GLAMR_sxtWVX_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/GLAMR_sxtWVX_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/GLAMR_sxtWVX_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/GLAMR_sxtWVX_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_GLAMRraw_sxtWVX_mm2:
    input:
        expand("mapping/reads/GLAMR_sxtWVX_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__fwd.fastq.gz").sample, ref_seqs = "GLAMRsxtWVX")

rule read_mapping_GLAMRraw_sxtA_mm2:
    input: 
        f_reads = "mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__fwd.fastq.gz",
        r_reads = "mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__rev.fastq.gz",
        ref = "mapping/reads/GLAMR_sxtA_mm2/database/{ref_seqs}.fa"
    output: 
        sam = temp("mapping/reads/GLAMR_sxtA_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/GLAMR_sxtA_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/GLAMR_sxtA_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/GLAMR_sxtA_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/GLAMR_sxtA_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_GLAMRraw_sxtA_mm2:
    input:
        expand("mapping/reads/GLAMR_sxtA_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/GLAMR_sxtAll_mm2/reads/{sample}__fwd.fastq.gz").sample, ref_seqs = "LE20WE8sxtA")

rule assemble_biosyntheticSPAdes_samp_471:
    input:
        decon_reads_fwd = "/geomicro/data2/kiledal/GLAMR/data/omics/metagenomes/{run}/reads/decon_fwd_reads_fastp.fastq.gz",
        decon_reads_rev = "/geomicro/data2/kiledal/GLAMR/data/omics/metagenomes/{run}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        touch("/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/SPAdes/biosynthetic/{run}/.done"),
        assembly_dir = directory("/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/SPAdes/biosynthetic/{run}/")
    conda: "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/config/main.yml"
    log: "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/SPAdes/biosynthetic/log/{run}.log"
    resources: cpus = 24, time_min=20000, mem_mb = 500000, partition = "largemem"
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus}

        metaspades.py \
            -t {resources.cpus} \
            --bio \
            --memory $(({resources.mem_mb}/1024)) \
            -1 {input.decon_reads_fwd} \
            -2 {input.decon_reads_rev} \
            -o {output.assembly_dir} 2>&1 | tee {log}
        """

rule run_assemble_biosyntheticSPAdes_samp_471:
    input: 
        "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL2/SPAdes/biosynthetic/samp_471/.done" 

rule read_mapping_GLAMR_metatranscriptomes_sxtAll_mm2:
    input: 
        f_reads = "mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/reads/{sample}_decon_fwd_reads_fastp.fastq.gz",
        r_reads = "mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/reads/{sample}_decon_rev_reads_fastp.fastq.gz",
        ref = "mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/database/{ref_seqs}.fa"
    output: 
        sam = temp("mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_GLAMR_metatranscriptomes_sxtAll_mm2:
    input:
        expand("mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/GLAMR_metatranscriptomes_sxtAll_mm2/reads/{sample}_decon_fwd_reads_fastp.fastq.gz").sample, ref_seqs = "GLAMRsxtAll")

rule read_mapping_GLAMR_metatranscriptomes_sxtXA_mm2:
    input: 
        f_reads = "mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/reads/{sample}_decon_fwd_reads_fastp.fastq.gz",
        r_reads = "mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/reads/{sample}_decon_rev_reads_fastp.fastq.gz",
        ref = "mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/database/{ref_seqs}.fasta"
    output: 
        sam = temp("mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/output/{sample}_{ref_seqs}_mapped.sam"),
        temp_bam = temp("mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/output/{sample}_{ref_seqs}_mapped_temp.bam"),
        unsorted_bam = temp("mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/output/{sample}_{ref_seqs}_mapped_unsorted.bam"),
        bam = "mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam"
    conda: "config/minimap2.yml"
    log: "mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/output/log/{sample}_{ref_seqs}_bam.log"
    resources: cpus=16
    shell: 
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --sam-hit-only \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS -@ {resources.cpus} {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_read_mapping_GLAMR_metatranscriptomes_sxtXA_mm2:
    input:
        expand("mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/output/bam/{sample}_{ref_seqs}_mapped.bam", sample = glob_wildcards("mapping/reads/GLAMR_metatranscriptomes_sxtXA_mm2/reads/{sample}_decon_fwd_reads_fastp.fastq.gz").sample, ref_seqs = "sxtXA")

rule blastn_coassembly_13_vs_sxt_notebook9b: 
    input: 
        blast_db = "blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "blast/database/nucl/{blast_name}.fasta.nin",
        raw_bin = "blast/query/coassembly_13_vs_sxt_notebook9b/{sample_binner_number}.fa"
    output:
        "blast/output/sxt/blastn/coassembly_13_vs_sxt_notebook9b/{sample_binner_number}__{blast_name}_coassembly_bin_blastn.tsv"
    log:
        "blast/logs/output/sxt/blastn/coassemblies/{sample_binner_number}__{blast_name}_coasssembly_notebook9b_bin_blastn.log"
    conda: 
        "config/blast.yml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000 
    shell: 
        """
        blastn -db {input.blast_db} -query {input.raw_bin} -out {output} -outfmt '6 std qcovs stitle' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_coassembly_13_vs_sxt_notebook9b:
    input: 
        expand("blast/output/sxt/blastn/coassembly_13_vs_sxt_notebook9b/{sample_binner_number}__{blast_name}_coassembly_bin_blastn.tsv", sample_binner_number = glob_wildcards("blast/query/coassembly_13_vs_sxt_notebook9b/{raw_bin}.fa",followlinks=True).raw_bin, blast_name = "blastn_db_sxt_AGHI_FINAL")
