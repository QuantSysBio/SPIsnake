UNIPROT = features["reference"]["reference_proteome_fasta"]

rule transcript_sequences:
    input: 
        gff_compare="results/stringtie_assemble/gffcompare/gffcmp.combined.gtf",
        genome_fasta=features["reference"]["genome_fasta"],
        merged_gtf="results/stringtie_assemble/merged.gtf"
    output: 
        transcripts_fa="results/tr_2_prot/gffcmp.transcripts.fasta",
        gtf_final="results/tr_2_prot/merged_reference.transcripts.gtf"
    benchmark: 
        "results/benchmarks/transcript_sequences.txt"
    log: 
        "results/logs/transcript_sequences.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "gffread {input.gff_compare} \
        -T -o {output.gtf_final} \
        --no-pseudo \
        --force-exons \
        -M -Q \
        -w {output.transcripts_fa} \
        -g {input.genome_fasta} 2> {log}"


rule tr_names:
    input:
        gffcmp="results/tr_2_prot/gffcmp.transcripts.fasta",
        reference_tr_fasta=features["reference"]["transcriptome_fasta"]
    output:
        "results/tr_2_prot/merged_reference.transcripts.fasta"
    benchmark: 
        "results/benchmarks/tr_names.json"
    log: 
        "results/logs/tr_names.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    shell:
        "seqkit rmdup --by-seq \
        {input.reference_tr_fasta} {input.gffcmp} > {output}"


rule transdecoder_LongOrfs:
    input:
        fasta = "results/tr_2_prot/merged_reference.transcripts.fasta"    
    output:
        "results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    benchmark: 
        "results/benchmarks/transdecoder_LongOrfs.json"
    log: 
        "results/logs/transdecoder_LongOrfs.txt"
    conda: 
        "transdecoder.yaml"
    shell:
        "bin/TransDecoder-master/TransDecoder.LongOrfs \
        -t {input.fasta} \
        -O results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir \
        -m 100 &> {log}"


rule transdecoder_copyfiles:
    input: 
        pep = "results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: 
        pep="results/tr_2_prot/longest_orfs.pep",
        gff3="results/tr_2_prot/longest_orfs.gff3",
        cds="results/tr_2_prot/longest_orfs.cds",
        dat="results/tr_2_prot/base_freqs.dat"
    benchmark: 
        "results/benchmarks/transdecoder_movefiles.json"
    log: 
        "results/logs/transdecoder_movefiles.txt"
    params:
        gff3 = "results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/longest_orfs.gff3",
        cds = "results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/longest_orfs.cds",
        dat = "results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/base_freqs.dat"
    shell:
        "cp -f {input.pep} {output.pep} && \
        cp -f {params.gff3} {output.gff3} && \
        cp -f {params.cds} {output.cds} && \
        cp -f {params.dat} {output.dat}"
        

rule makeblastdb:
    input: 
        UNIPROT=UNIPROT
    output:
        [UNIPROT+'.pin', UNIPROT+'.phr', UNIPROT+'.psq']
    benchmark: 
        "results/benchmarks/makeblastdb.json"
    log: 
        "results/logs/makeblastdb.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "makeblastdb \
        -in {input.UNIPROT} \
        -dbtype prot \
        -logfile {log}"


rule blast:
    input: 
        longest_orfs = "results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        blast_ref_prot_db=[UNIPROT+'.pin', UNIPROT+'.phr', UNIPROT+'.psq']
    output: 
        "results/tr_2_prot/blastp.outfmt6"
    benchmark: 
        "results/benchmarks/blastp.json"
    log: 
        "results/logs/blastp.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "blastp \
        -num_threads {params.n} \
        -query {input.longest_orfs}  \
        -db {UNIPROT}  \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-15 \
        > {output} 2> {log}"


rule hmmpress:
    input: 
        "data/hmm/Pfam-A.hmm"
    output: 
        h3f="data/hmm/Pfam-A.hmm.h3f",
        h3i="data/hmm/Pfam-A.hmm.h3i",
        h3m="data/hmm/Pfam-A.hmm.h3m",
        h3p="data/hmm/Pfam-A.hmm.h3p"
    benchmark: 
        "results/benchmarks/hmmpress.txt"
    log: 
        "results/logs/hmmpress.txt"
    conda: 
        "hmmscan.yaml"
    params: 
        h3f="results/tr_2_prot/Pfam-A.hmm.h3f",
        h3i="results/tr_2_prot/Pfam-A.hmm.h3i",
        h3m="results/tr_2_prot/Pfam-A.hmm.h3m",
        h3p="results/tr_2_prot/Pfam-A.hmm.h3p",
        hmm="results/tr_2_prot/Pfam-A.hmm"
    shell: 
        "hmmpress {input} \
        2> {log} && \
        cp -f -t {params.h3f} {output.h3f} && \
        cp -f -t {params.h3i} {output.h3i} && \
        cp -f -t {params.h3m} {output.h3m} && \
        cp -f -t {params.h3p} {output.h3p} && \
        cp -f -t {params.hmm} data/hmm/Pfam-A.hmm"


rule hmm_db_copyfiles:
    input: 
        h3p="data/hmm/Pfam-A.hmm.h3p"
    output: 
        h3f="results/tr_2_prot/Pfam-A.hmm.h3f",
        h3i="results/tr_2_prot/Pfam-A.hmm.h3i",
        h3m="results/tr_2_prot/Pfam-A.hmm.h3m",
        h3p="results/tr_2_prot/Pfam-A.hmm.h3p",
        hmm="results/tr_2_prot/Pfam-A.hmm"
    benchmark: 
        "results/benchmarks/hmm_db_copyfiles.json"
    log: 
        "results/logs/hmm_db_copyfiles.txt"
    params:
        h3f="data/hmm/Pfam-A.hmm.h3f",
        h3i="data/hmm/Pfam-A.hmm.h3i",
        h3m="data/hmm/Pfam-A.hmm.h3m",
        h3p="data/hmm/Pfam-A.hmm.h3p",
        hmm="data/hmm/Pfam-A.hmm"
    shell:
        "cp -f {input.h3p} {output.h3p} && \
        cp -f {params.h3f} {output.h3f} && \
        cp -f {params.h3i} {output.h3i} && \
        cp -f {params.h3m} {output.h3m} && \
        cp -f {params.hmm} {output.hmm}"


rule hmmscan:
    input: 
        longest_orfs="results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        hmm="results/tr_2_prot/Pfam-A.hmm"
    output: 
        "results/tr_2_prot/hmmscan/pfam.domtblout"
    benchmark: 
        "results/benchmarks/hmmscan.txt"
    log: 
        "results/logs/hmmscan.txt"
    conda: 
        "hmmscan.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "hmmscan --cpu {params.n} \
        --domtblout {output} \
        results/tr_2_prot/Pfam-A.hmm \
        {input.longest_orfs} \
        &> {log}"



rule transdecoder_Predict:
    input: 
        transcripts_fa="results/tr_2_prot/merged_reference.transcripts.fasta",
        blastp="results/tr_2_prot/blastp.outfmt6",
        pfam="results/tr_2_prot/hmmscan/pfam.domtblout",
        pep="results/tr_2_prot/longest_orfs.pep"
    output: 
        transdecoder_pep="results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder.pep",
        transdecoder_gff3="results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder.gff3",
        bed="results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder.bed",
        cds="results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder.cds"
    benchmark: 
        "results/benchmarks/transdecoder_Predict.json"
    log: 
        "results/logs/transdecoder_Predict.txt"
    conda: 
        "transdecoder.yaml"
    params: 
        n=config["max_cores"],
        out_pep="merged_reference.transcripts.fasta.transdecoder.pep",
        out_gff3="merged_reference.transcripts.fasta.transdecoder.gff3",
        out_bed="merged_reference.transcripts.fasta.transdecoder.bed",
        out_cds="merged_reference.transcripts.fasta.transdecoder.cds"
    shell:
        "bin/TransDecoder-master/TransDecoder.Predict \
        -t {input.transcripts_fa} \
        --single_best_only \
        --retain_blastp_hits {input.blastp} \
        --retain_pfam_hits {input.pfam}\
        -O results/tr_2_prot 2> {log} 1>&2 && \
        cp -f {params.out_pep} {output.transdecoder_pep} && \
        cp -f {params.out_gff3} {output.transdecoder_gff3} && \
        cp -f {params.out_bed} {output.bed} && \
        cp -f {params.out_cds} {output.cds}"


rule gtf_to_alignment_gff3:
    input: 
        gtf_final="results/tr_2_prot/merged_reference.transcripts.gtf"
    output: 
        gff3_final="results/tr_2_prot/merged_reference.transcripts.gff3"
    benchmark: 
        "results/benchmarks/gtf_to_alignment_gff3.txt"
    log: 
        "results/logs/gtf_to_alignment_gff3.txt"
    conda: 
        "transdecoder.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "bin/TransDecoder-master/util/gtf_to_alignment_gff3.pl \
        {input.gtf_final} > {output.gff3_final} 2> {log}"


rule transcript_sequences_alignment_orf_to_genome_orf:
    input: 
        gff3_final="results/tr_2_prot/merged_reference.transcripts.gff3",
        transcripts_fa="results/tr_2_prot/merged_reference.transcripts.fasta",
        transdecoder_gff3="results/tr_2_prot/merged_reference.transcripts.fasta.transdecoder.gff3"
    output: 
        "results/tr_2_prot/transcripts.genome.gff3"
    benchmark: 
        "results/benchmarks/cdna_alignment_orf_to_genome_orf.txt"
    log: 
        "results/logs/cdna_alignment_orf_to_genome_orf.txt"
    conda: 
        "transdecoder.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "bin/TransDecoder-master/util/cdna_alignment_orf_to_genome_orf.pl {input.transdecoder_gff3} {input.gff3_final} {input.transcripts_fa} > {output} 2> {log}"


rule gff3_file_to_bed:
    input: 
        "results/tr_2_prot/transcripts.genome.gff3"
    output: 
        "results/tr_2_prot/proteome.bed"
    benchmark: 
        "results/benchmarks/gff3_file_to_bed.txt"
    log: 
        "results/logs/gff3_file_to_bed.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | bin/TransDecoder-master/util/gff3_file_to_bed.pl /dev/stdin | tail -n +2 > {output} 2> {log}"


rule gff3_file_to_proteins:
    input: 
        "results/tr_2_prot/transcripts.genome.gff3"
    output: 
        "results/tr_2_prot/proteome.fasta"
    benchmark: 
        "results/benchmarks/gff3_file_to_proteins.txt"
    log: 
        "results/logs/gff3_file_to_proteins.txt"
    conda: 
        "gff3_file_to_proteins.yaml"
    params:
        genome_fasta=features["reference"]["genome_fasta"]
    shell: "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {params.genome_fasta} | egrep -o '^[^*]+' > {output} 2> {log}"


rule seqkit_rmdup_orf:
    """
    remove duplivate entries in the predicted orfs
    """
    input:
        "results/tr_2_prot/proteome.fasta"
    output:
        "results/tr_2_prot/proteome.unique.fasta"
    benchmark: 
        "results/benchmarks/seqkit_rmdup_orf.txt"
    log: 
        "results/logs/seqkit_rmdup_orf.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "seqkit rmdup {input} \
        --by-seq --ignore-case \
        -o {output}"


rule seqkit_rmdup_final:
    """
    keep only novel sequences, matched to the reference transcripts
    the others will be included as they are in a reference proteome 
    """
    input:
        rmdup_orf = "results/tr_2_prot/proteome.unique.fasta",
        protein_all = features["reference"]["gencode_translation_sequences_fasta"]
    output:
        "results/tr_2_prot/proteome.final.fasta"
    benchmark: 
        "results/benchmarks/seqkit_rmdup_final.txt"
    log: 
        "results/logs/seqkit_rmdup_final.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "seqkit rmdup {input.protein_all} {input.rmdup_orf} \
        --by-seq --ignore-case \
        -o {output}"


rule tx2gene:
    input:
        gffcomp_tracking = "results/stringtie_assemble/gffcompare/gffcmp.tracking",
        transcripts_fa = "results/tr_2_prot/merged_reference.transcripts.fasta",
        gencode_annotation_gtf = features["reference"]["reference_gtf"]
    output:
        tx2gene = "results/tr_2_prot/tx2gene.csv"
    benchmark: 
        "results/benchmarks/tx2gene.txt"
    log: 
        "results/logs/tx2gene.txt"
    conda: 
        "expression_cutoffs_for_proteome_db.yaml"
    params: 
        n=config["max_cores"]
    script: 
        "tx2gene.R"


rule expression_cutoffs_for_proteome_db:
    """
    filters the proteome.fasta to keep the proteins with a minimum
    supporting number of reads in min_filt_samples
    """
    input:
        prot = "results/tr_2_prot/proteome.final.fasta",
        pred_prot = "results/RNAsamba_coding_potential/pred_prot.fasta",
        tx2gene = "results/tr_2_prot/tx2gene.csv"
    output:
        expr_prot="results/tr_2_prot/expr_prot.fasta",
        expr_prot_incl_deep="results/tr_2_prot/expr_prot_incl_deep.fasta",
        gene_tr_prot="results/tr_2_prot/gene_tr_prot.csv",
        gene_tr_prot_SP="results/tr_2_prot/gene_tr_prot_SP.csv",
        transcript_expression="results/tr_2_prot/transcript_expression.csv",
        transcript_expression_tpm="results/tr_2_prot/transcript_expression_tpm.csv",
        gene_expression="results/tr_2_prot/gene_expression.csv",
        gene_expression_tpm="results/tr_2_prot/gene_expression_tpm.csv"
    benchmark: 
        "results/benchmarks/expression_cutoffs_for_proteome_db.txt"
    log: 
        "results/logs/expression_cutoffs_for_proteome_db.txt"
    conda: 
        "expression_cutoffs_for_proteome_db.yaml"
    params: 
        n=config["max_cores"],
        min_filt_samples=config["min_filt_samples"],
        min_filt_counts=config["min_filt_counts"]
    script: 
        "expression_cutoffs.R"


rule seqkit_rmdup_extended:
    """
    removes the sequence duplicates from the final database
    """
    input:
        "results/tr_2_prot/expr_prot_incl_deep.fasta"
    output:
        "results/tr_2_prot/expr_prot_incl_deep_nodup.fasta"
    benchmark: 
        "results/benchmarks/seqkit_rmdup_extended.txt"
    log: 
        "results/logs/seqkit_rmdup_extended.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "seqkit rmdup {input} \
        --by-seq --ignore-case \
        -o {output}"