import os
# workflow/rules/40_index.smk

rule generate_index:
    """
    Build a STAR genome index once (sjdbOverhang & genomeSAindexNbases come from config).
    """
    conda: os.path.join(config["envs_dir"], "star.yaml")
    input:
        genome_fasta = config["genome_fasta"],
        gtf          = config["gtf"]
    output:
        index_dir = directory(index_dir)
    threads: threads
    params:
        sjdb_overhang = sjdb_overhang,
        sa_nbases     = sa_nbases,
        script        = scripts["generate_index"]
    shell:
        r"""
        bash "{params.script}" \
          "{input.genome_fasta}" \
          "{input.gtf}" \
          "{output.index_dir}" \
          {params.sjdb_overhang} \
          {params.sa_nbases} \
          {threads}
        """
