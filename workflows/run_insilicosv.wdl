version 1.0

#
# Workflow to simulate genomes with SVs, simulate reads from them,
# map the reads to the reference, and visualize the results.
#

task sim_svs_insilicosv {
  input {
    File insilicosv_config_yaml
    String insilicosv_path
  }

  command <<<
    set -ex -o pipefail

    mkdir work
    cp ~{insilicosv_config_yaml} work/config.yaml
    PYTHONPATH=~{insilicosv_path} ~{insilicosv_path}/insilicosv/simulate.py \
        -c work/config.yaml
    cat work/sim.hapA.fa work/sim.hapB.fa > work/sim.fa
  >>>

  output {
    File sim_vcf = "work/sim.vcf"
    File sim_fa = "work/sim.fa"
  }
}

# sr: short reads
# lr: long reads

task sim_sr_dwgsim {
  input {
    File sim_fa
    String dwgsim_args
  }

  command <<< 
    set -ex -o pipefail
   
    dwgsim ~{dwgsim_args} -o 1 ~{sim_fa} sim_sr.dwgsim
  >>>

  output {
    File sim_sr_read1_fq_gz = "sim_sr.dwgsim.bwa.read1.fastq.gz"
    File sim_sr_read2_fq_gz = "sim_sr.dwgsim.bwa.read2.fastq.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/dwgsim:1.1.14--h50ea8bc_0"
  }
}

task samtools_sort {
  input {
    File mapped_sam
    Int threads = 1
    String mapped_basename = basename(mapped_sam, ".sam")
  }
  command <<<
    set -ex -o pipefail

    cat ~{mapped_sam} | samtools view -Sb - > "~{mapped_basename}.bam"
    samtools sort -@ ~{threads} "~{mapped_basename}.bam" > "~{mapped_basename}.sorted.bam"
    samtools index "~{mapped_basename}.sorted.bam"
  >>>
  output {
    File sorted_bam = "~{mapped_basename}.sorted.bam"
    File sorted_bam_bai = "~{mapped_basename}.sorted.bam.bai"
  }
  runtime {
    docker: "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    cpu: threads
  }
}

task map_sr_bwa {
  input {
    File ref_fa
    File sim_sr_read1_fq_gz
    File sim_sr_read2_fq_gz
    String bwa_args
    Int threads
  }

  command <<<
    set -ex -o pipefail
   
    bwa index ~{ref_fa}
    bwa mem ~{bwa_args} -t ~{threads} ~{ref_fa} ~{sim_sr_read1_fq_gz} ~{sim_sr_read2_fq_gz} \
        > sim_sr.bwa.sam
  >>>

  output {
    File sim_sr_bwa_sam = "sim_sr.bwa.sam"
  }
  runtime {
    docker: "quay.io/biocontainers/bwa:0.7.18--he4a0461_1"
    cpu: threads
  }
}

task sim_lr_pbsim3 {
  input {
    File sim_fa
    String pbsim3_args
  }

  command <<<
    set -ex -o pipefail
    
    pbsim --genome ~{sim_fa} --prefix sim_lr ~{pbsim3_args}
    cat sim_lr_*.fastq > sim_lr.fastq
  >>>

  output {
    File sim_lr_fastq = "sim_lr.fastq"
  }

  runtime {
    docker: "quay.io/biocontainers/pbsim3:3.0.4--h4ac6f70_0"
  }
}

task map_lr_minimap2 {
  input {
    File ref_fa
    File sim_lr_fastq
    String minimap2_args
    Int threads
  }

  command <<<
    set -ex -o pipefail

    minimap2 -eqx -ax map-hifi -t ~{threads} ~{minimap2_args} ~{ref_fa} ~{sim_lr_fastq} \
        > sim_lr.minimap2.sam
  >>>

  output {
    File sim_lr_minimap2_sam = "sim_lr.minimap2.sam"
  }

  runtime {
    docker: "quay.io/biocontainers/minimap2:2.28--he4a0461_3"
    cpu: threads
  }
}

task viz_mapped_reads_igvreports {
  input {
    File ref_fa

    File sim_vcf

    File sim_sr_bwa_bam
    File sim_sr_bwa_bam_bai

    File sim_lr_minimap2_bam
    File sim_lr_minimap2_bam_bai

    String igvreports_args = ""
  }

  command <<<
    set -eux -o pipefail

    create_report ~{sim_vcf} --fasta ~{ref_fa} \
         --tracks ~{sim_vcf} ~{sim_sr_bwa_bam} ~{sim_lr_minimap2_bam} \
        --output igvreport.html ~{igvreports_args}
  >>>

  output {
    File igvreport_html = "igvreport.html"
  }
  runtime {
    docker: "quay.io/biocontainers/igv-reports:1.12.0--pyh7cba7a3_0"
  }
}

workflow run_insilicosv_wf {
  input {
    File ref_fa

    File insilicosv_config_yaml
    String insilicosv_path

    String dwgsim_args
    String bwa_args

    String pbsim3_args
    String minimap2_args

    String igvreports_args = ""

    Int threads = 4
  }

  call sim_svs_insilicosv {
    input:
      insilicosv_config_yaml=insilicosv_config_yaml,
      insilicosv_path=insilicosv_path
  }

  call sim_sr_dwgsim {
    input:
    sim_fa = sim_svs_insilicosv.sim_fa,
    dwgsim_args = dwgsim_args
  }
  call map_sr_bwa {
    input:
    ref_fa = ref_fa,
    sim_sr_read1_fq_gz = sim_sr_dwgsim.sim_sr_read1_fq_gz,
    sim_sr_read2_fq_gz = sim_sr_dwgsim.sim_sr_read2_fq_gz,
    bwa_args = bwa_args,
    threads = threads
  }
  call samtools_sort as samtools_sort_sr {
    input:
    mapped_sam = map_sr_bwa.sim_sr_bwa_sam,
    threads = threads
  }

  call sim_lr_pbsim3 {
    input:
    sim_fa = sim_svs_insilicosv.sim_fa,
    pbsim3_args = pbsim3_args
  }
  call map_lr_minimap2 {
    input:
    ref_fa = ref_fa,
    sim_lr_fastq = sim_lr_pbsim3.sim_lr_fastq,
    minimap2_args = minimap2_args,
    threads = threads
  }
  call samtools_sort as samtools_sort_lr {
    input:
    mapped_sam = map_lr_minimap2.sim_lr_minimap2_sam,
    threads = threads
  }

  call viz_mapped_reads_igvreports {
    input:
    ref_fa = ref_fa,

    sim_vcf = sim_svs_insilicosv.sim_vcf,

    sim_sr_bwa_bam = samtools_sort_sr.sorted_bam,
    sim_sr_bwa_bam_bai = samtools_sort_sr.sorted_bam_bai,

    sim_lr_minimap2_bam = samtools_sort_lr.sorted_bam,
    sim_lr_minimap2_bam_bai = samtools_sort_lr.sorted_bam_bai,

    igvreports_args = igvreports_args
  }

  output {
     File sim_svs_sim_vcf = sim_svs_insilicosv.sim_vcf
     File sim_svs_sim_fa = sim_svs_insilicosv.sim_fa

     File sim_sr_bwa_bam = samtools_sort_sr.sorted_bam
     File sim_sr_bwa_bam_bai = samtools_sort_sr.sorted_bam_bai

     File sim_lr_minimap2_bam = samtools_sort_lr.sorted_bam
     File sim_lr_minimap2_bam_bai = samtools_sort_lr.sorted_bam_bai
    
     File viz_igvreport_html = viz_mapped_reads_igvreports.igvreport_html
  }
}
