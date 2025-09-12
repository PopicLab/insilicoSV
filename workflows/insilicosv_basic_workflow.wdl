version 1.0

task insilicosv {
  input {
    File yaml
    String outDir
  }
  String yamlBasename = basename(yaml)
  command <<<
    set -e
    cp ~{yaml} ~{outDir}/.
    insilicosv -c ~{outDir}/~{yamlBasename}
    cat ~{outDir}/sim.hapA.fa ~{outDir}/sim.hapB.fa > ~{outDir}/sim.fa
  >>>

  output {
    File vcf = "~{outDir}/sim.vcf"
    File fa = "~{outDir}/sim.fa"
  }
}

task dwgsim {
  input {
    String outDir
    File fa
    Int hapCoverage
    String defaultArgs = "-1 150 -2 150 -c 0 -S 0 -y 0 -r 0"
    String? customArgs 
  }

  String args = if defined(customArgs) then customArgs else defaultArgs 
  command <<< 
    set -e
    mkdir -p ~{outDir}/short
    dwgsim ~{args} -C ~{hapCoverage} -H ~{fa} ~{outDir}/short/sim.dwgsim
  >>>

  output {
    File R1 = "~{outDir}/short/sim.dwgsim.bwa.read1.fastq.gz"
    File R2 = "~{outDir}/short/sim.dwgsim.bwa.read2.fastq.gz"
  }
}

task pbsim {
  input {
    File fa
    String outDir
    Int hapCoverage
    String platform
    String? customArgs
  }

  command <<<
    set -e 
    mkdir -p ~{outDir}/~{platform}

    if [[ ~{platform} == "hifi" ]]; then
       args="--method qshmm --qshmm $CONDA_PREFIX/data/QSHMM-RSII.model --length-mean 20000 --accuracy-mean 0.999"
    elif [[ ~{platform} == "ont" ]]; then
       args="--method qshmm --qshmm ${CONDA_PREFIX}/data/QSHMM-ONT-HQ.model --length-mean 50000 --accuracy-mean 0.98"
    elif [[ ~{platform} == "custom" ]]; then
       args="~{customArgs}"
    fi
 
    pbsim --strategy wgs $args --depth ~{hapCoverage} --genome ~{fa} --prefix ~{outDir}/~{platform}/sd
    zcat ~{outDir}/~{platform}/sd*.fq.gz > ~{outDir}/~{platform}/sim.pbsim.fastq
    rm ~{outDir}/~{platform}/sd*
  >>>

  output {
    File reads = "~{outDir}/~{platform}/sim.pbsim.fastq"
  }
}

task minimap2 {
  input {
    File ref
    File? R1
    File? R2
    File? R
    String platform
    String? customPreset
    String outDir
    Int threads
  }

  command <<<
    set -e

    if [[ ~{platform} == "hifi" ]]; then
        minimap2 --eqx -Y -ax map-hifi -t ~{threads} ~{ref} ~{R} | samtools sort -@ ~{threads} -o ~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam -
    elif [[ ~{platform} == "ont" ]]; then
        minimap2 --eqx -Y -ax map-ont -t ~{threads} ~{ref} ~{R} | samtools sort -@ ~{threads} -o ~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam -
    elif [[ ~{platform} == "short" ]]; then
        minimap2 --eqx -Y -ax sr -t ~{threads} ~{ref} ~{R1} ~{R2} | samtools sort -@ ~{threads} -o ~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam -
    elif [[ ~{platform} == "custom" ]]; then
        minimap2 --eqx -Y ~{customPreset} -t ~{threads} ~{ref} ~{R} | samtools sort -@ ~{threads} -o ~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam -
    fi
    samtools index ~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam
  >>>

  output {
    File bam = "~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam"
    File bai = "~{outDir}/~{platform}/sim.~{platform}.minimap2.sorted.bam.bai"
  }
}

task igvreports {
  input {
    File ref
    File vcf
    Array[File?]+ bams
    Array[File?]+ bai
    String outDir
    String args = "--flanking 1000"
  }

  Array[File] valid_bams = select_all(bams)
  command <<<
    set -e

    mkdir -p ~{outDir}/igv
    create_report ~{vcf} --fasta ~{ref} \
        --tracks ~{vcf} ~{sep=" " valid_bams} \
        --output ~{outDir}/igv/igvreport.html ~{args}
  >>>

  output {
    File report = "~{outDir}/igv/igvreport.html"
  }
}

workflow insilicosv_basic {
  input {
    # genome simulation
    String outDir
    File? configYAML
    File? insilicosvFa
    
    # read simulation
    Int hapCoverage
    
    # --- short reads
    Boolean short

    # --- PacBio long reads
    Boolean hifi
    
    # --- ONT long reads
    Boolean ont

    # --- custom read simulation
    String? customSrPreset 
    String? customLrPreset
    
    # alignment
    File ref
    String? customLrAlnPreset
    
    # visualization
    String? igvArgs

    # resources
    Int threads
  }

  if (defined(configYAML)) {
    call insilicosv {
        input:
            yaml=select_first([configYAML]),
            outDir=outDir
    }
  }
  File syntheticFa = select_first([insilicosvFa, insilicosv.fa])
  if (short) {
    call dwgsim {
      input:
          fa = syntheticFa,
          hapCoverage = hapCoverage,
          outDir = outDir,
          customArgs = customSrPreset
    }
  
    call minimap2 as srAln {
      input:
          ref = ref,
          R1 = dwgsim.R1,
          R2 = dwgsim.R2,
          platform = "short",
          outDir = outDir,
          threads = threads
    }
  }
  
  if (hifi) {   
    call pbsim as hifiSim{
      input:
          fa = syntheticFa,
          hapCoverage = hapCoverage,
          outDir = outDir,
          platform = "hifi"
    }
  
    call minimap2 as hifiAln {
      input:
          ref = ref,
          R = hifiSim.reads,
          platform = "hifi",
          outDir = outDir,
          threads = threads 
    }
  }

  if (ont) {   
    call pbsim as ontSim {
      input:
          fa = syntheticFa,
          hapCoverage = hapCoverage,
          outDir = outDir,
          platform = "ont"
    }
  
    call minimap2 as ontAln {
      input:
          ref = ref,
          R = ontSim.reads,
          platform = "ont",
          outDir = outDir,
          threads = threads
    }
  }

  if (defined(customLrPreset)) {
    call pbsim as customSim{
      input:
          fa = syntheticFa,
          platform = "custom",
          hapCoverage = hapCoverage,
          outDir = outDir,
          customArgs = customLrPreset
    }

    call minimap2 as customAln {
      input:
          ref = ref,
          R = customSim.reads,
          platform = "custom",
          customPreset = customLrAlnPreset,
          outDir = outDir,
          threads = threads
    }
  } 

  if (defined(configYAML)) {
      call igvreports {
        input:
            ref = ref,
            vcf = select_first([insilicosv.vcf]),
            bams = [srAln.bam, hifiAln.bam, ontAln.bam, customAln.bam],
            bai = [srAln.bai, hifiAln.bai, ontAln.bai, customAln.bai],
            outDir = outDir,
            args = igvArgs
      }
  }

  output {
     File? vcf = insilicosv.vcf
     File? genome = insilicosv.fa

     File? srBAM = srAln.bam
     File? srBAI = srAln.bai

     File? hifiBAM = hifiAln.bam
     File? hifiBAI = hifiAln.bai
      
     File? ontBAM = ontAln.bam
     File? ontBAI = ontAln.bai

     File? customBAM = customAln.bam
     File? customBAI = customAln.bai

     File? igvHTML = igvreports.report
  }
}
