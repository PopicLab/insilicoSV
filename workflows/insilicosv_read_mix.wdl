version 1.0
import "insilicosv_tree_level_workflow.wdl" as insilico_tree
import "insilicosv_basic_workflow.wdl" as basic

task MergeBams {
  input {
    Array[File]+ bams
    String outDir
    String sampleName
  }

  command <<<
    set -e
    samtools merge -f ~{outDir}/~{sampleName}.bam ~{sep=" " bams}
    samtools index ~{outDir}/~{sampleName}.bam
  >>>

  output {
    File bam = "~{outDir}/~{sampleName}.bam"
    File bai = "~{outDir}/~{sampleName}.bam.bai"
  }
}

workflow SimulateReadMix {
    input {
        String outDir # top-level output directory
        String originalRef # starting reference genome
        Array[Array[Int]] genomes # list of config indices (i.e. variant sets) in each genome
        Array[File] treeRefs
        Array[Float] prevalence # cellular prevalence of each genome
        Int totalCoverage # combined read depth
        Boolean? short
        Boolean? hifi
        Boolean? ont
        Int? threads
    }

    scatter (i in range(length(genomes))) {
        call insilico_tree.PathNames {
            input:
                path = genomes[i]
        }
        call basic.insilicosv_basic as ReadSim {
            input:
                outDir = outDir + "/" + PathNames.pathName,
                insilicosvFa = outDir + "/" + PathNames.pathName + "/" + PathNames.pathName + ".fa",
                hapCoverage = ceil(totalCoverage * prevalence[i] / 2),
                ref = originalRef,
                short = select_first([short, false]),
                hifi = select_first([hifi, false]),
                ont = select_first([ont, false]),
                threads = select_first([threads, 1])
        }
    }
    if (select_first([short, false])) {
        call MergeBams as mergeShort {
            input:
                bams  = select_all(ReadSim.srBAM),
                outDir = outDir,
                sampleName = "mix_short"
        }
    }
    if (select_first([hifi, false])) {
        call MergeBams as mergeHifi {
            input:
                bams  = select_all(ReadSim.hifiBAM),
                outDir = outDir,
                sampleName = "mix_hifi"
        }
    }
    if (select_first([ont, false])) {
        call MergeBams as mergeOnt {
            input:
                bams  = select_all(ReadSim.ontBAM),
                outDir = outDir,
                sampleName = "mix_ont"
        }
    }
    output {
        File? srBAM = mergeShort.bam
        File? srBAI = mergeShort.bai
        File? hifiBAM = mergeHifi.bam
        File? hifiBAI = mergeHifi.bai
        File? ontBAM = mergeOnt.bam
        File? ontBAI = mergeOnt.bai
    }
}
