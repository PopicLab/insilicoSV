version 1.0

task GetPathsAtLevel {
    input {
        Array[Array[Int]] genomes
        Int treeLevel
    }
    File genomesJson = write_json(object{genomes: genomes})

    command <<<
    python3 <<CODE
    import json
    with open("~{genomesJson}") as f: genomes = json.load(f)["genomes"]
    paths = set()
    for path in genomes:
        if len(path) > ~{treeLevel}:
            paths.add(tuple(path[:~{treeLevel}+1]))
    paths = [list(p) for p in sorted(paths)]
    with open("paths.json", "w") as f: json.dump({"paths": paths}, f)
    CODE
    >>>

    output {
        Array[Array[Int]] paths = read_json("paths.json")["paths"]
    }
}

task PathNames {
    input {
        Array[Int] path
    }
    File pathJson = write_json(object{path: path})

    command <<<
    python3 <<CODE
    import json
    with open("~{pathJson}") as f: path = json.load(f)["path"]
    with open("name.txt", "w") as f: f.write("_".join([str(p) for p in path]))
    with open("parent.txt", "w") as f: f.write("_".join([str(p) for p in path[:-1]]) if len(path) > 1 else "")
    CODE
    >>>

    output {
        String pathName = read_string("name.txt")
        String parentPathName = read_string("parent.txt")
    }
}

task insilicosv {
  input {
    File yaml
    String outDir
    String outPrefix
    String ref
    Int treeLevel
  }
  String yamlBasename = basename(yaml)

  command <<<
    set -e
    mkdir -p ~{outDir}
    cp ~{yaml} ~{outDir}/~{yamlBasename}
    echo -e "\nreference: ~{ref}" >> ~{outDir}/~{yamlBasename}
    if (( ~{treeLevel} > 0 )); then
      echo -e "\nhomozygous_only: True" >> ~{outDir}/~{yamlBasename}
    fi
    insilicosv -c ~{outDir}/~{yamlBasename}
    if (( ~{treeLevel} == 0 )); then
      cat ~{outDir}/sim.hapA.fa ~{outDir}/sim.hapB.fa > ~{outDir}/sim.fa
    else
      mv ~{outDir}/sim.hapA.fa ~{outDir}/sim.fa
      rm ~{outDir}/sim.hapB.fa
    fi
    mv  ~{outDir}/sim.fa ~{outDir}/~{outPrefix}.fa
    mv  ~{outDir}/sim.vcf ~{outDir}/~{outPrefix}.vcf
    samtools faidx  ~{outDir}/~{outPrefix}.fa
  >>>

  output {
    File vcf = "~{outDir}/~{outPrefix}.vcf"
    File fa = "~{outDir}/~{outPrefix}.fa"
  }
}

workflow SimulateLineageTreeLevel {
    input {
        String outDir # top-level output directory
        String originalRef # starting reference genome
        Array[File] insilicoConfigs # list of insilicoSV YAML files for each variant set/time point
        Array[Array[Int]] genomes # list of config indices (i.e. variant sets) in each genome
        Int treeLevel # lineage tree level ( 0=root/normal)
        Array[File]? previousRef # genome set generated in the previous level
    }

    call GetPathsAtLevel {
        input:
            genomes = genomes,
            treeLevel = treeLevel
    }

    scatter (path in GetPathsAtLevel.paths) {
        call PathNames {
            input:
                path = path
        }
        String ref = if treeLevel > 0 then outDir + "/" + PathNames.parentPathName + "/" + PathNames.parentPathName + ".fa" else originalRef
        Int configIndex = path[treeLevel]
        call insilicosv {
            input:
                ref = ref,
                yaml = insilicoConfigs[configIndex],
                outDir = outDir + "/" + PathNames.pathName,
                outPrefix = PathNames.pathName,
                treeLevel = treeLevel
        }
    }
    output {
        Array[File] ref = insilicosv.fa
        Array[File] vcf = insilicosv.vcf
    }
}



