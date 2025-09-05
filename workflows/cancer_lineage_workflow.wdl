version 1.0

import "./insilicosv_workflow.wdl" as insilicosv_workflow

struct GenomeParameter {
  File config_path
  Int coverage
  Boolean generate_reads
}

task ClonalGenomeGenerator {
input {
String config_path
}

command <<<
    # Run the Python script
    python3 <<CODE

import yaml
import os
import shutil
import json

config_list = []
folder_path = os.path.dirname(config_path)
with open(config_path) as config_yaml:
    config = yaml.safe_load(config_yaml)
reference = config['reference']
hap_coverage = config['coverage']
seen_config = []
for clone_name, dependencies in config['clones'].items():
    current_path = folder_path
    previous_vcf_path = ''
    for dependency in dependencies:
        current_path = current_path + '/dependency_' + str(dependency) + '/'
        config_name = config['config_files'][dependency].split('/')[-1]

        if os.path.exists(current_path):
            previous_vcf_path = current_path
            continue

        os.makedirs(current_path)
        shutil.copy(config['config_files'][dependency], current_path)

        # Add the previous config path to the genomes to generate, if not already queued
        if previous_vcf_path:
            if previous_vcf_path + config_name not in seen_config:
                config_list.append((previous_vcf_path + config_name, 0, False))
                seen_config.append(previous_vcf_path + config_name)

        # use the previous genome as the reference
        with open(current_path + config_name, 'a') as file:
            if previous_vcf_path:
                file.write('\nreference: ' + previous_path  + 'sim.hapA.fa\n')
            else:
                file.write('\nreference: ' + reference + "\n")
            # The workflow generate a single haplotype, can be ran twice for diploid
            file.write('\nhomozygous_only: True\n')
            file.write('\nhaploid: True\n')

        previous_vcf_path = current_path

    # Terminal config, corresponds to a clone to generate
    config_params = (previous_vcf_path + config_name, ceil(config['frequency'][clone_name] * hap_coverage / 100), True)
    if previous_vcf_path + config_name not in seen_config:
        config_list.append(config_params)
    else:
        # The config was already in the list, we need to update the frequency and the boolean to indicate reads are needed.
        pos_config = seen_config.index(previous_vcf_path + config_name)
        config_list[pos_config] = config_params

    seen_config.append(previous_vcf_path + config_name)

# Print outputs for WDL to capture
print(f"WDL_OUTPUT_REFERENCE={reference}")
json_output = json.dumps(config_list)
print(f"WDL_OUTPUT_CONFIG_LIST={json_output}")



CODE
>>>
output {
    String reference = read_string(stdout())
    Array[GenomeParameter] config_list = read_json(stdout())
}

workflow cancer_lineage_workflow {
  input {
    # genome simulation
    File configYAML

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
    String? customLrAlnPreset

    # visualization
    String? igvArgs

    # resources
    Int threads
  }

  String outDir = sub(configYAML, basename(configYAML), "")

  call ClonalGenomeGenerator {
    input:
         config_path=configYAML
  }
  scatter( item in ClonalGenomeGenerator.config_list) {
    config_path=item.config_path
    coverage=item.frequency
    generate_reads=item.generate_reads

    current_short=short
    current_ont=ont
    current_hifi=hifi

    if (generate_reads) {
        current_short=false
        current_ont=false
        current_hifi=false
    }

    call insilicosv_workflow.insilicosv_workflow(
        input:
            configYAML=config_path
            hapCoverage=coverage
            short=current_short
            hifi=current_hifi
            ont=current_ont
            customSrPreset=customSrPreset
            customLrPreset=customLrPreset
            ref=ClonalGenomeGenerator.reference
            customLrAlnPreset=customLrAlnPreset
            igvArgs=igvArgs
            threads=threads
    )
    # merge the reads
    if (short) {
        command <<<
            find ~{outDir}  -type f -name "*short.minimap2.sorted.bam" |
            samtools merge -@ ~{threads} -o ~{outDir}.short.cancer.bam -
            samtools index -@ ~{threads} ~{outDir}.short.cancer.bam
        >>>
    }
    if (ont) {
        command <<<
            find ~{outDir} -type f -name "*ont.minimap2.sorted.bam" |
            samtools merge -@ ~{threads} -o ~{outDir}.ont.cancer.bam -
            samtools index -@ ~{threads} ~{outDir}.ont.cancer.bam
        >>>
    }
    if (hifi) {
        command <<<
            find ~{outDir} -type f -name "*hifi.minimap2.sorted.bam" |
            samtools merge -@ ~{threads} -o ~{outDir}.hifi.cancer.bam -
            samtools index -@ ~{threads} ~{outDir}.hifi.cancer.bam
        >>>
    }
  }
  output {
  }

}
