class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: phyllite
baseCommand:
  - python
  - /phyllite/phyllite/phyllite.py
inputs:
  - id: dna_normal_vaf
    type: File
    inputBinding:
      position: 0
      prefix: '--dna-normal-vaf'
  - id: dna_tumor_vaf
    type: File
    inputBinding:
      position: 0
      prefix: '--dna-tumor-vaf'
  - id: rna_normal_vaf
    type: File
    inputBinding:
      position: 0
      prefix: '--rna-normal-vaf'
  - id: rna_tumor_vaf
    type: File
    inputBinding:
      position: 0
      prefix: '--rna-tumor-vaf'
  - id: min_depth
    type: int?
    inputBinding:
      position: 0
      prefix: '--min-depth'
  - id: max_dna_minor_count
    type: int?
    inputBinding:
      position: 0
      prefix: '--max-dna-minor-count'
  - id: max_dna_minor_vaf
    type: float?
    inputBinding:
      position: 0
      prefix: '--max-dna-minor-vaf'
  - id: min_rna_minor_count
    type: int?
    inputBinding:
      position: 0
      prefix: '--min-rna-minor-count'
  - id: min_rna_minor_vaf
    type: float?
    inputBinding:
      position: 0
      prefix: '--min-rna-minor-vaf'
  - id: no_vaf_header
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--no-vaf-header'
outputs:
  - id: table_output
    type: File?
    outputBinding:
      glob: output.tsv
  - id: json_output
    type: File?
    outputBinding:
      glob: output.json
label: phyllite
arguments:
  - position: 0
    prefix: '--table-output'
    valueFrom: output.tsv
  - position: 0
    prefix: '--json-output'
    valueFrom: output.json
requirements:
  - class: DockerRequirement
    dockerPull: 'estorrs/phyllite:0.0.3'
