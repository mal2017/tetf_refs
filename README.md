# REFS2: generation of references and resources for Ellison lab TETF project


## dependencies

- `snakemake` (tested with v7.21.0 on Ubuntu 22 LTS and 7.18.2 on CentOS 7)
- `mamba` (tested with version 1.0.0 on Ubuntu 22 LTS and 0.9.1 on CentOS 7). `conda` alone should work just fine.
- `singularity` (tested with v3.1.0-1 on CentOS 7) or `apptainer` (tested with v1.1.3 on Ubuntu 22 LTS)

## usage

```
snakemake --use-conda --use-singularity --cores <threads>
```

## Description

Pulls various reference sequences and annotations to produce the final references, annotations, and indices used by the Ellison lab's TETF project.