# Google DeepVariant Subworkflow

A self-contained Nextflow (DSL2) subworkflow that runs [Google DeepVariant](https://github.com/google/deepvariant) on one or more aligned BAM files to produce small-variant VCFs.

This directory is intentionally minimal: just `main.nf`, `nextflow.config`, and this `README.md`. It has no dependencies on any private registry, shared configuration, or organization-specific reference bucket, so it can be run by external users out of the box.

## Implementation

The workflow runs DeepVariant as its three native stages, which lets you scale each stage independently (e.g. run `call_variants` on a GPU):

```
DEEPVARIANT_MAKE_EXAMPLES_ONLY  →  DEEPVARIANT_CALL_VARIANTS  →  DEEPVARIANT_POSTPROCESS
          (CPU, sharded)              (CPU or GPU)                 (CPU)
```

## Requirements

- [Nextflow](https://www.nextflow.io/) `>= 23.10.0`
- [Docker](https://www.docker.com/) or [Singularity / Apptainer](https://apptainer.org/)
- A prepared reference directory containing a FASTA (`genome.fa` by default) and its `.fai` index
- Aligned, coordinate-sorted, indexed BAM files
- (Optional) BioSkryb custom DeepVariant model + population VCFs — see [Resource Access](#resource-access)

Container images are pulled from the public Docker Hub repository [`google/deepvariant`](https://hub.docker.com/r/google/deepvariant).

## Resource Access

The BioSkryb-trained DeepVariant model (`bioskryb-af-20241102`) and matching population VCFs are hosted in a public-access BioSkryb bucket on Wasabi. To obtain credentials, email **basejumper@bioskryb.com**.

Once you have the provided access key and secret, export them and sync the resources to a local directory:

```bash
export AWS_ACCESS_KEY_ID=<provided_access_key>
export AWS_SECRET_ACCESS_KEY=<provided_secret_key>
export AWS_DEFAULT_REGION=us-east-1

# DeepVariant model
aws s3 sync \
    s3://bioskryb-public-data/pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/DeepVariant/Illumina/bioskryb-af-20241102/model/ \
    <YOUR_LOCAL_DIRECTORY>/bioskryb-af-20241102/model/ \
    --endpoint-url https://s3.us-east-1.wasabisys.com

# Matching population VCFs
aws s3 sync \
    s3://bioskryb-public-data/pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/DeepVariant/Illumina/bioskryb-af-20241102/population_vcfs/ \
    <YOUR_LOCAL_DIRECTORY>/bioskryb-af-20241102/population_vcfs/ \
    --endpoint-url https://s3.us-east-1.wasabisys.com
```

Then point the workflow at the downloaded paths:

```bash
--deepvariant_model <YOUR_LOCAL_DIRECTORY>/bioskryb-af-20241102/model \
--population_vcfs   <YOUR_LOCAL_DIRECTORY>/bioskryb-af-20241102/population_vcfs
```

The default `--checkpoint_filename` (`sc_af.model.ckpt`) matches the file shipped inside the model directory, so no extra configuration is required.

## Layout

```
google_deepvariant_wf/
├── main.nf            # Three DEEPVARIANT_* processes + GOOGLE_DEEPVARIANT_WF + standalone entry
├── nextflow.config    # Defaults, containers, profiles (docker, singularity, gpu, test)
└── README.md
```

## Inputs

### Input CSV (`--input_csv`)

Headered CSV with one row per sample:

```
sample_name,bam,bai
HG001,/data/HG001.bam,/data/HG001.bam.bai
HG002,/data/HG002.bam
```

The `bai` column is optional; if omitted, the workflow derives the index path as `<bam>.bai`.

### Parameters

| Parameter                    | Default       | Description                                                                    |
|------------------------------|---------------|--------------------------------------------------------------------------------|
| `--input_csv`                | _required_    | CSV of samples (`sample_name,bam[,bai]`)                                       |
| `--reference`                | _required_    | Directory containing the reference FASTA and `.fai`                            |
| `--reference_fasta_name`     | `genome.fa`   | Filename of the FASTA within `--reference`                                     |
| `--regions`                  | `null`        | Optional BED file to restrict calling regions (pass this for exome/targeted)   |
| `--deepvariant_model`        | `null`            | Optional directory with a custom model checkpoint (e.g. the BioSkryb `bioskryb-af-20241102/model` — see [Resource Access](#resource-access)); leave unset to use the container's default |
| `--checkpoint_filename`      | `sc_af.model.ckpt`| Filename of the checkpoint inside `--deepvariant_model` (default matches the BioSkryb model) |
| `--population_vcfs`          | `null`            | Optional directory of per-chromosome `*.vcf.gz` used as `--population_vcfs` (BioSkryb `bioskryb-af-20241102/population_vcfs`) |
| `--deepvariant_version`      | `1.8.0`       | Tag of the public `google/deepvariant` image                                   |
| `--call_variants_batch_size` | `1024`        | `call_variants --batch_size`                                                   |
| `--outdir`                   | `./results`   | Output directory                                                               |
| `--max_cpus`                 | `32`          | CPU ceiling for `MAKE_EXAMPLES` and `POSTPROCESS`                              |
| `--max_memory`               | `120.GB`      | Memory ceiling for `MAKE_EXAMPLES` and `POSTPROCESS`                           |
| `--max_time`                 | `24.h`        | Per-task time ceiling                                                          |

> **Note on `--deepvariant_model`**: the default `--checkpoint_filename` is `sc_af.model.ckpt`, matching the BioSkryb `bioskryb-af-20241102` model. If you bring your own model, set `--checkpoint_filename` to the name of your `.ckpt` file inside that directory.

## Usage

### Standalone

```bash
nextflow run ./main.nf \
    -profile docker \
    --input_csv samples.csv \
    --reference /data/reference/GRCh38 \
    --outdir results
```

GPU-accelerated `call_variants` (uses `google/deepvariant:<version>-gpu`):

```bash
nextflow run ./main.nf \
    -profile docker,gpu \
    --input_csv samples.csv \
    --reference /data/reference/GRCh38
```

Singularity / Apptainer:

```bash
nextflow run ./main.nf \
    -profile singularity \
    --input_csv samples.csv \
    --reference /data/reference/GRCh38
```

Exome / targeted (restrict calling to a BED):

```bash
nextflow run ./main.nf \
    -profile docker \
    --input_csv samples.csv \
    --reference /data/reference/GRCh38 \
    --regions /data/panels/exome.bed
```

With the BioSkryb custom model + population VCFs (after following [Resource Access](#resource-access)):

```bash
nextflow run ./main.nf \
    -profile docker \
    --input_csv samples.csv \
    --reference /data/reference/GRCh38 \
    --deepvariant_model /data/bioskryb-af-20241102/model \
    --population_vcfs   /data/bioskryb-af-20241102/population_vcfs
```

### As a subworkflow in another pipeline

```groovy
include { GOOGLE_DEEPVARIANT_WF } from './shared/subworkflow/google_deepvariant_wf/main.nf'

workflow {
    ch_bam               = Channel.of( tuple('HG001', file('HG001.bam'), file('HG001.bam.bai')) )
    ch_reference         = file('/data/reference/GRCh38')
    ch_deepvariant_model = []         // or file('/data/custom_model')
    ch_population_vcfs   = []         // or file('/data/population_vcfs')
    ch_regions           = []         // or file('regions.bed')

    GOOGLE_DEEPVARIANT_WF(
        ch_bam,
        ch_deepvariant_model,
        ch_population_vcfs,
        ch_reference,
        ch_regions
    )

    GOOGLE_DEEPVARIANT_WF.out.vcf.view()
}
```

To reuse the process defaults (containers, CPUs, memory) from the parent pipeline:

```groovy
includeConfig './shared/subworkflow/google_deepvariant_wf/nextflow.config'
```

## Outputs

Emitted from `GOOGLE_DEEPVARIANT_WF`:

- `vcf` — tuple of `(sample_name, [*.vcf.gz, *.vcf.gz.tbi])`
- `version` — `deepvariant_version.yml` capturing the DeepVariant image tag used

The VCF and its index are named `<sample_name>_deepvariant.vcf.gz{,.tbi}`.

## Notes

- `DEEPVARIANT_MAKE_EXAMPLES_ONLY` shards across `task.cpus` using GNU `parallel`; more CPUs ⇒ more shards ⇒ faster `make_examples`.
- `DEEPVARIANT_CALL_VARIANTS` is the stage that benefits most from a GPU. Enable the `gpu` profile to switch to the `-gpu` image and pass `--gpus all`.
- The `--reference` directory is staged as a single path so the FASTA, `.fai`, and any companion files are co-located at runtime.

## License

DeepVariant is distributed under the BSD-3-Clause License by Google. See the [DeepVariant repository](https://github.com/google/deepvariant) for details.
