
# BaseJumper BJ-DNA-QC

The BioSkryb BJ-DNA-QC pipeline evaluates the quality of the single-cell library and it provides several qc metrics to assess the quality of the sequencing reads.

One way that users can ensure that a single-cell library is uniformly amplified with low allelic dropouts, is by first sequencing using “low-pass” or low throughput sequencing of around 2M reads per sample. Data from the “low-pass” can used to estimate the genome coverage if the single-cell libraries were to be used for high-depth sequencing. Users can then select only quality libraries for high-depth sequencing.

The BJ-DNA-QC pipeline uses low-pass sequencing data and generates several QC metrics that help assess whether the single-cell libraries are ready for high-depth sequencing.

# Pipeline Overview
Following are the steps and tools that pipeline uses to perform the analyses:

- Subsample the reads to 2 million using SEQTK SAMPLE to compare metrics across samples. Bypass with `--skip_subsampling` to use all reads instead
- Evaluate sequencing quality control using FASTP and trim/clip reads
- *Map reads to reference genome using BWA MEM
- *Remove duplicate reads using PICARD MARKDUPLICATES
- *Evaluate metrics using PICARD METRICS which includes Alignment, GC Bias, Insert Size, and Coverage metrics
- Evaluate the BAM quality control using QUALIMAP BAMQC
- Evaluate the library complexity using PRESEQ BAM2MR and PRESEQ GC EXTRAP
- Evaluate the CNV using a custom Ginkgo impelmentation
- Evaluate taxonomic classification with Kraken
- Aggregate the metrics across biosamples and tools to create overall pipeline statistics summary using MULTIQC

*BWA and PICARD run by default. SENTIEON optimized equivalents can be enabled for these steps by passing `--run_sentieon`: SENTIEON BWA MEM, LOCUSCOLLECTOR, DEDUP, and DRIVER METRICS

# Running Locally

Following are instructions for running BJ-DNA-QC in a local Ubuntu server

## Install Java

```
sudo apt-get install default-jdk

java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```

## Download pipeline
```
git clone https://github.com/OpenOmics/bj-dna-qc.git
cd bj-dna-qc

```

## Sentieon License Setup

The Sentieon license is a "localhost" license that starts a lightweight license server on the localhost. This type of license is very easy to use and get started with. However, because it can be used anywhere, we restrict this license to short-term testing/evaluation only. To use this type of license, you need to set the environment variable SENTIEON_LICENSE to point to the license file on the compute nodes:
```
export SENTIEON_LICENSE=</path/to/sentieon_eval.lic>
```
The license file should be saved at the base directory of the pipeline eg: `bj-dna-qc/sentieon_eval.lic`
All users will need to  [ submit helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156) to get an evaluation/full pass-through BioSkryb's Sentieon license.

## Resources Required

For running the pipeline, a typical dataset (less than 8 million reads) requires 4 CPU cores and 14 GB of memory. For larger datasets, you may need to increase the resources to 8 CPU cores. You can specify these resources in the command as follows:
```
--max_cpus 4 --max_memory 14.GB
```

<br>

# Pipeline Execution

All pipeline resources are publically available at `s3://bioskryb-public-data/pipeline_resources` users need not have to download this, and will be downloaded during nextflow run.

**Run Command**

** Example using test data input **

```
nextflow run main.nf \
    --input_csv $PWD/tests/data/inputs/input.csv \
    --publish_dir test \
    --max_cpus 4 \
    --max_memory 14.GB
```

**Input Options**

The input for the pipeline can be passed via a samplesheet CSV.

- **CSV Metadata Input**: The CSV file should have 3 columns: `biosampleName`, `read1` and `read2`.
The `biosampleName` column contains the name of the biosample, `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,read1,read2
DNAQC-test1-100reads,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test1-100reads_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test1-100reads_S1_L001_R2_001.fastq.gz
DNAQC-test2-1000reads,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test2-1000reads_S2_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test2-1000reads_S2_L001_R2_001.fastq.gz
```

- **Optional Groups Column**: The CSV file can include an optional `groups` column containing sample group information. This column is mandatory when using the QC_Plot and Mutational Signature profile modules.

**Optional Skip Modules**

This pipeline includes several optional modules. You can choose to include or exclude these modules by adjusting the following parameters:

- `--skip_subsampling`: Set to `true` to skip subsampling of input reads. Default: `false`
- `--skip_kraken` : Set to `true` to exclude the KRAKEN2 module. Default: `true`
- `--skip_fastqc` : Set to `true` to exclude the FastQC module. Default: `true`
- `--skip_qualimap` : Set to `true` to exclude the Qualimap module. Default: `true`
- `--skip_ginkgo`: Set to `true` to exclude the CNV Ginkgo module. Default: `true`. **Note**: When the Ginkgo module is enabled, the QC_Plot module will also be automatically executed.
- `--skip_mapd`: Set to `true` to exclude the MAPD module. Default: `false`
- `--skip_sigprofile` : Set to `false` to include the Mutational Signature profile module. Default: `true`.
        The mutational signature runs require additional metadata in the input.csv file under the "groups" column. Add the group name for samples that need to be combined by the pseudobulk process for mutational signature profiling. This field can be left blank for other samples.


**Outputs**

The pipeline saves its output files in the designated "publish_dir" directory. The different QC metrics files are stored in the "secondary_analyses/metrics/<sample_name>/" subdirectory. For details: [BJ-DNA-QC outputs](https://docs.basejumper.bioskryb.com/pipelines/secondary/bj-dna-qc/1.9.1/docs/#output-directories)

**Command Options**
```
    Usage:
        nextflow run main.nf --input_csv <input.csv> --publish_dir <output_dir> [options]

    Script Options: see nextflow.config

        [Required Arguments]
        --------------------

        --input_csv         FILE    Path to input CSV samplesheet.

        --genome            STR     Reference genome to use. Available options - GRCh38, GRCm39
                                    DEFAULT: GRCh38

        --publish_dir       DIR     Path to run output directory. If the provided output directory does not exist, it will be created automatically.

        

        [Optional Arguments]
        --------------------
        
        --genomes_base      STR     Path to the genomes directory
                                    DEFAULT: /data/OpenOmics/references/bioskryb
                                    
        --timestamp         STR     User can specify timestamp otherwise uses runtime generated timestamp 

        --n_reads           VAL     Number of reads to sample for analysis eg. 2.5M == 5M paired reads
                                    DEFAULT: 2000000

        --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                    DEFAULT: 50

        --min_reads         VAL     Minimum number of reads required for analysis. Samples with fewer reads will be flagged.
                                    DEFAULT: 1000

        --run_sentieon      BOOL    Will default to all Sentieon modules if set to true, including BWA, Deduplication, Driver Metrics, and DNAscope
                                    DEFAULT: false
        
        --run_deepvariant   BOOL    Run Google DeepVariant to call germline variants
                                    DEFAULT: false


        [Optional Skip Modules]: Set to `true` to enable skip and `false` to disable skip

        --skip_subsampling  BOOL    Skip subsampling of input reads
                                    DEFAULT: false
                                    
        --skip_kraken       BOOL    Skip KRAKEN2 module
                                    DEFAULT: true

        --skip_fastqc       BOOL    Skip FastQC module
                                    DEFAULT: true
                                    
        --skip_qualimap     BOOL    Skip Qualimap module
                                    DEFAULT: true

        --skip_ginkgo       BOOL    Skip CNV Ginkgo and QC_plot modules
                                    DEFAULT: true
                                    
        --skip_mapd         BOOL    Skip MAPD module. MAPD is a measurement of the bin-to-bin variation in read coverage that is robust to the presence of CNVs, and is an indicator of the evenness of whole genome amplification (WGA)
                                    DEFAULT: false

        --skip_sigprofile   BOOL    Skip Mutational Signature Profiling
                                    DEFAULT: true
        
        --email_on_fail     STR     Email to receive upon failure
                                    
        --help              BOOL    Display help message
```

**Example Usage on Low-pass Data**
```
nextflow run main.nf \
    --input_csv $PWD/tests/data/inputs/input.csv \
    --publish_dir test \
    --skip_kraken false \
    --skip_fastqc false \
    --skip_ginkgo false
```

**Example Usage on High-pass Data**
```
nextflow run main.nf \
    --input_csv $PWD/tests/data/inputs/input.csv \
    --publish_dir test \
    --skip_subsampling \
    --skip_kraken false \
    --skip_fastqc false \
    --skip_ginkgo false \
    --run_deepvariant
```

**Tool versions**

- `Seqtk: 1.3-r106`
- `fastp: 0.20.1`
- `Kraken2: 2.1.3`
- `FastQC: v0.11.9`
- `Sentieon: 202308.01`
- `BWA: v0.7.17-r1188`
- `Picard: 3.4.0`
- `bam-lorenz-coverage: 2.3.0 GNU`
- `QualiMap: v.2.2.2-dev`
- `Preseq: 2.0.3`
- `Ginkgo: 0.0.2`
- `bedtools: v2.28.0`


**nf-test**

The BioSkryb BJ-DNA-QC nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
sudo mv nf-test /usr/local/bin/
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("preseq_test") {
        when {
            params {
                // define parameters here. Example: 
                publish_dir = "${outputDir}/results"
                timestamp = "test"
            }
        }

        then {
            assertAll(
                
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the multiqc report HTML file
                { 
                    assert new File("${outputDir}/results_test/multiqc/multiqc_report.html").exists()
                },

                // Check for a match in the all metrics MQC text file
                { 
                    assert snapshot(path("${outputDir}/results_test/secondary_analyses/metrics/nf-preseq-pipeline_all_metrics_mqc.txt"))
                        .match("all_metrics_mqc")
                },

                // Check for a match in the selected metrics MQC text file
                { 
                    assert snapshot(path("${outputDir}/results_test/secondary_analyses/metrics/nf-preseq-pipeline_selected_metrics_mqc.txt"))
                        .match("selected_metrics_mqc")
                },

                // Verify existence of the fastp JSON file
                { 
                    assert new File("${outputDir}/results_test/primary_analyses/metrics/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/fastp/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64_no_qc_fastp.json").exists()
                },

                // Check for a match in the kraken report text file
                { 
                    assert snapshot(path("${outputDir}/results_test/primary_analyses/metrics/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/kraken2/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64_kraken2_report.txt"))
                        .match("kraken_report")
                },

                // Verify existence of the qualimap HTML report
                { 
                    assert new File("${outputDir}/results_test/primary_analyses/metrics/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/qualimap/QUALIMAP_BAMQC_WF_H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/qualimapReport.html").exists()
                },

                // Check for a match in the Ginkgo segment summary text file
                { 
                    assert snapshot(path("${outputDir}/results_test/tertiary_analyses/cnv_ginkgo/AllSample-GinkgoSegmentSummary.txt"))
                        .match("ginkgo_summary")
                }
            )
            
        }

    }
```

# Need Help?
If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References
For more information, you can refer to the following publications:

- Chung, C., Yang, X., Hevner, R. F., Kennedy, K., Vong, K. I., Liu, Y., Patel, A., Nedunuri, R., 
Barton, S. T., Noel, G., Barrows, C., Stanley, V., Mittal, S., Breuss, M. W., Schlachetzki,
J. C. M., Kingsmore, S. F., & Gleeson, J. G. (2024). Cell-type-resolved mosaicism 
reveals clonal dynamics of the human forebrain. Nature, 629(8011), 384–392. [https://doi.org/10.1038/s41586-024-07292-5](https://doi.org/10.1038/s41586-024-07292-5)

- Zhao, Y., Luquette, L. J., Veit, A. D., Wang, X., Xi, R., Viswanadham, V. V, Shao, D. D., 
Walsh, C. A., Yang, H. W., Johnson, M. D., & Park, P. J. (2024).
High-resolution 
detection of copy number alterations in single cells with HiScanner. BioRxiv, 
2024.04.26.587806. [https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full](https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full)

- Zawistowski, J. S., Salas-González, I., Morozova, T. V, Blackinton, J. G., Tate, T., Arvapalli, 
D., Velivela, S., Harton, G. L., Marks, J. R., Hwang, E. S., Weigman, V. J., & West, J. A. 
A. (n.d.). Unifying genomics and transcriptomics in single cells with ResolveOME amplification
chemistry to illuminate oncogenic and drug resistance mechanisms. [https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full](https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full)

NOTE: Several studies have utilized BaseJumper pipelines as part of the standard 
quality control processes implemented through ResolveServices<sup>SM</sup>. While these 
pipelines may not be explicitly cited, they are integral to the methodologies 
described.

