nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  ===================================
  BJ-DNA-QC   P I P E L I N E  v${workflow.manifest.version}
  ===================================
  glob fastq files    : ${ params.reads }
  csv fastq files     : ${ params.input_csv }
  publish_dir         : ${ params.publish_dir }
  timestamp           : ${ params.timestamp }
  genome              : ${ params.genome }
  ${ params.skip_subsampling ? "skip_subsampling    : ${ params.skip_subsampling }" : "subsampling n_reads : ${ params.n_reads }"}
  skip_kraken         : ${ params.skip_kraken }
  skip_fastqc         : ${ params.skip_fastqc }
  skip_qualimap       : ${ params.skip_qualimap }
  skip_ginkgo         : ${ params.skip_ginkgo }
  skip_mapd           : ${ params.skip_mapd }
  ${ params.run_deepvariant ? "run_deepvariant     : ${ params.run_deepvariant }" : ""}
  \n
  """

}

def helpMessage() {

  def yellow = "\033[0;33m"
  def blue = "\033[0;34m"
  def white = "\033[0m"
  def red = "\033[0;31m"

  log.info """\
${blue}
    -------------------------------
    bj-dna-qc pipeline v${workflow.manifest.version}
    -------------------------------
    Evaluate the quality of single-cell libraries

    Usage:
        nextflow run main.nf --input_csv <input.csv> --publish_dir <output_dir> [options]

    Script Options: see nextflow.config

${red}
        [required]
        --input_csv         FILE    Path to input csv file

        --genome            STR     Reference genome to use. Available options - GRCh38, GRCm39
                                    DEFAULT: ${params.genome}
        
        --publish_dir       DIR     Path to run output directory

${yellow}
        [optional]
        
        --genomes_base      STR     Path to the genomes directory
                                    DEFAULT: ${params.genomes_base}
                                    
        --timestamp         STR     User can specify timestamp otherwise uses runtime generated timestamp 

        --n_reads           VAL     Number of reads to sample for analysis eg. 2.5M == 5M paired reads
                                    DEFAULT: ${params.n_reads}

        --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                    DEFAULT: ${params.read_length}

        --min_reads         VAL     Minimum number of reads required for analysis. Samples with fewer reads will be flagged.
                                    DEFAULT: ${params.min_reads}
        
        --run_sentieon      BOOL    Will default to all Sentieon modules if set to true, including BWA, Deduplication, Driver Metrics, and DNAscope
                                    DEFAULT: ${params.run_sentieon}
        
        --run_deepvariant   BOOL    Run Google DeepVariant to call germline variants
                                    DEFAULT: ${params.run_deepvariant}

        --skip_subsampling  BOOL    Skip subsampling of input reads
                                    DEFAULT: ${params.skip_subsampling}
                                    
        --skip_kraken       BOOL    Skip KRAKEN2 module
                                    DEFAULT: ${params.skip_kraken}

        --skip_fastqc       BOOL    Skip FastQC module
                                    DEFAULT: ${params.skip_fastqc}
                                    
        --skip_qualimap     BOOL    Skip Qualimap module
                                    DEFAULT: ${params.skip_qualimap}

        --skip_ginkgo       BOOL    Skip CNV Ginkgo and QC_plot modules
                                    DEFAULT: ${params.skip_ginkgo}
                                    
        --skip_mapd         BOOL    Skip MAPD module. MAPD is a measurement of the bin-to-bin variation in read coverage that is robust to the presence of CNVs, and is an indicator of the evenness of whole genome amplification (WGA)
                                    DEFAULT: ${params.skip_mapd}

        --skip_sigprofile   BOOL    Skip Mutational Signature Profiling
                                    DEFAULT: ${params.skip_sigprofile}

        --email_on_fail     STR     Email to receive upon failure
                                    
        --help              BOOL    Display help message
                                    
${white}
    """.stripIndent()
}

workflow{
  printHeader()
  helpMessage()
}
