version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human exome sequencing data.
##
## Requirements/expectations :
## - Array of Human exome sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "https://raw.githubusercontent.com/talkowski-lab/gatk4-exome-analysis-pipeline/master/tasks/GermlineVariantDiscovery.wdl" as Calling
import "https://raw.githubusercontent.com/talkowski-lab/gatk4-exome-analysis-pipeline/master/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/talkowski-lab/gatk4-exome-analysis-pipeline/master/tasks/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/talkowski-lab/gatk4-exome-analysis-pipeline/master/tasks/VariantCalling.wdl" as ToGvcf
import "https://raw.githubusercontent.com/talkowski-lab/gatk4-exome-analysis-pipeline/master/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow ExomeGermlineSingleSample {

  input {
    Array[File] input_bams
    Array[File] input_bais
    Array[String] base_file_name
    Array[String] final_gvcf_base_name

    Int preemptible_tries             # default: 3
    Int agg_preemptible_tries         # defult: 3

    Int break_bands_at_multiples_of   # default: 1000000
    File calling_interval_list        # gs://broad-references/hg38/v0/exome_calling_regions.v1.interval_list
    File evaluation_interval_list     # gs://broad-references/hg38/v0/exome_evaluation_regions.v1.interval_list
    File dbsnp_vcf                    # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    File dbsnp_vcf_index              # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
    File ref_fasta                    # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    File ref_fasta_index              # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
    File ref_dict                     # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dict
    File ref_alt                      # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
    File ref_ann                      # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
    File ref_sa                       # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa
    File ref_amb                      # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
    File ref_bwt                      # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
    File ref_pac                      # gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
    Int haplotype_scatter_count       # 50

  }

  scatter (scatter_index in range(length(input_bams))){

    call ToGvcf.VariantCalling as BamToGvcf {
      input:
        calling_interval_list = calling_interval_list,
        evaluation_interval_list = evaluation_interval_list,
        haplotype_scatter_count = haplotype_scatter_count,
        break_bands_at_multiples_of = break_bands_at_multiples_of,
        input_bam = input_bams[scatter_index],
        input_crai = input_bais[scatter_index],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        base_file_name = base_file_name[scatter_index],
        final_vcf_base_name = final_gvcf_base_name[scatter_index],
        agg_preemptible_tries = agg_preemptible_tries
    }


  }

  


  # Outputs that will be retained when execution is complete
  output {
    Array[File] gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    Array[File] gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics
    Array[File] output_vcf = BamToGvcf.output_vcf
    Array[File] output_vcf_index = BamToGvcf.output_vcf_index
  }

  meta {
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }
}
