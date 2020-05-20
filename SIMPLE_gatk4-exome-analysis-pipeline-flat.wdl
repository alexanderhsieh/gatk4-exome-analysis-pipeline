version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human exome sequencing data.
##
## Requirements/expectations :
## - Human exome sequencing data in unmapped BAM (uBAM) format
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

import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/UnmappedBamToAlignedBam.wdl" as ToBam
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/AggregatedBamQC.wdl" as AggregatedQC
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/GermlineVariantDiscovery.wdl" as Calling
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/BamToCram.wdl" as ToCram
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/tasks/VariantCalling.wdl" as ToGvcf
import "https://raw.githubusercontent.com/alexanderhsieh/gatk4-exome-analysis-pipeline/master/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow ExomeGermlineSingleSample {

  input {
    PapiSettings papi_settings

    String base_file_name
    String final_gvcf_base_name

    File input_cram
    File input_crai

    String sample_name
    String unmapped_bam_suffix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File ref_alt
    File ref_sa
    File ref_amb
    File ref_bwt
    File ref_ann
    File ref_pac

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    Int haplotype_scatter_count
    Int break_bands_at_multiples_of

    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices

    File dbsnp_vcf
    File dbsnp_vcf_index

    File calling_interval_list
    File evaluation_interval_list
    File target_interval_list
    File bait_interval_list

    Boolean provide_bam_output = false
    File? haplotype_database_file
  }

  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      calling_interval_list = calling_interval_list,
      evaluation_interval_list = evaluation_interval_list,
      haplotype_scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      input_bam = input_cram,
      input_crai = input_crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      base_file_name = base_file_name,
      final_vcf_base_name = final_gvcf_base_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries
  }


  # Outputs that will be retained when execution is complete
  output {


    File gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics


    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
}
