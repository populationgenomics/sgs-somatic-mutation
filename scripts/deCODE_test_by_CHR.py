#!/usr/bin/env python3

# Import library for analysis-runner
import click
import hail as hl
from cpg_utils.hail_batch import output_path

# Import library for deCODE job
import pandas as pd
from pandas import Series,DataFrame
import numpy as np
import math
import collections

@click.command()
@click.option('--dataset', help="data to query")
@click.option('--chrom', help="chromsome")
@click.option('--output', help="output name")
@click.option('--rerun', help='Whether to overwrite cached files', default=False)
def query(dataset, chrom, output, rerun):
    p_out = output_path(output)
    log_out = output_path(output + ".log")

    hl.init(default_reference='GRCh38')
    hl.context.warning("inited")
    hl.copy_log(log_out)

    if rerun or not hl.hadoop_exists(p_out):
        # Section 1 - Read dataset
        mt = hl.read_matrix_table(dataset)
        mt = mt.filter_rows(mt.locus.contig == chrom) 
        mt = hl.experimental.densify(mt)
        mt = hl.variant_qc(mt)        

        # Section 2 - Sample-level QC
        mt2 = hl.sample_qc(mt)
       
        # imputed sex equals to XX or XY
        # Call rate >= 0.99
        # mean coverage >= 20X
        filter_conditions = ((mt2.meta.sex_karyotype == "XX") | (mt2.meta.sex_karyotype == "XY")) & \
                            (mt2['sample_qc'].call_rate >= 0.99) & \
                            (mt2['sample_qc'].dp_stats.mean >= 20)         
        mt3 = mt2.filter_cols(filter_conditions, keep=True)
 
        # Exclude related samples
        mt4 = mt3.filter_cols(mt3['meta'].related, keep=False)
        
        # Skip ancestry check

        # Section 3 - Variant-level QC
        
        mt4 = hl.variant_qc(mt4) 

        # Apply AS_VQSR filters
        # Restricted to bi-alleleic variants
        # Exclude variants with inbreeding Coefficient < -0.3
        filter_conditions = (hl.is_missing(mt4['allele_type'])) | \
                            ((hl.is_defined(mt4['allele_type'])) & (mt4['allele_type']=='snv') & \
                             (mt4['AS_VQSLOD'] < mt4['filtering_model'].snv_cutoff.min_score)) | \
                            ((hl.is_defined(mt4['allele_type'])) & (mt4['allele_type']=='ins') & \
                             (mt4['AS_VQSLOD'] < mt4['filtering_model'].indel_cutoff.min_score)) | \
                            ((hl.len(mt4.alleles) < 2) | (hl.len(mt4.alleles) > 2)) | \
                            ((hl.len(mt4.alleles)==2) & (mt4.n_unsplit_alleles != 2)) | \
                            (hl.is_missing(mt4['InbreedingCoeff'])) | \
                            ((hl.is_defined(mt4['InbreedingCoeff'])) & (mt4['InbreedingCoeff'] < -0.3))
        mt5 = mt4.filter_rows(filter_conditions, keep=False)
         
        # Restricted to high quality variants (GQ>=20, DP>=10)
        filter_condition = (mt5.DP >= 10) & (mt5.GQ >= 20)
        mt6 = hl.variant_qc(mt5.filter_entries(filter_condition, keep=True))

        # Section 4 - deCODE specific filter
       
        # Exclude variants with call rate < 0.99 (not in deCODE paper)
        mt7 = mt6.filter_rows(mt6.variant_qc.call_rate >= 0.99)

        # Identify singleton mutations without filters 
        mt8 = mt7.filter_rows(mt7.variant_qc.n_non_ref == 1)

        # Identify singleton mutations with filters
          # sequence depth >= 16 & GQ >= 90
          # >=3 indepedent reads containing a variant allele required
        filter_condition = (mt8.DP >= 16) & (mt8.GQ >= 90) & (mt8.AD[1] >= 3)
        mt9 = hl.variant_qc(mt8.filter_entries(filter_condition, keep=True))
        mt10 = mt9.filter_rows(mt9.variant_qc.n_non_ref == 1) 
        
        # gnomAD allele frequency 
        ref_ht = ('gs://cpg-reference/seqr/v0-1/combined_reference_data_grch38-2.0.4.ht')
        ref_ht = hl.read_table(ref_ht)
      
        # Annotate variants with CADD scores, gnomAD etc.
        mt11 = mt10.annotate_rows(
                    cadd=ref_ht[mt10.row_key].cadd,
                    gnomad_genomes=ref_ht[mt10.row_key].gnomad_genomes,
                    gnomad_genome_coverage=ref_ht[mt10.row_key].gnomad_genome_coverage
               )
        
        del ref_ht

        # Apply gnomAD AF filter
        AF_cutoff = 1/(11262*2)
        mt12 = mt11.filter_rows(mt11.gnomad_genomes.AF_POPMAX_OR_GLOBAL <= AF_cutoff)
         
        # Exclude variants in simple repeat regions 
        # simple repeat regions - combining the entire Simple Tandem Repeats by TRF track in UCSC hg38 with all homopolymer regions in hg38 of length 6bp or more
 
        # Read the (Combined) Simple Repeat Regions 
        simple_repeat_regions =('gs://cpg-sgs-somatic-mtn-test-upload/Simple_Repeat_Regions_GRCh38_Excluded_Unmapped_Regions.bed')
        interval_table = hl.import_bed(simple_repeat_regions, reference_genome='GRCh38')
         
        mt13 = hl.variant_qc(mt12.filter_rows(hl.is_defined(interval_table[mt12.locus]), keep=False))
        
        # Section 5 - Write out to pVCF file (much smaller than mt)
        output = mt13.select_rows(mt13.rsid, mt13.qual)
        output = output.select_entries(output.GT, output.DP, output.AD, output.GQ)

        metadata = {'filter': {'InbreedingCoeff':{'Description':''},
            'VQSR':{'Description':''},
            'AC0':{'Description':''},
            'PASS':{'Description':''}},
            'format': {'AD': {'Description':'AD',
                'Number':'R',
                'Type':'Integer'}}
        }

        hl.export_vcf(output, p_out, metadata=metadata)
        hl.copy_log(log_out)
        
    hl.context.warning("done")  
    hl.copy_log(log_out)  

if __name__ == "__main__":
    query()
