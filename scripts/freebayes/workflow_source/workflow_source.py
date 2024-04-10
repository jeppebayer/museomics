#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def freebayes_vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`VCF` using :script:`freebayes`.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    SAMPLE_LIST: list = config['sample_list']
    REFERENCE_GENOME: str = config['reference_genome_path']
    OUTPUT_DIR: str = config['output_directory_path']
    WORK_DIR: str = config['working_directory_path']
    PARTITION_SIZE: int = config['partition_size']
    PLOIDY: int = config['sample_ploidy']
    BESTN: int = config['best_n_alleles']
    ALT_FRACTION: float | int = config['min_alternate_fraction']
    ALT_COUNT: int = config['min_alternate_count']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
        # If files exists reads data directly from files
        # Loads reference genome partitioning
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
            partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
            npadding = len(str(sum(1 for line in partitions)))
        # Loads list of contigs in reference genome
        with open('reference_sequences.txt', 'r') as infile:
            contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
    else:
        # If files don't exist, generate data and write files
        sequences = parse_fasta(REFERENCE_GENOME)
        with open(f'reference_sequences.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
        # Partitions reference genome
        npadding = padding_calculator(parse_fasta=sequences, size=PARTITION_SIZE)
        partitions = partition_chrom(parse_fasta=sequences, size=PARTITION_SIZE, npad=npadding)
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
        # Creates list of contigs in reference genome
        contigs = [{'contig': contig['sequence_name']} for contig in sequences]

    top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes'
    output_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes'

    freebayes_parts = gwf.map(
        template_func=freebayes_chrom,
        inputs=partitions,
        name=name_freebayes_chrom,
        extra={'reference_genome_file': REFERENCE_GENOME,
               'bam_file_list': SAMPLE_LIST,
               'output_directory': top_dir,
               'species_name': SPECIES_NAME,
               'ploidy': PLOIDY,
               'best_n_alleles': BESTN,
               'min_alternate_fraction': ALT_FRACTION,
               'min_alternate_count': ALT_COUNT}
        )
    
    concat_freebayes = gwf.target_from_template(
            name='concatenate_freebayes',
            template=concat_vcf(
                files=collect(freebayes_parts.outputs, ['vcf'])['vcfs'],
                output_name=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{BESTN}_p{PLOIDY}_minaltfrc{ALT_FRACTION}_minaltcnt{ALT_COUNT}',
                output_directory=output_dir,
                compress=True
            )
        )

    return gwf

def freebayes_population_set_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`VCF` file for each sample in configuration.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    REFERENCE_GENOME: str = config['reference_genome_path']
    SAMPLE_LIST: list = config['sample_list']
    PARTITION_SIZE: int = config['partition_size']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
        # If files exists reads data directly from files
        # Loads reference genome partitioning
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
            partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
            npadding = len(str(sum(1 for line in partitions)))
        # Loads list of contigs in reference genome
        with open('reference_sequences.txt', 'r') as infile:
            contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
    else:
        # If files don't exist, generate data and write files
        sequences = parse_fasta(REFERENCE_GENOME)
        with open(f'reference_sequences.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
        # Partitions reference genome
        npadding = padding_calculator(parse_fasta=sequences, size=PARTITION_SIZE)
        partitions = partition_chrom(parse_fasta=sequences, size=PARTITION_SIZE, npad=npadding)
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
        # Creates list of contigs in reference genome
        contigs = [{'contig': contig['sequence_name']} for contig in sequences]

    for sample in SAMPLE_LIST:
        NAME: str = sample['sample_name']
        GROUP: str = sample['sample_group'].lower()
        BAM: str = sample['bam_file']
        WORK_DIR: str = sample['working_directory_path']
        OUTPUT_DIR: str = sample['output_directory_path']
        FREEBAYES_SETTINGS: dict = sample['freebayes_settings']

        top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes/{GROUP}/{NAME}'
        if not OUTPUT_DIR:
            output_dir = f'{top_dir}'
        else:
            output_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/{GROUP}/{NAME}'

        freebayes_parts = gwf.map(
            name=name_freebayes_chrom,
            template_func=freebayes_chrom_single_sample,
            inputs=partitions,
            extra={'reference_genome_file': REFERENCE_GENOME,
                   'bam_file': BAM,
                   'output_directory': top_dir,
                   'sample_name': NAME,
                   'ploidy': FREEBAYES_SETTINGS['sample_ploidy'],
                   'best_n_alleles': FREEBAYES_SETTINGS['best_n_alleles'],
                   'min_alternate_fraction': FREEBAYES_SETTINGS['min_alternate_fraction'],
                   'min_alternate_count': FREEBAYES_SETTINGS['min_alternate_count']}
        )

        concat_freebayes = gwf.target_from_template(
                name=f'concatenate_freebayes_{NAME.replace("-", "_")}',
                template=concat_vcf(
                    files=collect(freebayes_parts.outputs, ['vcf'])['vcfs'],
                    output_name=f'{NAME}.freebayes_n{FREEBAYES_SETTINGS['best_n_alleles']}_p{FREEBAYES_SETTINGS['sample_ploidy']}_minaltfrc{FREEBAYES_SETTINGS['min_alternate_fraction']}_minaltcnt{FREEBAYES_SETTINGS['min_alternate_count']}',
                    output_directory=output_dir,
                    compress=True
                )
        )
    
    return gwf
