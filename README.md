## Complimentary code for analyses used in the publication: 
# Genomic Identification and Characterization of Saxitoxin Producing Cyanobacteria in Western Lake Erie Harmful Algal Blooms (https://pubs.acs.org/doi/10.1021/acs.est.4c10888)

### Approach:
- Analysis utilized the Great Lakes Atlas for Multi-omics Research analysis products (GLAMR - https://greatlakesomics.org)
- R Markdown files direct the workflow and are divided by topic/approach
- Snakemake file is directly called through R Markdown code

### Data Availability:
Raw sequence files associated with this study can be accessed at the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) under BioProjects PRJNA1190886 (accession numbers SRR31514638 - SRR31514807) and PRJNA1191405 (accession numbers SRR31534945 - SRR31535021). Metagenome-assembled genomes (MAGs) from these Bioprojects are listed under accession numbers SAMN45874233 - SAMN45875077 and SAMN45872214 - SAMN45872681. Metadata, summary statistics, and sequence data are available for WLE MAGs, dereplicated at 99% ANI within each individual sample, though Zenodo (10.5281/zenodo.14803595). Environmental and sample metadata are accessible at the National Centers for Environmental Information (NCEI): 10.25921/11da-3x54.

### Notebooks:

**WLE_sxt_notebook1_FindMAGs.Rmd**
  - Identify the presence of saxitoxin genes in metagenomic data
  - Refine to only include metagenome amplified genomes (>90% completeness, <10% contamination) with saxitoxin genes
    
**WLE_sxt_notebook2_ADAMAGRefine.Rmd**
  - Investigate the presence of Anabaena-Dolichospermum-Aphanizomenon (ADA)-clade MAGs, regardless of saxitoxin gene presence
    
**WLE_sxt_notebook3_sxtPresenceAbsence.Rmd**
  - Characterize the spatiotemproal distribution of saxitoxin production potential
    
**WLE_sxt_notebook4_MAGPresenceAbsence.Rmd**
  - Characterize the spatiotemproal distribution of ADA-clade MAGs
    
**WLE_sxt_notebook6_KEGGHeatmap.Rmd**
  - Run KEGGDecoder functional analysis and create a heatmap visualization of pathways for study MAGs of interest (Figure S7)
    
**WLE_sxt_notebook8_sxtA_qPCR_BasicStats.Rmd**
  - Basic inventory of *sxtA* gene concentrations from Phytoxigene CyanoDTec Toxin Assay (205-0051) 
  
**WLE_sxt_notebook10_sxt2PresenceAbsence.Rmd**
  -  Characterize the spatiotemproal distribution of saxitoxin production potential for an operon containing the *sxtX* gene
  
**WLE_sxt_notebook11_Transcriptomics.Rmd**
  - Analyses metatranscriptomic data for the presence of identified saxitoxin operon(s)
  
**WLE_sxt_notebook12_sxtA_qPCR_HurdleModel.Rmd**
  - Analyses *sxtA* gene qPCR concentrations using a hurdle model
  - Relationships between *sxtA* gene concentration and environmental parameters examined
