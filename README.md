# CanSAR - Druggable Disease Networks

## Summary
Code for creating protein druggable disease networks. Protein interaction networks and druggability derived from CanSAR (https://cansar.ai/) 

## How to run

### Requirements

#### Perl

Tested on Perl version >5.20 but may work on earlier versions.

Perl modules used:

* Graph::Undirected
* Parallel::ForkManager
* Getopt::Long

#### Interactome
CanSAR interactome translated into SIF format. Interactions included are:

* OncoKB drugtarget to biomarker mapped associations extracted from OncoKB (https://www.oncokb.org/actionableGenes#sections=Tx)
* Signor signalling interactions derived from the Signor database (https://signor.uniroma2.it/)
* complex protein interactions within the same complex
* direct protein interactions directly bound to each other (includes crystallographic information and ligandable binding interfaces)
* reaction post-translational modifications indicating enzyme and substrate
* transcriptional 

#### Druggability information
Proteome wide assessed best ligandability status derived from CanSAR. Best ligandability reflects the most advanced compound status. FDA approved > In clinical trials > Quality probes from Chemical Probes Portal (https://www.chemicalprobes.org/) and ProbeMiner (https://probeminer.icr.ac.uk/#/)  > Predicted ligandable based on 3D structure. 

### generate_ddn.pl
This script creates the DDN.

```
perl generate_ddn.pl --hits <file with Uniprot primary accession list> --sif 20210427_interactome.sif --drug 20210527_best_druggability_status.txt --threads 8 --out generated_ddn.1.sif 
```

### find_best_sif.pl
The creation of a DDN is not deterministic, so the best thing to do is create several and then using this script to pick the best candidate. The script will order candidate networks based on the most advanced ligandability status followed by minimum network size footprint.

```
perl find_best_sif.pl --hits <file with Uniprot primary accession list> --dir <directory containing sif files created by generate_ddn.pl> --drug 20210527_best_druggability_status.txt --sif generated_ddn.1.sif --sif generated_ddn.2.sif --sif generated_ddn.3.sif 
```

### Example files
Example versions of the interactome and druggability files are located in the example_files folder.

## Citation

canSAR: update to the cancer translational research and drug discovery knowledgebase. DOI: 10.1093/nar/gkac1004

