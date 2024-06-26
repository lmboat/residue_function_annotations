# Residue Site Annotations

## About
Counts of how many proteins had UniProt annotations for active sites, binding sites, catalytic activity, disulfide bonds and redox potentials were calculated. Further parsing of UniProt active and binding site annotations were extracted to obtain specific residues and amino acid numbers. Positions of binding and active sites that were not specified residues were discarded. Exact amino acid positions of UniProt active and binding sites were cross-referenced with residue identifiers. 

## Set Up
1. ```
   pip3 install biopython
   ```
2. Download input file <2401_uniprot.fasta> from UniProt for all canonical human protein sequences.
   https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28Human%29+AND+%28model_organism%3A9606%29+AND+%28reviewed%3Atrue%29
4. Input files <aa_ids/2401_uniprot_residueid.csv> can be calculated using get_aaids_from_fasta.py
5. Input files <compiled_identifiers.csv> are required and can be obtained from the MS-CpDAA Analysis Suite.

## Usage
1. Create main amino acid identifiers for the entire proteome for a specific residue
2. Default input file is 2401_uniprot.fasta
   
   ```
   python3 get_aaids_from_fasta.py -i <input.fasta> -o <output.csv>
   ```

   ```
   python3 get_aaids_from_fasta.py
   ```
   
4. Create main active site and binidng site identifiers from UniProtKB annotations
   ```
   active_binding_site_function_annotation.ipynb
   ```
6. Merge custom input data to identify annotated active site and binding sites within custom dataset
   ```
   disulfide_redox_site_function_annotation.ipynb
   ```
8. Create main disulfide bond and redox active identifiers from UniProtKB annotations
9. Merge custom input data to identify annotated disulfide bond and redox active sites within custom dataset


## Citation
1. Boatner LM, Palafox MF, Schweppe DK, Backus KM. CysDB: a human cysteine database based on experimental quantitative chemoproteomics. Cell Chem Biol. 2023 Jun 15;30(6):683-698.e3. doi: 10.1016/j.chembiol.2023.04.004. Epub 2023 Apr 28. PMID: 37119813; PMCID: PMC10510411.
2. Yan T, Boatner LM, Cui L, Tontonoz PJ, Backus KM. Defining the Cell Surface Cysteinome Using Two-Step Enrichment Proteomics. JACS Au. 2023 Dec 13;3(12):3506-3523. doi: 10.1021/jacsau.3c00707. PMID: 38155636; PMCID: PMC10751780.
