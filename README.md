# IB_Bioinf_tools
## Tools to work with nucleic acids, protein sequences, fastq reads.

Every bioinfomatician at some point of his/her career deals with nucleic acid sequences, protein sequences. Nowadays it's also hard to imagine bioinformatics without New Generation Sequecning methods and data, that needs to be analyzed with bioinformatician hands. For those reasons as well as for practical experience and better theorethical understanding of bioinformatics I create this repo.
`bioinf_tools.py` is an open-source program that facilitates working with bioinformatics data.

## Usage and options
The programm is based on three main functions:

**1.** `run_protein_tools` - facilitates work with protein sequences, has 5 procedures:
  
- Search for conserved amino acids residues in protein sequence: **search_for_motifs**
- Search for alternative frames in a protein sequences: **search_for_alt_frames**
- Convert protein sequences to RNA or DNA sequences: **convert_to_nucl_acids**
- Reverse the protein sequences from one-letter to three-letter format and vice-versa: **three_one_letter_code**
- Define molecular weight of the protein sequences: **define_molecular_weight**

Function works with *one-letter amino acid sequences*,  a name of procedure and a relevant argument. If you have three-letter amino acids sequences please convert it with `three_one_letter_code` procedure before using any other procedures on them.

To start with the program run the following command:

`run_protein_tools(sequences, procedure: str ="procedure", ...)`

Where:
- sequences: positional argument, a *str | list[str] | tuple[str]* of protein sequences
- procedure: keyword argument, a type of procedure to use that is inputed in *string* type
- ... - an additional keyword arguments. Those are:
  For "search_for_motif" procedure:
- motif (str): desired motif to check presense in every given sequence. Example: motif = "GA"
- overlapping (bool): count (True) or skip (False) overlapping matches. Example: overlapping = False (Optional)\
  For "search_for_alt_frames" procedure:
- alt_start_aa (str): the name of an amino acid that is encoded by alternative start codon. Example: alt_start_aa = "I" (Optional)\
  For "convert_to_nucl_acids" procedure:
- nucl_acids (str): the nucleic acid to convert to. Example: nucl_acids = "both" | nucl_acids = "RNA"

**2.** `run_dna_rna_tools` - facilitates work with nucleic acids sequences, has 4 procedures:
- Transcribe nucleic acid sequence: **transcribe**
- Convert nucleic acid sequence into its reverse counterpart: **reverse**
- Convert nucleic acid sequence into its complement counterpart: **complement**
- Convert nucleic acid sequence into its reverse-complement counterpart: **reverse_complement**

Function works with nucleic acids sequences (DNA/RNA) and a name of procedure.

To start with the program run the following command:

`run_dna_rna_tools(sequences, procedure = "procedure")`

Where:
- sequences *(str | list[str] | tuple[str])*: positional argument, a  of nucleic acids sequences
- procedure *(str)*: keyword argument, a name of procedure to use that is inputed in *string* type

**3.** `run_fastq_filter` - facilitates filtering of sequncing reads by several parameters:
- GC content: **gc_bounds**
- read length: **length_bounds** 
- mean nucleotide quality: **quality_threshold**

Function works with default arguments. 
To start with the program run the following command:

`run_fastq_filter(seqs, gc_bounds, length_bounds, quality_threshold, verbose=True)`

Where:
- seqs (dict[str, tuple[str] | list[str]]): fastq reads to be filtered. Default value: *an example dictionary*
- gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds. Default value: *(0:100)*
- length_bounds (int | tuple[int] | list[int]): read length thresholds. Default value: *(1,2**32)*
- quality_thresholds (int | float): read Phred-33 scaled quality thresholds. Default value: *0*
- verbose (bool): add detailed statistics for each read. Defaukt value: *True*

Please provide GC content bounds in percentages.

For *gc_bounds* and *length_bounds* if one number is provided bounds from 0 to number are considered.

Please provide threshold quality in Phred-33 scale.

## Examples
```python
### run_protein_tools


## three_one_letter_code
run_protein_tools(['met-Asn-Tyr', 'Ile-Ala-Ala'], procedure='three_one_letter_code')  # ['mNY', 'IAA']
run_protein_tools(['mNY','IAA'], procedure='three_one_letter_code')  # ['met-Asn-Tyr', 'Ile-Ala-Ala']


## define_molecular_weight
run_protein_tools(['MNY','IAA'], procedure='define_molecular_weight')  # {'MNY': 426.52, 'IAA': 273.35}


## check_for_motifs
run_protein_tools(['mNY','IAA'], procedure='search_for_motifs', motif='NY')
#Sequence: mNY
#Motif: NY
#Motif is present in protein sequence starting at positions: 1

#Sequence: IAA
#Motif: NY
#Motif is not present in protein sequence

#{'mNY': [1], 'IAA': []}


## search_for_alt_frames
run_protein_tools(['mNYQTMSPYYDMId'], procedure='search_for_alt_frames')  # {'mNYQTMSPYYDMId': ['MSPYYDMId']}
run_protein_tools(['mNYTQTSP'], procedure='search_for_alt_frames', alt_start_aa='T')  # {'mNYTQTSP': ['TQTSP']}


## convert_to_nucl_acids
run_protein_tools('MNY', procedure='convert_to_nucl_acids', nucl_acids = 'RNA')  # {'RNA': ['AUGAACUAU']}
run_protein_tools('MNY', procedure='convert_to_nucl_acids', nucl_acids = 'DNA')  # {'DNA': ['TACTTGATA']}
run_protein_tools('MNY', procedure='convert_to_nucl_acids', nucl_acids = 'both') # {'RNA': ['AUGAACUAU'], 'DNA': ['TACTTGATA']}

### run_dna_rna_tools
run_dna_rna_tools(("ATGC", "TGCA"), procedure='transcribe') # {'ATGC': 'AUGC', 'TGCA': 'UGCA'}
run_dna_rna_tools("ATGC", procedure='reverse') # 'CGTA'
run_dna_rna_tools("ATGC", procedure='complement') # 'TACG'
run_dna_rna_tools("ATGC", procedure='reverse_complement') # 'GCAT'

### run_fastq_filter
run_fastq_filter(gc_bounds=50)
#Read: ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA
#GC Content: 38.2022
#Read Length: 89
#Mean Nucleotide Quality: 36.1011

#Read: GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC
#GC Content: 49.4253
#Read Length: 87
#Mean Nucleotide Quality: 33.2989

#{'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
# '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')}
```
