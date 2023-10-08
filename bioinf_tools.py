from src import dna_rna_tools
from src import fastq_tools
from src import protein_tools


def run_dna_rna_tools(
    sequences: str or list[str] or tuple[str], procedure: str
) -> str or dict:
    """
    Process nucleic acid sequence by one of the developed tools.\n
    Run one procedure at a time:
    - Transcribe nucleic acid sequence.
    - Convert nucleic acid sequence into its reverse counterpart.
    - Convert nucleic acid sequence into its complement counterpart.
    - Convert nucleic acid sequence into its reverse-complement counterpart.

    All functions are letter case sensitive.\n
    If only one sequence provided - *sequences* can be string.\n
    If more - please provide *sequences* as list or tuple.\n
    If more information needed please see README or desired docstring.

    Arguments:
    - sequences (str or list[str] or tuple[str]): sequences to process.
    - procedure (str]: desired procedure:
        - "transcribe"
        - "reverse"
        - "complement"
        - "reverse_complement"

    Return:
    - processed_sequence (dict or str): Sequences, processed with desired tool.
    """
    sequences = dna_rna_tools.check_user_input(sequences, procedure)
    processed_sequences = {}
    for sequence in sequences:
        processed_sequences[sequence] = dna_rna_tools.DNA_RNA_PROCEDURES_FUNCTIONS[
            procedure
        ](sequence)
    if len(processed_sequences) == 1:
        return list(processed_sequences.values())[0]
    return processed_sequences


def run_protein_tools(sequences: (str, tuple[str] or list[str]), **kwargs: str) -> dict:
    """
    Process protein sequence by one of the developed tools.\n
    Run one procedure at a time:
    - Search for conserved amino acids residues in protein sequence
    - Search for alternative frames in a protein sequences
    - Convert protein sequences to RNA or DNA sequences
    - Reverse the protein sequences from one-letter to three-letter format and vice-versa
    - Define molecular weight of the protein sequences

    All functions except *search_for_alt_frames* are letter case sensitive\n
    If only one sequence provided - *sequences* can be string.\n
    If more - please provide *sequences* as list or tuple.\n
    Provide protein sequence in one letter code.\n
    You can obtain one letter code from three letter code with *three_one_letter_code*\n
    If more information needed please see README or desired docstring

    Arguments:
    - sequences (str, list[str] or tuple[str]): sequences to process
    - procedure (str]: desired procedure:
        - "search_for_motifs"
        - "search_for_alt_frames"
        - "convert_to_nucl_acids"
        - "three_one_letter_code"
        - "define_molecular_weight"

    For "search_for_motif" procedure provide:
    - motif (str]: desired motif to check presense in every given sequence\n
            Example: motif = "GA"
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)\n
            Example: overlapping = False

    For "search_for_alt_frames" procedure provide:
    - alt_start_aa (str]: the name of an amino acid that is encoded by alternative start codon (Optional)\n
            Example: alt_start_aa = "I"

    For "convert_to_nucl_acids" procedure provide:
    - nucl_acids (str]: the nucleic acid to convert to\n
            Example: nucl_acids = "RNA"\n
                           nucl_acids = "DNA"\n
                           nucl_acids = "both"

    Return:
    - dict: Dictionary with processed sequences. Depends on desired tool\n
            Please see Readme or desired docstring
    """
    procedure_arguments, procedure = protein_tools.check_and_parse_user_input(
        sequences, **kwargs
    )
    return protein_tools.PROTEINS_PROCEDURES_TO_FUNCTIONS[procedure](
        **procedure_arguments
    )


def run_fastq_filter(
    seqs: dict[str, tuple[str] | list[str]] = {
        "@SRX079804:1:SRR292678:1:1101:21885:21885": (
            "ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA",
            "FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD",
        ),
        "@SRX079804:1:SRR292678:1:1101:30161:30161": (
            "GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC",
            "DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD",
        ),
    },
    gc_bounds: (int | float | tuple[int | float] | list[int | float]) = (0, 100),
    length_bounds: (tuple[int]) = (0, 2**32),
    quality_threshold: (int) = 0,
    verbose=True,
) -> dict:
    """
    Filter out fastq reads by several parameters:
    - GC content
    - length
    - mean nucleotide quality
    Please provide GC content bounds in percentages.\n
    Please provide threshold quality in Phred-33 scale.\n
    For GC content and length bounds: if one number is provided, bounds from 0 to number are considered.\n

    Arguments:
    - seqs (dict[str, tuple[str] | list[str]]): fastq reads to be filtered
    - gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds
    - length_bounds (int | tuple[int] | list[int]): read length thresholds
    - quality_thresholds: read Phred-33 scaled quality thresholds
    Examples are provided as default values for each argument

    Return:
    - seqs_filtered (dict): similar dictionary as input, bad reads are filtered out
    """
    (
        seqs,
        gc_bounds,
        length_bounds,
        quality_threshold,
        verbose,
    ) = fastq_tools.check_user_input(
        seqs, gc_bounds, length_bounds, quality_threshold, verbose
    )
    return fastq_tools.fastq_filter(
        seqs, gc_bounds, length_bounds, quality_threshold, verbose
    )
