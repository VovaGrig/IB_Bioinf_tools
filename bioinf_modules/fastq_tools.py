import dictionaries


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

    Return:
    - seqs_filtered (dict): similar dictionary as input, bad reads are filtered out
    """
    seqs, gc_bounds, length_bounds, quality_threshold = check_user_input(
        seqs, gc_bounds, length_bounds, quality_threshold
    )
    return fastq_filter(seqs, gc_bounds, length_bounds, quality_threshold)


def check_user_input(
    seqs: dict[str, tuple[str] | list[str]],
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
    length_bounds: (int | tuple[int] | list[int]),
    quality_threshold: int,
):
    """
    Check if user input can be correctly processed\n

    Arguments:
    - seqs (dict[str, tuple[str] | list[str]]): fastq reads to be filtered
    - gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds
    - length_bounds (int | tuple[int] | list[int]): read length thresholds
    - quality_thresholds: read Phred-33 scaled quality thresholds

    Return:
    - same arguments as input, checked for correctness
    """
    if not isinstance(seqs, dict):
        raise ValueError("Please provide sequences info with a dictionary")
    for seq_name in seqs.keys():
        if not isinstance(seq_name, str):
            raise ValueError("Invalid sequence name given")
        if not set(seqs[seq_name][0].upper()).issubset({"A", "T", "G", "C"}):
            raise ValueError("Invalid read sequence given")
        if not set(seqs[seq_name][1]).issubset(dictionaries.QUALITY_SYMBOLS):
            raise ValueError("Invalid quality sequence given")
    return seqs, gc_bounds, length_bounds, quality_threshold


def fastq_filter(
    seqs: dict[str, tuple[str] | list[str]],
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
    length_bounds: (int | tuple[int] | list[int]),
    quality_threshold: (int | float),
) -> dict:
    """
    Parse checked input and filter out bad reads.\n

    Arguments:
    - seqs (dict[str, tuple[str] | list[str]]): fastq reads to be filtered
    - gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds
    - length_bounds (int | tuple[int] | list[int]): read length thresholds
    - quality_thresholds: read Phred-33 scaled quality thresholds

    Return:
    - seqs_filtered (dict): similar dictionary as input, bad reads are filtered out
    """
    seqs_filtered = {}
    for seq_name in seqs.keys():
        seq = seqs[seq_name][0]
        seq_qual = seqs[seq_name][1]
        gc_result = is_gc_good(seq, gc_bounds)
        len_result = is_len_good(seq, length_bounds)
        qual_result = is_qual_good(seq_qual, quality_threshold)
        if gc_result and len_result and qual_result:
            seqs_filtered[seq_name] = seqs[seq_name]
    return seqs_filtered


NEW_LINE = "\n"  # needed for output in f-strings


def is_gc_good(
    seq: str,
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
) -> bool:
    """
    Check GC content of a given read.\n
    Please provide bounds in percentages.\n
    If one number is provided, bounds from 0 to number are considered.\n

    Arguments:
    - seq (str): read to check it's length
    - gc_bounds (int | float | tuple[int] | list[int]): GC content thresholds, by which reads are filtered

    Return:
    - condition (bool): True - GC content is within bounds, False - read is to be filtered
    """
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    gc_content = seq.count("G") + seq.count("C")
    gc_content = gc_content / len(seq) * 100
    print(f"Read: {seq}")
    print(f"GC Content: {round(gc_content, 4)}")
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def is_len_good(seq: str, length_bounds: (int | tuple[int] | list[int])) -> bool:
    """
    Check length of a given read.\n
    If one number is provided, bounds from 1 to number are considered.\n

    Arguments:
    - seq (str): read to check it's length
    - length_bounds (int | tuple[int] | list[int]): length thresholds, by which reads are filtered

    Return:
    - condition (bool): True - length is within bounds, False - read is to be filtered
    """
    if isinstance(length_bounds, int):
        length_bounds = (1, length_bounds)
    print(f"Read Length: {round(len(seq), 4)}")
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_qual_good(seq_qual: str, quality_threshold: int | float) -> bool:
    """
    Check mean nucleotide sequencing quality for a given read.\n
    Reads with mean quality less then threshold are filtered out.\n
    Please provide Phred-33 Q-Scores as input.\n

    Arguments:
    - seq_qual (str): quality sequence in Phred-33 scale for a read
    - quality_threshold (int | float): threshold, by which reads are filtered

    Return:
    - condition (bool): True - quality is above threshold, False - read is to be filtered
    """
    mean_quality = sum(ord(symbol) - 33 for symbol in seq_qual) / len(seq_qual)
    print(f"Mean Nucleotide Quality: {round(mean_quality, 4)}{NEW_LINE}")
    return mean_quality > quality_threshold
