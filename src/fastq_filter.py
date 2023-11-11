if __name__ == "__main__":
    import dictionaries
else:
    from src import dictionaries


def check_user_input(
    seqs: dict[str, tuple[str] | list[str]],
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
    length_bounds: (int | tuple[int] | list[int]),
    quality_threshold: int,
    verbose,
):
    """
    Check if user input can be correctly processed\n

    Arguments:
    - seqs (dict[str, tuple[str] | list[str]]): fastq reads to be filtered
    - gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds
    - length_bounds (int | tuple[int] | list[int]): read length thresholds
    - quality_thresholds: read Phred-33 scaled quality thresholds
    - verbose (bool): add detailed statistics for each read

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
    if verbose != True and verbose != False:
        raise ValueError("Invalid *verbose* argument given")
    return seqs, gc_bounds, length_bounds, quality_threshold, verbose


def fastq_filter(
    seqs: dict[str, tuple[str] | list[str]],
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
    length_bounds: (int | tuple[int] | list[int]),
    quality_threshold: (int | float),
    verbose,
) -> dict:
    """
    Parse checked input and filter out bad reads.\n

    Arguments:
    - seqs (dict[str, tuple[str] | list[str]]): fastq reads to be filtered
    - gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds
    - length_bounds (int | tuple[int] | list[int]): read length thresholds
    - quality_thresholds: read Phred-33 scaled quality thresholds
    - verbose (bool): add detailed statistics for each read

    Return:
    - seqs_filtered (dict): similar dictionary as input, bad reads are filtered out
    """
    seqs_filtered = {}
    for seq_name in seqs.keys():
        seq = seqs[seq_name][0]
        seq_qual = seqs[seq_name][1]
        gc_result = is_gc_good(seq, gc_bounds, verbose)
        if gc_result:
            len_result = is_len_good(seq, length_bounds, verbose)
            if len_result:
                qual_result = is_qual_good(seq_qual, quality_threshold, verbose)
                if qual_result:
                    seqs_filtered[seq_name] = seqs[seq_name]
    return seqs_filtered


NEW_LINE = "\n"  # needed for output in f-strings


def is_gc_good(
    seq: str, gc_bounds: (int | float | tuple[int | float] | list[int | float]), verbose
) -> bool:
    """
    Check GC content of a given read.\n
    Please provide bounds in percentages.\n
    If one number is provided, bounds from 0 to number are considered.\n

    Arguments:
    - seq (str): read to check it's length
    - gc_bounds (int | float | tuple[int] | list[int]): GC content thresholds, by which reads are filtered
    - verbose (bool): add detailed statistics for each read

    Return:
    - condition (bool): True - GC content is within bounds, False - read is to be filtered
    """
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    gc_content = seq.count("G") + seq.count("C")
    gc_content = gc_content / len(seq) * 100
    if verbose:
        print(f"Read: {seq}")
        print(f"GC Content: {round(gc_content, 4)}")
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def is_len_good(
    seq: str, length_bounds: (int | tuple[int] | list[int]), verbose
) -> bool:
    """
    Check length of a given read.\n
    If one number is provided, bounds from 1 to number are considered.\n

    Arguments:
    - seq (str): read to check it's length
    - length_bounds (int | tuple[int] | list[int]): length thresholds, by which reads are filtered
    - verbose (bool): add detailed statistics for each read

    Return:
    - condition (bool): True - length is within bounds, False - read is to be filtered
    """
    if isinstance(length_bounds, int):
        length_bounds = (1, length_bounds)
    if verbose:
        print(f"Read Length: {round(len(seq), 4)}")
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_qual_good(seq_qual: str, quality_threshold: int | float, verbose) -> bool:
    """
    Check mean nucleotide sequencing quality for a given read.\n
    Reads with mean quality less then threshold are filtered out.\n
    Please provide Phred-33 Q-Scores as input.\n

    Arguments:
    - seq_qual (str): quality sequence in Phred-33 scale for a read
    - quality_threshold (int | float): threshold, by which reads are filtered
    - verbose (bool): add detailed statistics for each read

    Return:
    - condition (bool): True - quality is above threshold, False - read is to be filtered
    """
    mean_quality = sum(ord(symbol) - 33 for symbol in seq_qual) / len(seq_qual)
    if verbose:
        print(f"Mean Nucleotide Quality: {round(mean_quality, 4)}{NEW_LINE}")
    return mean_quality > quality_threshold
