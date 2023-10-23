if __name__ == "__main__":
    import dictionaries
else:
    from src import dictionaries
import os


def parse_and_check_user_input(
    sequences_path: str,
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
    length_bounds: (int | tuple[int] | list[int]),
    quality_threshold: int,
    verbose: bool,
    output_filename: "str",
):
    """
    Parse input fasta file to dictionary[sequence_name: [sequence, sequence_quality]] and check if input can be correctly processed\n

    Arguments:
    - sequences_path(str): absolute or relative path to desired file, containing sequences in fasta format
    - gc_bounds (int | float | tuple[int | float] | list[int | float]): GC content thresholds
    - length_bounds (int | tuple[int] | list[int]): read length thresholds
    - quality_thresholds: read Phred-33 scaled quality thresholds
    - verbose (bool): add detailed statistics for each read
    - output_filename (str): output name of the filtered fasta file

    Return:
    - same arguments as input, checked for correctness
    """
    if output_filename == "":
        output_filename = os.path.basename(sequences_path)
    seqs = {}
    with open(sequences_path, "r") as seqs_file:
        count = 0
        for line in seqs_file:
            count += 1
            if count == 1 or (count - 1) % 4 == 0:
                seqs[line.strip()] = []
            if count % 2 == 0 and count % 4 != 0:
                seqs[list(seqs)[-1]].append(line.strip())
            if count % 4 == 0:
                seqs[list(seqs)[-1]].append(line.strip())
    for seq_name in seqs.keys():
        if not isinstance(seq_name, str):
            raise ValueError("Invalid sequence name given")
        if not set(seqs[seq_name][0].upper()).issubset({"A", "T", "G", "C"}):
            raise ValueError("Invalid read sequence given")
        if not set(seqs[seq_name][1]).issubset(dictionaries.QUALITY_SYMBOLS):
            raise ValueError("Invalid quality sequence given")
    if verbose != True and verbose != False:
        raise ValueError("Invalid *verbose* argument given")
    return seqs, gc_bounds, length_bounds, quality_threshold, verbose, output_filename


def fastq_filter(
    seqs: dict[str, list[str]],
    gc_bounds: (int | float | tuple[int | float] | list[int | float]),
    length_bounds: (int | tuple[int] | list[int]),
    quality_threshold: (int | float),
    verbose,
) -> dict:
    """
    Filter out bad reads.\n

    Arguments:
    - seqs (dict[str, list[str]]): fastq reads to be filtered
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


def save_filtered_seqs(
    filtered_seqs: dict[str, list[str]],
    output_filename: str,
    save_to_dir: str,
):
    """
    Save filtered reads in fasta format.

    Arguments:
    - filtered_seqs (dict[str, list[str]]): filtered fastq_reads
    - output_filename (str): output name of the filtered fasta file
    - save_to_dir (str): absolute or realtive path to directory to save to
    """
    os.makedirs(save_to_dir, exist_ok=True)
    output_path = os.path.join(save_to_dir, output_filename)
    for seq_name in filtered_seqs.keys():
        filtered_seqs[seq_name].append("+" + seq_name[1:])
        filtered_seqs[seq_name] = [x + "\n" for x in filtered_seqs[seq_name]]
    filtered_seqs[list(filtered_seqs)[-1]][1] = filtered_seqs[list(filtered_seqs)[-1]][
        1
    ].strip()  # get rid of last \n, to avoid empty line in filtered fasta file
    with open(output_path, "w") as output_file:
        for key, value in filtered_seqs.items():
            output_file.write("%s\n%s%s%s" % (key, value[0], value[2], value[1]))
