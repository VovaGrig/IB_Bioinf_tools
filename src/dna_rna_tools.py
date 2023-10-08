import dictionaries


def check_user_input(sequences: str or list[str] or tuple[str], procedure):
    """
    Check if user input can be correctly processed.\n

    Arguments:
    - sequences (str list[str] or tuple[str]): sequences to process.

    Return:
    - sequences (list): sequences to process.
    """
    if procedure not in DNA_RNA_PROCEDURES_FUNCTIONS.keys():
        raise ValueError("Wrong procedure")
    if isinstance(sequences, str):
        sequences = sequences.split()
    for sequence in sequences:
        if not set(sequence.upper()).issubset(
            dictionaries.NUCL_ACIDS_COMPLEMENT_RULE.keys()
        ):
            raise ValueError("Invalid sequence given")
        if "T" in sequence.upper() and "U" in sequence.upper():
            raise ValueError("Neither DNA or RNA sequence given")
    return sequences


def isrna(sequence: str) -> bool:
    """
    Check if given sequence is RNA sequence.\n

    Arguments:
    - sequence (str): sequence to check.

    Return:
    - boolean True/False.
    """
    return "U" in sequence or "u" in sequence


def transcribe(sequence):
    """
    Transcribe nucleic acid sequence.\n
    Also works for reverse transcription.\n
    Letter case sensitive.\n

    Arguments:
    - sequence (str): sequence to transcribe

    Return:
    - trancript (sequence): transcribed sequence
    """
    if isrna(sequence):
        transcript = sequence.replace("U", "T").replace("u", "t")
    else:
        transcript = sequence.replace("T", "U").replace("t", "u")
    return transcript


def reverse(sequence: str) -> str:
    """
    Convert nucleic acid sequence into its reverse counterpart.\n
    Letter case sensitive.\n

    Arguments:
    - sequence (str): sequence to convert

    Return:
    - str: reversed sequence
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Convert nucleic acid sequence into its complement counterpart.\n
    Letter case sensitive.\n

    Arguments:
    - sequence (str): sequence to convert

    Return:
    - str: comlepement sequence
    """
    if isrna(sequence):
        dictionaries.NUCL_ACIDS_COMPLEMENT_RULE["A"] = "U"
    complement_strand = ""
    for letter in sequence:
        if letter.islower():
            complement_strand += dictionaries.NUCL_ACIDS_COMPLEMENT_RULE[letter].lower()
        else:
            complement_strand += dictionaries.NUCL_ACIDS_COMPLEMENT_RULE[letter]
    dictionaries.NUCL_ACIDS_COMPLEMENT_RULE["A"] = "T"
    return complement_strand


def reverse_complement(sequence: str) -> str:
    """
    Convert nucleic acid sequence into its reverse-complement counterpart.\n
    Letter case sensitive.\n

    Arguments:
    - sequence (str): sequence to convert

    Return:
    - reverse_complement_strand (str): reverse-comlepement sequence
    """
    reverse_complement_strand = complement(sequence)
    reverse_complement_strand = reverse(reverse_complement_strand)
    return reverse_complement_strand


DNA_RNA_PROCEDURES_FUNCTIONS = {
    "transcribe": transcribe,
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
}
