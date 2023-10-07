complements_rule = {"A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a",
            "g": "c", "c": "g", "U": "A", "u": "a"}


def check_user_input(sequence):
    if not all(i in "".join(complements_rule.keys()) for i in sequence):
        raise ValueError("Invalid sequence given")
    if "T" in sequence.upper() and "U" in sequence.upper():
        raise ValueError("New nucleic acid discovered, go get Nobel")


def isrna(sequence):
    return ("U" in sequence or "u" in sequence)


def transcribe(sequence):
    if isrna(sequence):
        transcript = sequence.replace("U", "T").replace("u", "t")
    else:
        transcript = sequence.replace("T", "U").replace("t", "u")
    return transcript


def reverse(sequence):
    return sequence[::-1]


def complement(sequence):
    if isrna(sequence):
        complements_rule["A"] = "U"
        complements_rule["a"] = "u"
    complement_strand = ""
    for letter in sequence:
        complement_strand += complements_rule[letter]
    complements_rule["A"] = "T"
    complements_rule["a"] = "t"
    return complement_strand


def reverse_complement(sequence):
    if isrna(sequence):
        complements_rule["A"] = "U"
        complements_rule["a"] = "u"
    reverse_complement_strand = ""
    complement_strand = sequence[::-1]
    for letter in complement_strand:
        reverse_complement_strand += complements_rule[letter]
    complements_rule["A"] = "T"
    complements_rule["a"] = "t"
    return reverse_complement_strand


procedures_functions = {"transcribe": transcribe,
                        "reverse": reverse,
                        "complement": complement,
                        "reverse_complement": reverse_complement}


def run_dna_rna_tools(*args):
    procedure = args[-1]
    if procedure not in procedures_functions.keys():
        raise ValueError("Wrong procedure")
    sequences = list(args[:-1])
    processed_sequences = []
    for sequence in sequences:
        check_user_input(sequence)
        processed_sequences.append(procedures_functions[procedure](sequence))
    if len(processed_sequences) == 1:
        return processed_sequences[0]
    else:
        return processed_sequences
