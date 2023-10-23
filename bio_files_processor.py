import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = ""):
    """
    Converts a multiline FASTA file to a one-line FASTA file.

    Arguments:
    - input_fasta (str): Path to the input FASTA file.
    - output_fasta (str): Path to the output FASTA file.\\
    If output_fasta is not provided, output file will be saved in the current directory with the same name as input file but with "_oneline.fasta" suffix.

    Return:
        None

    Example:
        convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
    """
    if output_fasta == "":
        output_fasta = os.path.basename(input_fasta)
        output_fasta = output_fasta.replace(".fasta", "_oneline.fasta")
    output_path = os.path.join("./", output_fasta)
    with open(input_fasta, "r") as input:
        with open(output_path, "w"):
            read = []
            while True:
                line = input.readline().strip()
                if not line:
                    break
                if line.startswith(">"):
                    line += "\n"
                    if read:
                        with open(output_path, "a") as output:
                            output.write("".join(read) + "\n")
                    read = [line]
                else:
                    read.append(line)
            with open(output_path, "a") as output:
                output.write("".join(read))


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: str or tuple[str] or list[str],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "",
):
    """
    Extract the translations of neighbour genes in amino acids from a GenBank (gbk) file for a list of desired genes and saves them in a FASTA file format ready for blasting.

    Arguments:
    - input_gbk (str): Path to the input GenBank (gbk) file.
    - genes (str or tuple[str] or list[str]): Gene name or list of gene names to extract neigbhour genes translations for.
    - n_before (int, optional): Number of genes to include before the desired gene. Defaults to 1.
    - n_after (int, optional): Number of genes to include after the desired gene. Defaults to 1.
    - output_fasta (str, optional): Path to the output FASTA file.\\
    If output_fasts is not provided, output file will be saved in the current directory with the same name as input file but with "_trans_for_blast.fasta" suffix.

    Return:
        None

    Example:
        select_genes_from_gbk_to_fasta("input.gbk", ["gene1", "gene2"], n_before=2, n_after=2, output_fasta="output.fasta")
    """
    if output_fasta == "":
        output_fasta = os.path.basename(input_gbk)
        output_fasta = output_fasta.replace(".gbk", "_trans_for_blast.fasta")
    output_path = os.path.join("./", output_fasta)
    if isinstance(genes, str):
        genes = [genes]
    with open(input_gbk, "r") as input:
        translations = []
        ind = 0
        prev_gene = "undefined"
        was_printed = set()
        while True:
            line = input.readline().strip()
            if not line:
                break
            if line.startswith("/gene="):
                prev_gene = line[7:-1]
            if line.startswith("/translation="):
                translation = line[14:]
                while not translation.endswith('"'):
                    line = input.readline().strip()
                    translation += line
                translation = translation.rstrip('"')
                translations.append((prev_gene, translation))
                prev_gene = "undefined"
        with open(output_path, "w") as out:
            for gene in genes:
                for ind, cur in enumerate(translations):
                    if gene in cur[0]:
                        for i in range(
                            max(0, ind - n_before),
                            min(len(translations), ind + n_after + 1),
                        ):
                            if i == ind or i in was_printed:
                                continue
                            was_printed.add(i)
                            out_gene, translation = translations[i]
                            out.write(f">{out_gene}\n{translation}\n")


select_genes_from_gbk_to_fasta(
    "/mnt/c/users/vovag/Documents/IB/Python_1sem/HW5_Grigoriants/example_gbk.gbk",
    genes="ybgD_2",
    n_before=10,
    n_after=10,
)


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = ""):
    """
    Shift the starting position of each sequence in a fasta file by a specified number of nucleotides and saves the modified sequences in a new fasta file.

    Arguments:
    - input_fasta (str): Path to the input fasta file.
    - shift (int): Number of nucleotides to shift the starting position of each sequence.
    - output_fasta (str, optional): Path to the output fasta file.\\
    If output_fasta is not provided, output file will be saved in the current directory with the same name as input file but with "_shifted.fasta" suffix.

    Return:
        None

    Example:
        change_fasta_start_pos("input.fasta", 3, "output.fasta")
    """
    if output_fasta == "":
        output_fasta = os.path.basename(input_fasta)
        output_fasta = output_fasta.replace(".fasta", "_shifted.fasta")
    output_path = os.path.join("./", output_fasta)
    with open(input_fasta, "r") as input:
        with open(output_path, "w"):
            while True:
                line = input.readline().strip()
                print(line)
                if not line:
                    break
                if not line.startswith(">"):
                    line = line[shift:] + line[:shift]
                with open(output_path, "a") as output:
                    output.write(line + "\n")
