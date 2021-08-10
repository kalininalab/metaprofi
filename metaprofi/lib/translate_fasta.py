"""MetaProFi six-frame translation module
"""

from metaprofi.lib.utilities_cython import reverse_complement

# fmt: off
DNA_translation_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}
# fmt: on


def six_frame_translation(seq, seq_name, min_len, outfile_name=None, out_file=False):
    """Performs six-frame translation and writes all translated sequences of size greater than minimum length to the output file

    Args:
      seq: Sequence
      seq_name: Name of the sequence
      min_len: K-mer size
      outfile_name: Filename/Path to store the six-frame translated sequence (Default value = None)
      out_file: Boolean, whether to write output to a file or not (Default value = False)

    """
    six_frames = dict()

    # Forward strand translation
    for idx, frame in enumerate(translate_seq(seq), start=0):
        if len(frame) >= min_len:
            six_frames[f"{seq_name}_forward_reading_frame_{idx}"] = frame

    # Reverse strand translation
    # Get the reverse complement of the sequence
    reverse_seq = reverse_complement(seq)
    for idx, frame in enumerate(translate_seq(reverse_seq), start=0):
        if len(frame) >= min_len:
            six_frames[f"{seq_name}_reverse_reading_frame_{idx}"] = frame

    # Output
    if not six_frames:
        return None
    if out_file and six_frames:
        write_to_file(six_frames, outfile_name)
        return "Written to file"
    return six_frames


def translate_seq(seq):
    """Three-frame translation of any given sequence

    Args:
      seq: Nucleotide (DNA) sequence (forward/reverse strand)

    Note:
      1. Does not terminate translation in any frame when a stop codon is hit, instead the translation continues by replacing stop codon with '_'
      2. Ambiguous codons are translated to 'X'

    Returns: Three translated sequences as a list

    """
    reading_frames = []

    for frame in range(3):
        reading_frames.append(
            "".join(
                [
                    DNA_translation_table.get(seq[pos : pos + 3], "X")
                    for pos in range(frame, len(seq) - 2, 3)
                ]
            )
        )
    return reading_frames


def write_to_file(translated_seqs, outfile_name):
    """Write the translated sequences to a fasta file (uncompressed)

    Args:
      translated_seqs: Six-frame translated sequences
      outfile_name: Output filename to store the six-frame translated sequences

    """
    with open(outfile_name, "a") as writer:
        for seq_name, seq in translated_seqs.items():
            writer.write(f">{seq_name}\n{seq}\n")
