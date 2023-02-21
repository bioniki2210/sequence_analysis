from search_for_sample import seq_array
import numpy as np

def degen_bases_sym():
    """Nucleotide ambiguity code (IUPAC)"""

    return {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G',
            'AT': 'w', 'CG': 's',
            'AC': 'm', 'GT': 'k',
            'AG': 'r', 'CT': 'y',
            'CGT': 'b', 'AGT': 'd', 'ACT': 'h', 'ACG': 'v',
            'ACGT': 'n'}

def get_primer_lyb(infile, primer_size=21, occurrence=0.8, degen_nucleo_percent=20):
    """Generate degenerative primers

    Return path to txt file with potential degenerative primers;
    Structure of outfile:
    _________________________________
    Sequence    Start   Stop    Score
    _________________________________

    infile - path to aligned fasta file with continuous sequences (type 'str');
    primer_size - length of degenerative primers (default 21 nb; type 'int');
    occurrence - frequency threshold. If the occurrence of a nucleotide is greater than
    or equal to the threshold value, then this nucleotide is considered
    valid for this position in the degenerate primer (type 'float'; default 0.8);
    degen_nucleo_percent - percentage of degenerate nucleotides in primers (type 'int', default 20%)"""

    sequences = seq_array(infile)
    c = np.array([list(i) for i in sequences], dtype=np.chararray)
    with open(infile.replace('.fasta', '_primer.txt'), 'w') as out:
        for i in range(len(sequences[0])):
            seq_slice = c[:, i:i + primer_size]

            # Primer coordinates
            start, stop = str(i), str(i + primer_size)
            primer = str()

            # Score of primer;
            # score = sum(max_frequency_i), where i is a position in primer;
            # max(score) = primer_size (highly conserved, no degenerate nucleotides);
            # min(score) = 0 (not conserved, all nucleotides are degenerate)
            score = 0
            if seq_slice.shape[1] == primer_size:
                # Filter out all fragments shorter than the set primer size

                # Generate nucleotide frequency model for fragments
                seq_model = [{j: round(np.count_nonzero(seq_slice[:, i] == j) / len(sequences), 3)
                              for j in sorted(set(seq_slice[:, i]))} for i in range(len(seq_slice[0, :]))]
                for n in seq_model:
                    max_freq_char = max(n, key=n.get)
                    score += n.get(max_freq_char)
                    if n.get(max_freq_char) >= occurrence:
                        # Add conserved nucleotide to primer
                        primer += str(max_freq_char)
                    else:
                        # Add degenerate nucleotide to primer
                        base_key = ''.join(sorted(n.keys())).replace('-', '')
                        primer += str(degen_bases_sym().get(base_key))
            else:
                continue
            if '--' not in primer:
                # Remove primers with two or more gaps
                deg_base_num = (len([c for c in primer if c.islower()]) * 100) / len(primer)
                if deg_base_num < degen_nucleo_percent:
                    # This condition is valid for primers containing
                    # no more than N% degenerate nucleotides (N = degen_nucleo_percent);
                    out.write(';'.join([primer, start, stop, str(round(score, 3))]) + '\n')
                else:
                    continue
            else:
                continue
    return infile.replace('.fasta', '_primer.txt')
