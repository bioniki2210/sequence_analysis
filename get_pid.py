from search_for_sample import seq_array, seq_name
import sys

def get_pid_table(infile, outfile):
    """Generate symmetric matrix of percentage identity (PID) between aligned sequences

    Return path to txt file with matrix of PID values;
    PID = Neq * 100 / (N1 + N2 - Neq), %
    where Neq is the number of identical aligned non-gapped characters;
    N1 and N2 are the length of two aligned sequences;

    infile - path to initial aligned fasta file with continuous sequences (type 'str');
    outfile - path to txt file containing symmetric matrix of PID values (type 'str')"""

    # Start progress counter
    counter = 0
    sequences = seq_array(infile)
    org_names = seq_name(infile)
    with open(outfile, 'w') as out:
        for i in range(len(sequences)):
            # Result string for seq_i
            res = [org_names[i]]
            seq_i = sequences[i]
            len_seq_i = len(seq_i.replace('-', ''))
            seq_i_list = list(seq_i)
            for j in range(len(sequences)):
                seq_j = sequences[j]
                len_seq_j = len(seq_j.replace('-', ''))
                seq_j_list = list(seq_j)

                # Initial identity of seq_1 and seq_2
                identical = 0
                for n in range(len(seq_i_list)):
                    if seq_i_list[n] == seq_j_list[n]:
                        a = list(zip(seq_i_list[n], seq_j_list[n]))
                        if a != [('-', '-')]:
                            # This condition is valid for non-gapped pair of characters
                            identical += 1
                        else:
                            continue
                    else:
                        continue
                pid = (identical / (len_seq_i + len_seq_j - identical)) * 100
                res.append(str(pid))

            # Write result string to txt outfile;
            # Delimiter ';', decimal separator '.'
            out.write(';'.join(res) + '\n')

            # Progress tracking
            counter += 1
            progress = str(round(counter * 100 / len(sequences), 1))

            # Progress bar
            sys.stdout.write('\rProgress: ' + progress + '%')
            sys.stdout.flush()
    return outfile