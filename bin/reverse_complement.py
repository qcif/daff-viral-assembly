#!/usr/bin/env python
import argparse

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")
    
    # All the required arguments #
    parser.add_argument("--ids_to_rc", type=str)
    parser.add_argument("--sample", type=str)
    parser.add_argument("--fasta", type=str)
    args = parser.parse_args()

    ids_to_rc = args.ids_to_rc
    sample = args.sample
    fasta = args.fasta


    
    #raw_data = pd.read_csv(results_path, header=0, sep="\t",index_col=None)
    with open(ids_to_rc, 'r') as f:
        lines = []
        for line in f:
            lines.append(line.strip())
    

    contig_dict = {}
    with open(fasta) as file:
         
         for line in file:
            if line.startswith(">"):
                header = line.strip().replace(">","")
                seq_fasta = next(file).strip()

                if header in lines:
                    contig_dict[header] = reverse_complement(seq_fasta)
                else:
                    contig_dict[header] = seq_fasta

    ofile =  open(sample + "_final_polished_consensus_rc.fasta", "w")
    
    for seq_id,fasta in contig_dict.items():
        identifier_line = ">" + seq_id + "\n"
        ofile.write(identifier_line)
        sequence_line = fasta + "\n"
        ofile.write(sequence_line)

    ofile.close()


def reverse_complement(seq):
    alt_map = {'N':'0'}
    complement = {'A': 'T',
                  'T': 'A',
                  'C': 'G', 
                  'G': 'C', 
                  'N': 'N',
                  'R': 'Y',
                  'Y': 'R',
                  'K': 'M',
                  'M': 'K',
                  'B': 'V',
                  'V': 'B',
                  'H': 'D',
                  'D': 'H',
                  'S': 'S',
                  'W': 'W',
                  'a': 't',
                  't': 'a',
                  'c': 'g',
                  'g': 'c',
                  'n': 'n',
                  'r': 'y',
                  'y': 'r',
                  'k': 'm',
                  'm': 'k',
                  'b': 'v',
                  'v': 'b',
                  'h': 'd',
                  'd': 'h',
                  's': 's',
                  'w': 'w'}
        
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

if __name__ == "__main__":
    main()