#! /usr/bin/env python3

import sys
import translate
import find_orf

#definig start-stop codons and genetic codes

def translate_first_orf(
    start_codons = ['AUG'],
    stop_codons = ['UAA', 'UGA', 'UAG'],
    genetic_code = {
        'GUC': 'V', 'ACC': 'T', 'GUA': 'V','GUG': 'V', 'ACU': 'T', 'AAC': 'N', 
        'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 
        'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 
        'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 
        'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 
        'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 
        'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 
        'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 
        'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 
        'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'
        },
    ):
    """
   Translation of the orf in a RNA sequence. 

    Parameters:
    -------------------
    sequence : string containing RNA bases building the orf to be translated
    start_codons : string containing start codons
    stop_codons : strings containing stop codons
    genetic_code : triplet of bases composing a codon to specify amino acid abbreviation 
    """
    orf = find_orf.find_first_orf(sequence,
            start_codons = start_codons,
            stop_codons = stop_codons) 
    amino_acid_seq = translate.translate_sequence(orf, genetic_code)
    return amino_acid_seq


def main():
    import argparse

    # a command line parser object 
    parser = argparse.ArgumentParser()

    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']

    # specifying command-line arguments to the parser
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A start codon. This option can be used multiple times '
                    'if there are multiple start codons. '
                    'Default: {0}.'.format(" ".join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A stop codon. This option can be used multiple times '
                    'if there are multiple stop codons. '
                    'Default: {0}.'.format(" ".join(default_stop_codons))))

    args = parser.parse_args()

    # Checking whether the path option was set to true by the caller. Checking whether if start/stop codons were provided by the caller.
    if args.path:
        sequence = find_orf.parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons

    orf = find_orf.find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)
    sys.stdout.write('{}\n'.format(orf))

if __name__ == '__main__':
    main()
