#!/usr/bin/env python

import argparse
import subprocess

import textwrap as _textwrap
class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 1000)

        
parser = argparse.ArgumentParser(prog="phasstools",formatter_class=LineWrapRawTextHelpFormatter)
subparsers = parser.add_subparsers(help="command")
parser_1 = subparsers.add_parser("chromosome", help='phase, break, and scaffold and input assembly',formatter_class=LineWrapRawTextHelpFormatter)
parser_1.add_argument("--output", "-o", required=True, type=str, help="output directory")
parser_1.add_argument("--fasta", "-f", required= True, type=str, help="contigs to phase / phasing aware scaffold")
parser_1.add_argument("--hic_reads", "-i", required=True, type=str, help="hic fofn (file of file names) with format R1_filename\nR2_filename etc")
parser_1.add_argument("--kmer_size", "-k", required=False,type=int, default=31, help="kmer size (default 31)")
parser_1.add_argument("--linked_reads", "-l", required=False, type=str, help="linked reads fofn (file of file names) in format R1_filename\ttrim\nR1_filename\ttrim etc")
parser_1.add_argument("--ccs_reads", "-c", required=False, type=str, help="ccs reads fofn")
parser_1.add_argument("--short_reads", "-s", required=False, type=str, help= "short read data fofn alternating read1 and read2 filenames")
parser_1.add_argument("--kmer_data", "-d", required=False, type=str, default="linked_reads",help="which data to use for kmer spectrum and het kmer detection, must be linked_reads, short_reads, or ccs_reads")

parser_1.set_defaults(parser1=True)

parser_2 = subparsers.add_parser('phase', help='phase an individual with an existing reference or assembly')
parser_2.set_defaults(parser2=True)


args = parser.parse_args()

def chromosome():
    print("works")

print(args)
if args.parser1:
    chromosome()

