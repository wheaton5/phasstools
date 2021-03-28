#!/usr/bin/env python

import argparse
import subprocess
import os

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
parser_1.add_argument("--barcode_whitelist", required=False, type=str, help="10x barcode whitelist")
parser_1.add_argument("--ccs_reads", "-c", required=False, type=str, help="ccs reads fofn")
parser_1.add_argument("--short_reads", "-s", required=False, type=str, help= "short read data fofn alternating read1 and read2 filenames")
parser_1.add_argument("--kmer_data", "-d", required=False, type=str, default="linked_reads",help="which data to use for kmer spectrum and het kmer detection, must be linked_reads, short_reads, or ccs_reads")
parser_1.add_argument("--threads", "-t", required=False, type=int, default=8, help="number of threads to use at maximum")
parser_1.add_argument("--mem", "-m", required=False, type=int, default=20, help="max memory in Gb")
parser_1.set_defaults(parser1=True)

parser_2 = subparsers.add_parser('phase', help='phase an individual with an existing reference or assembly')
parser_2.set_defaults(parser2=True)


args = parser.parse_args()

directory = os.path.dirname(os.path.realpath(__file__))
if args.linked_reads:
    assert args.barcode_whitelist, "requires barcode whitelist if you have linked reads"


def het_kmers():
    if args.kmer_data == "linked_reads":
        fofn = args.linked_reads
    elif args.kmer_data == "short_reads":
        fofn = args.short_reads
    elif args.kmer_data == "ccs_reads":
        fofn = args.ccs_reads
    else:
        assert False, "kmer_data must be one of linked_reads, short_reads, or ccs_reads" 
    cmd = [directory + "/het_kmers.py", "-i", fofn, "-o", args.output, "-k", 
        str(args.kmer_size), "-t", str(args.threads), "-m", str(args.mem)]
    subprocess.check_call(cmd)

def purge_dups():
    cmd = [directory + "/purge_dups/bin/calcuts", args.output + "/hist.tsv"]
    with open(args.output+"/calcuts.out",'w') as out:
        with open(args.output+"/calcuts.err",'w') as err:
            subprocess.check_call(cmd, stdout=out, stderr = err)

def load_cutoffs():
    with open(args.output+"/calcuts.out") as cuts:
        for line in cuts:
            if line.startswith("["):
                continue
            else:
                toks = line.strip().split()
                minimum = int(toks[1])
                maximum = int(toks[2])
                return [minimum, maximum]

def het_kmer_molecules(cutoffs):
    cmd = [directory + "/het_snp_molecule_kmers.py", "--hic_reads", args.hic_reads,
        "--fasta", args.fasta, "--kmer_size", str(args.kmer_size), "--output", args.output,
        "-t", str(args.threads), "-m", str(args.mem), "--min_coverage", str(cutoffs[0]), 
        "--max_coverage", str(cutoffs[1]), "--max_total_coverage", str(cutoffs[1]*2)]
    if args.linked_reads:
        cmd.extend(["--txg_reads", args.linked_reads, "--whitelist", args.barcode_whitelist])
    if args.ccs_reads:
        cmd.extend(["--long_reads", args.ccs_reads])
    with open(args.output+"/het_kmer_molecules.out",'w') as out:
        with open(args.output+"/het_kmer_molecules.err",'w') as err:
            subprocess.check_call(cmd, stdout= out, stderr=err)
    
    

def phasing():
    cmd = [directory + "/phasst_phase/target/release/phasst_phase", "--output", args.output, 
        "--het_kmers", args.output + "/het_kmers.tsv", "--hic_mols", args.output + "/hic.fofn",
        "--threads", str(args.threads), "--assembly_fasta", args.fasta, "--assembly_kmers", 
        args.output + "/fasta_kmers.bin"]
    if args.linked_reads:
        cmd.extend(["--linked_read_mols", args.output + "/txg.fofn"])
    if args.ccs_reads:
        cmd.extend(["--long_read_mols", args.output + "/ccs.fofn"])

    with open(args.output + "/phasing.out", 'w') as out:
        with open(args.output + "/phasing.err", 'w') as err:
            subprocess.check_call(cmd, stdout = out, stderr = err)


def scaffolding():
    if not os.path.exists(args.output + "/breaks_kmers"):
        subprocess.check_call(["mkdir", args.output + "/breaks_kmers"])
    if not os.path.exists(args.output + "/breaks_kmers/fasta_kmers.bin"):
        cmd = [directory + "/molecule_kmers/target/release/molecule_kmers", "-o", args.output + "/breaks_kmers", 
            "--paired_kmers", args.output + "/het_kmers.tsv", "--fasta", args.output + "/breaks.fa", 
            "--kmer_size", str(args.kmer_size), "--threads", "1"]
        with open(args.output + "/breaks_kmers/fasta_kmers.err", 'w') as err:
            with open(args.output + "/breaks_kmers/fasta_kmers.out", 'w') as out:
                subprocess.check_call(cmd, stderr=err, stdout = out)
    cmd = [directory + "/phasst_scaff/target/release/phasst_scaff", "-o", args.output, "--het_kmers",
        args.output + "/het_kmers.tsv", "--linked_read_barcodes", args.output + "/txg.fofn",
        "--hic_mols", args.output + "/hic.fofn", "--assembly_fasta", args.output + "/breaks.fa",
        "--assembly_kmers", args.output + "/breaks_kmers/fasta_kmers.bin", "--phased_vcf", args.output + "/phasing_breaks.vcf"]
    
    with open(args.output + "/scaffolding.out", 'w') as out:
        with open(args.output + "/scaffolding.err", 'w') as err:
            err.write(" ".join(cmd)+"\n")
            subprocess.check_call(cmd, stdout = out, stderr = err)


def chromosome():
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(args.output + "/hist.tsv"):
        het_kmers()
    else:
        print("using previously generated kmer count data")
    if not os.path.exists(args.output + "/calcuts.out"):
        print("determining haploid coverage cutoffs")
        purge_dups()
    else:
        print("using previously generaged kmer coverage cutoffs")
    if not os.path.exists(args.output + "/fasta_kmers.bin"):
        cutoffs = load_cutoffs()
        print("finding het kmers on read data")
        het_kmer_molecules(cutoffs)
    else:
        print("using previously generated kmer molecule data")
    if not os.path.exists(args.output + "/breaks.fa"):
        print("phasing het kmers")
        phasing()
    else:
        print("using previously generated phasing")
    if not os.path.exists(args.output + "/scaff.tsv"):
        print("scaffolding with phased data")
        scaffolding()
    else:
        print("previous pipeline pointed at this directory is complete. If you want to rerun, specify a new output directory.")
    print("done.")



    




    

print(args)
if args.parser1:
    chromosome()

