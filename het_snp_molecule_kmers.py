#!/usr/bin/env python

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description = "find heterozygous kmers from short read data")

parser.add_argument("--txg_reads", required=False, help = "10x genomics fofn")
parser.add_argument("--hic_reads", required=False, help = "Hi-C reads fofn")
parser.add_argument("--long_reads", required=False, help= "long reads fofn")
parser.add_argument("--fasta","-f", required=False, help= "fasta file")
parser.add_argument("-o", "--output", required=True, help = "prefix output")
parser.add_argument("-t", "--threads", required=True, help="threads")
parser.add_argument("--min_coverage", required=True, type= int, help="min coverage for het kmer")
parser.add_argument("--max_coverage", required=True, type=int, help= "max coverage for het kmer")
parser.add_argument("--max_total_coverage", required=True, type=int, help ="max total coverage for het kmer pair")
parser.add_argument("--max_error", required=False, default=4, type=int, help="max coverage of 3rd kmer with same outer bases")
parser.add_argument("--whitelist", required=False, help="10x barcode whitelist")
parser.add_argument("--ploidy", required=False, type=int, default=2, help="ploidy")
parser.add_argument("-m", "--memory", required=False, type=int, default = 24, help="memory in GB, default 24")
parser.add_argument("--unpaired_het_modimizer", required=False, default="23")
parser.add_argument("--hom_modimizer",required=False, default="53")
parser.add_argument("--kmer_size",required=True, type=int, help = "kmer length")
args = parser.parse_args()

mypath = os.path.dirname(os.path.realpath(__file__))
cmd = [mypath+"/het_snp_kmers/target/release/het_snp_kmers", "--kmer_counts", args.output+"/kmer_counts.tsv"]
cmd.extend(["--max_coverage", str(args.max_coverage), "--min_coverage", str(args.min_coverage), 
            "--max_error", str(args.max_error), "--max_total_coverage", str(args.max_total_coverage),
            "--unpaired_het_modimizer", args.unpaired_het_modimizer, "--hom_modimizer", args.hom_modimizer])
print(" ".join(cmd))
with open(args.output+"/het_kmers.tsv", 'w') as out:
    subprocess.check_call(cmd,stdout=out)

cmd = [mypath+"/molecule_kmers/target/release/molecule_kmers", "--output", args.output, "--kmer_size", str(args.kmer_size)]
cmd.extend(["--txg_barcodes", args.whitelist])
if not(args.txg_reads is None):
    cmd.extend(["--txg_reads", args.txg_reads])
if args.hic_r1s:
    cmd.extend(["--hic_reads", args.hic_reads])

if args.long_reads:
    cmd.extend(["--long_reads", args.long_reads])
cmd.extend(["--paired_kmers", args.output + "/het_kmers.tsv"])

print("running molecule kmers now")
print(" ".join(cmd))
with open(args.output + "/molecule_kmers.err", 'w') as err:
    subprocess.check_call(cmd, stderr = err)

mol_out_txg = [args.output+"/molecules_txg_"+"{0:0=3d}".format(x)+".bin" for x in range(len(args.txg_r1s))]
mol_out_longreads = [args.output+"/molecules_longreads_"+"{0:0=3d}".format(x)+".bin" for x in range(len(args.long_reads))]
if not(args.hic_r1s is None):
    mol_out_hic = [args.output+"/molecules_hic_"+"{0:0=3d}".format(x)+".bin" for x in range(len(args.hic_r1s))]

print("\t".join(mol_out_txg))
print("\t".join(mol_out_longreads))
print("\t".join(mol_out_hic))


cmd = [mypath + "/mixhap/target/release/mixhap"]
cmd.extend(["--kmers", args.output+"/het_kmers.tsv"])
cmd.extend(["--txg_mols"]+mol_out_txg)
cmd.extend(["--longread_mols"]+mol_out_longreads)
if not(args.hic_r1s is None):
    cmd.extend(["--hic_mols"]+mol_out_hic)

cmd.extend(["--longread_fqs"]+args.long_reads)
cmd.extend(["--ploidy", str(args.ploidy), "--variants", args.output + "/molecule_kmers.custom_binary"])
cmd.extend(["--barcode_whitelist", args.whitelist])
print("command for phasing")
print(" ".join(cmd))
#with open(args.output+"/something_to_be_determined.out", 'w') as out:
#    subprocess.check_call(cmd, stdout = out)

