#!/usr/bin/env python

import argparse
import subprocess
import os
from pyfaidx import Fasta
import time
import sys
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import numpy as np

FASTK = False


import textwrap as _textwrap
class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 1000)

        
parser = argparse.ArgumentParser(prog="phasstools",formatter_class=LineWrapRawTextHelpFormatter)
subparsers = parser.add_subparsers(help="command")
parser_1 = subparsers.add_parser("scaff", help='phase, break, and scaffold and input assembly',formatter_class=LineWrapRawTextHelpFormatter)
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
parser_1.add_argument("--hom_modimizer", required=False, type=int, default=11, help="use hom kmers when hom % modimizer == 0")
parser_1.add_argument("--mem", "-m", required=False, type=int, default=20, help="max memory in Gb")
parser_1.add_argument("--min_kmer_count", required = False, type = int, default = 10, help = "min kmer count to use in kmer spectrum")
parser_1.set_defaults(parser1=True)

parser_2 = subparsers.add_parser('phase', help='phase an individual with an existing reference or assembly')
parser_2.set_defaults(parser2=True)


args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

directory = os.path.dirname(os.path.realpath(__file__))
if args.linked_reads:
    assert args.barcode_whitelist, "requires barcode whitelist if you have linked reads"


def check_call(cmd, base_out_name, shell = False):
    print("running "+base_out_name)
    print(" ".join(cmd))
    with open(args.output+"/"+base_out_name+".out",'w') as out:
        with open(args.output+"/"+base_out_name+".err", 'w') as err:
            subprocess.check_call(cmd, shell = shell, stdout = out, stderr = err)

def fofn_to_list(fofn):
    with open(fofn) as infile:
        files = [line.strip() for line in infile]
    return(files)


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
        str(args.kmer_size), "-t", str(args.threads), "-m", str(args.mem), "--min_count", str(args.min_kmer_count)]
    check_call(cmd, "kmc_call")

def count_kmers_FASTK():
    r1s = fofn_to_list(args.ccs_reads)
    max_p = 0
    threads = args.threads

    mem = args.mem
    name = "ccs_spectrum"
    cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), "-t"+str(args.min_kmer_count), 
    "-N"+args.output+"/"+name, "-M"+str(mem), "-T"+str(threads)] + r1s 
   
    
    check_call(cmd, name)
    name = "fasta_spectrum"
    cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), "-t1", 
    "-N"+args.output+"/"+name, "-M"+str(mem), "-T"+str(threads), args.fasta]
    check_call(cmd, name)

    name = "fasta_ccs_profile"
    cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), 
          "-N"+args.output+"/"+name, "-p:"+args.output+"/ccs_spectrum", "-M"+str(mem), "-T"+str(threads), args.fasta]
    check_call(cmd, name)

    name = "fasta_self_profile"
    cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), 
          "-N"+args.output+"/"+name, "-p:"+args.output+"/fasta_spectrum", "-M"+str(mem), "-T"+str(threads), args.fasta]
    check_call(cmd, name)
    
    fasta = Fasta(args.fasta)
    contigs = len(fasta.keys())


    #for contig in range(contigs):
    name = "fasta_ccs_profile_text"
    cmd = [directory+"/FASTK/Profex", args.output+"/fasta_ccs_profile"] + [str(x+1) for x in range(contigs)]
    check_call(cmd, name)

    name = "fasta_self_profile_text"
    cmd = [directory+"/FASTK/Profex", args.output+"/fasta_self_profile"] + [str(x+1) for x in range(contigs)]
    check_call(cmd, name)




    total = 0
    denom = 0
    contig_kmer_cov = []
    contig = 0
    with open(args.output+"/fasta_ccs_profile_text.out") as ccs_prof:
        with open(args.output + "/fasta_self_profile_text.out") as self_prof:
            for (line1, line2) in zip(ccs_prof, self_prof):
                toks1 = line1.strip().split()
                toks2 = line2.strip().split()
                if len(toks1) == 0:
                    contig += 1
                    continue
                elif toks1[0] == "Read":
                    if toks1[1] == "1:":
                        continue
                    if denom == 0:
                        contig_kmer_cov.append(0.0)
                        continue
                    contig_kmer_cov.append(total/denom)
                    total = 0
                    denom = 0
                else:
                    if toks2[1] == "1":
                        x = int(toks1[1])
                        if x < 300:
                            total += x
                            denom += 1
            if denom != 0:
                contig_kmer_cov.append(total/denom)
            else:
                contig_kmer_cov.append(0.0)
    with open(args.output+"/contig_kmer_cov.tsv",'w') as out:
        for (i,  cov) in enumerate(contig_kmer_cov):
            out.write("\t".join([str(i+1), str(cov)])+"\n")
            print("contig " + str(i+1) + " " + str(i+1) + " cov " + str(cov))
    return contig_kmer_cov
    


    
def het_kmers_FASTK():
    bc_trim = 0
    r1s = []
    r2s = []
    if args.kmer_data == "linked_reads":
        bc_trim = 23
        with open(args.linked_reads) as fofnin:
            for (index, line) in enumerate(fofnin):
                if index % 2 == 0:
                    r1s.append(line.strip().split()[0])
                else:
                    r2s.append(line.strip().split()[0])
    elif args.kmer_data == "short_reads":
        r1s = fofn_to_list(args.short_reads)
    elif args.kmer_data == "ccs_reads":
        r1s = fofn_to_list(args.ccs_reads)
    else:
        assert False, "kmer_data must be one of linked_reads, short_reads, or ccs_reads"

    cmds = []
    max_p = 0
    threads = args.threads

    mem = args.mem
    if len(r2s) > 0:
        mem = args.mem // 2
        threads = args.threads // 2

    
    if len(r2s) > 0:
        name = "fastk_spectrum_R1"
        cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), "-t"+str(args.min_kmer_count), 
        "-bc"+str(bc_trim), "-N"+args.output+"/"+name, "-M"+str(mem), "-T"+str(threads)] + r1s 
        cmds.append(cmd)
        name = "fastk_spectrum_R2"
        cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), "-t"+str(args.min_kmer_count), 
            "-N"+args.output+"/"+name, "-M"+str(mem),
            "-T"+str(threads)] + r2s
        cmds.append(cmd)
    else:
        name = "fastk_spectrum"
        cmd = [directory+"/FASTK/FastK", "-k"+str(args.kmer_size), "-t"+str(args.min_kmer_count), 
        "-bc"+str(bc_trim), "-N"+args.output+"/"+name, "-M"+str(mem), "-T"+str(threads)] + r1s 
        cmds.append(cmd)
    
    for cmd in cmds:
        print(" ".join(cmd))
    print(threads)
    print(mem)
    
    with ThreadPoolExecutor(max_workers = 2) as executor:
        procs = []
        for (index, cmd) in enumerate(cmds):
            procs.append(executor.submit(check_call, cmd, "spectra_"+str(index)+"_proc"))
        for proc in concurrent.futures.as_completed(procs):
            print("waiting for proc")


    if len(cmds) > 1:
        cmd = [directory+"/FASTK/Logex", "-T"+str(args.threads),  '-h'+str(args.min_kmer_count)+":1000",
            "'"+args.output+"/fastk_spectrum=(A|+B)'", 
            args.output+"/fastk_spectrum_R1", args.output+"/fastk_spectrum_R2"]
        print(" ".join(cmd))
        base_out_name = "fastk_spectrum"
        with open(args.output+"/"+base_out_name+".out",'w') as out:
            with open(args.output+"/"+base_out_name+".err", 'w') as err:
                subprocess.check_call(" ".join(cmd), shell = True, stdout = out, stderr = err)
        #check_call(" ".join(cmd), "fastk_spectrum", shell= True)
        
    # histogram
    cmd = [directory+"/FASTK/Histex", "-A", "-h"+str(args.min_kmer_count)+":256", args.output+"/fastk_spectrum"]
    print(" ".join(cmd))
    check_call(cmd, "histex")
    
    
    


def purge_dups(fn):
    cmd = [directory + "/purge_dups/bin/calcuts", fn]
    check_call(cmd, "calcuts")

def load_cutoffs():
    with open(args.output+"/calcuts.out") as cuts:
        for line in cuts:
            if line.startswith("["):
                continue
            else:
                toks = line.strip().split()
                minimum = int(toks[1])
                middle = int(toks[3])
                maximum = int(toks[5])
                maxhom = int(toks[5])
                return [minimum, middle, maximum, maxhom]

def het_kmer_molecules(cutoffs):
    cmd = [directory + "/het_snp_molecule_kmers.py", "--hic_reads", args.hic_reads,
        "--fasta", args.fasta, "--kmer_size", str(args.kmer_size), "--output", args.output,
        "-t", str(args.threads), "-m", str(args.mem), "--min_coverage", str(cutoffs[0]), 
        "--max_coverage", str(cutoffs[1]), "--max_total_coverage", str(cutoffs[1]*2),
        "--hom_modimizer", str(args.hom_modimizer)]
    if args.linked_reads:
        cmd.extend(["--txg_reads", args.linked_reads, "--whitelist", args.barcode_whitelist])
    if args.ccs_reads:
        cmd.extend(["--long_reads", args.ccs_reads])
    check_call(cmd, "het_kmer_molecules")

def het_kmer_molecules_FASTK(cutoffs):
    print("running phasemer once for het kmer grouped output")
    cmd = [directory + "/PHASE-MERS/Phasemer", "-h"+str(cutoffs[0])+':'+str(cutoffs[1]),
        "-m5.0", "-d"+str(cutoffs[1])+":"+str(cutoffs[2]), "-Ls", args.output+'/fastk_spectrum']
    print(" ".join(cmd))
    check_call(cmd, "haplex")
    print("running phasemer again for sorted output")
    cmd = [directory + "/PHASE-MERS/Phasemer", "-h"+str(cutoffs[0])+':'+str(cutoffs[1]),
        "-m5.0", "-N"+args.output+"/phasemer", "-d"+str(cutoffs[1])+":"+str(cutoffs[2]), args.output+'/fastk_spectrum']
    print(" ".join(cmd))
    check_call(cmd, "phasemer")
    
    cmds = []
    threads = args.threads // 2
    mem = args.mem // 2
    ccs_files = []
    with open(args.ccs_reads) as ccs:
        for line in ccs:
            ccs_files.append(line.strip())
    print("running fastk profiles on phasemer.U")
    cmd = [directory + "/FASTK/FastK", "-T"+str(threads), "-M"+str(mem), "-k"+str(args.kmer_size), 
        "-p:"+args.output+"/phasemer.U", "-N"+args.output+"/phasemap.U"] + ccs_files
    #check_call(cmd, "fastk_phasemer.U")
    cmds.append(cmd)
    print("running fastk profiles on phasemer.L")
    cmd = [directory + "/FASTK/FastK", "-T"+str(threads), "-M"+str(mem), "-k"+str(args.kmer_size), 
        "-p:"+args.output+"/phasemer.L", "-N"+args.output+"/phasemap.L"] + ccs_files
    print(" ".join(cmd))
    cmds.append(cmd)
    #check_call(cmd, "fastk_phasemer.L")


    with ThreadPoolExecutor(max_workers = 2) as executor:
        procs = []
        for (index, cmd) in enumerate(cmds):
            procs.append(executor.submit(check_call, cmd, "fastk_phasemer_"+str(index)))
        for proc in concurrent.futures.as_completed(procs):
            print("waiting for proc")

    print("running phasemap to combine phasemer.U and phasemer.L")
    cmd = [directory + "/PHASE-MERS/Phasemap",
        args.output+"/phasemap.U", args.output+"/phasemap.L"] + ccs_files
    print(" ".join(cmd))
    check_call(cmd, "ccs_phasemap")

    hic_files = []
    with open(args.hic_reads) as hic:
        for line in ccs:
            ccs_files.append(line.strip())


        cmd = [directory + "/FASTK/FastK", "-k"+str(args.kmer_size), "-t1", "-N"+args.output+"/het_kmers",
            "-M"+str(args.mem), "-T"+str(args.threads), args.output + "/het_kmers.fasta"]
        check_call(cmd, "het_kmers_fastk")
        # binary profiles ccs
        ccs_files = []
        with open(args.ccs) as ccs:
            for line in ccs:
                ccs_files.append(line.split())
        txg_files = []
        with open(args.linked_reads) as txg:
            for line in txg:
                txg_files.append(line)
        hic_files = []
        with open(args.hic_reads) as hic:
            for line in hic:
                hic_files.append(hic.split())
        
        cmds = []



    
    
    

def phasing():
    cmd = [directory + "/phasst_phase/target/release/phasst_phase", "--output", args.output, 
        "--het_kmers", args.output + "/het_kmers.tsv", "--hic_mols", args.output + "/hic.fofn",
        "--threads", str(args.threads), "--assembly_fasta", args.fasta, "--assembly_kmers", 
        args.output + "/fasta_kmers.bin"]
    if args.linked_reads:
        cmd.extend(["--linked_read_mols", args.output + "/txg.fofn"])
    if args.ccs_reads:
        cmd.extend(["--long_read_mols", args.output + "/ccs.fofn"])
    check_call(cmd, "phasing")

def new_phasing():
    cmd = [directory + "/phaser/target/release/phaser", "--output", args.output,
        "--het_kmers", args.output + "/het_kmers.tsv", "--hic_mols", args.output + "/hic.fofn",
        "--threads", str(args.threads), "--assembly_fasta", args.fasta, "--assembly_kmers", 
        args.output + "/fasta_kmers.bin", "--contig_kmer_depths", args.output+"/contig_kmer_cov.tsv"]
    if args.linked_reads:
        cmd.extend(["--linked_read_mols", args.output + "/txg.fofn"])
    if args.ccs_reads:
        cmd.extend(["--long_read_mols", args.output + "/ccs.fofn"])
    check_call(cmd, "phasing")


def phasst_scaff():
    if not os.path.exists(args.output + "/breaks_kmers"):
        subprocess.check_call(["mkdir", args.output + "/breaks_kmers"])
    if not os.path.exists(args.output + "/breaks_kmers/fasta_kmers.bin"):
        cmd = [directory + "/molecule_kmers/target/release/molecule_kmers", "-o", args.output + "/breaks_kmers", 
            "--paired_kmers", args.output + "/het_kmers.tsv", "--fasta", args.output + "/breaks.fa", 
            "--kmer_size", str(args.kmer_size), "--threads", "1"]
        check_call(cmd, "breaks_kmers/fasta_kmers")
    Fasta(args.output+"/breaks.fa") # make sure there is a .fai index
    cmd = [directory + "/phasst_scaff/target/release/phasst_scaff", "-o", args.output, "--het_kmers",
        args.output + "/het_kmers.tsv", "--linked_read_barcodes", args.output + "/txg.fofn",
        "--hic_mols", args.output + "/hic.fofn", "--assembly_fasta", args.output + "/breaks.fa",
        "--assembly_kmers", args.output + "/breaks_kmers/fasta_kmers.bin", "--phased_vcf", args.output + "/phasing_breaks.vcf"]
    check_call(cmd, "scaffolding")


def scaffolding():
    
    print("welcome to phasstools scaffolding")

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(args.output+"/contig_kmer_cov.tsv"):
        print("measuring unique kmer depth on contigs for sex/autosomal categorization")
        start = time.time()/60
        count_kmers_FASTK()
        end = time.time()/60
        print("unique kmer depth calculation took "+str(end-start)+"min")
    else:
        print("using previously calculated unique kmer depth for sex/autosomal categorization")


    hist_check = args.output + "/hist.tsv"
    if FASTK:
        hist_check = args.output + "/fastk_spectrum.hist"
    if not os.path.exists(hist_check):
        start = time.time()/60
        if FASTK:
            het_kmers_FASTK()
        else:
            het_kmers()
        end = time.time()/60
        print("het kmers took "+str(end-start)+"min")
    else:
        print("using previously generated kmer count data")
    if not os.path.exists(args.output + "/calcuts.out"):
        print("determining haploid coverage cutoffs")
        hist = args.output+"/hist.tsv"
        if FASTK:
            hist = args.output+"/histex.out"
        purge_dups(hist)
    else:
        print("using previously generaged kmer coverage cutoffs")
    if not os.path.exists(args.output + "/fasta_kmers.bin"):
        cutoffs = load_cutoffs()
        print("finding het kmers on read data")
        start = time.time()/60
        if FASTK:
            het_kmer_molecules_FASTK(cutoffs)
        else:
            het_kmer_molecules(cutoffs)
        end = time.time()/60
        print("het kmer molecules took "+str(end-start)+"min")
    else:
        print("using previously generated kmer molecule data")
    if not os.path.exists(args.output + "/breaks.fa"):
        print("phasing het kmers")
        start = time.time()/60
        new_phasing()
        end = time.time()/60
        print("phasing took "+str(end-start)+"min")
    else:
        print("using previously generated phasing")
    if not os.path.exists(args.output + "/scaff.tsv"):
        print("scaffolding with phased data")
        start = time.time()/60
        phasst_scaff()
        end = time.time()/60
        print("scaffolding took "+str(end-start)+"min")
    else:
        print("previous pipeline pointed at this directory is complete. If you want to rerun, specify a new output directory.")
    print("done.")



    




    

print(args)
if args.parser1:
    scaffolding()

