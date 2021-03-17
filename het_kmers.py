#!/usr/bin/env python

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description = "find heterozygous kmers from short read data")

parser.add_argument("-i","--input", required=True, help = "input files one per line, all same file type, will use extensions to determine type (fa, fq, fa.gz, fq.gz, bam)")
parser.add_argument("-o","--output", required=True, help = "prefix output")
parser.add_argument("-k","--kmer_size", required = False, type = int, default = 21, help = "kmer size")
parser.add_argument("-t","--threads", required=True, help="threads")
parser.add_argument("--tmp",required=False, default = "/tmp", help = "temp directory, default /tmp")
parser.add_argument("-m","--memory", required=False, type=int, default = 24, help="memory in GB, default 24")
args = parser.parse_args()

subprocess.check_call(['mkdir',args.output])

mypath = os.path.dirname(os.path.realpath(__file__))
cmd = ["singularity", "exec", mypath+"/kmc.sif", "kmc", "-k"+str(args.kmer_size)]
if args.input[-3:] == 'bam':
    cmd.append('-fbam')
cmd.extend(["-ci4", "-t"+str(args.threads), "-m"+str(args.memory), "-sm", "@"+args.input, args.output+"/"+args.output, args.tmp])
with open(args.output+"/kmc.out",'w') as out:
    with open(args.output+"/kmc.err",'w') as err:
        subprocess.check_call(cmd, stdout=out, stderr=err)
cmd = ["singularity", "exec",  mypath+"/kmc.sif", "kmc_tools", "transform", args.output+"/"+args.output, "dump","-s", args.output+"/kmer_counts.tsv"]
with open(args.output+"/dump.out",'w') as out:
    with open(args.output+"/dump.err",'w') as err:
        subprocess.check_call(cmd, stdout = out, stderr = err)

# potentially replace awk hist with kmc tools hist
with open(args.output+"/awk.sh", 'w') as out:
    out.write("awk 'NF{ count[ $2 ]++} END{ for ( name in count ) { print name \"\t\" count[ name ] };} ' "+args.output+"/kmer_counts.tsv > "+args.output+"/hist.tsv")
subprocess.check_call(["chmod","777",args.output+"/awk.sh"])

with open(args.output + "/hist.err",'w') as err:
    subprocess.check_call([args.output+"/awk.sh"], stderr = err, shell=True)



    
