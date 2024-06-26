import argparse
import os
from datetime import datetime
from tqdm import tqdm
import openpyxl
import pandas as pd

parser = argparse.ArgumentParser(description = 'Script for 16s analysis')
parser.add_argument("-i", "--InPath", help = "Path to input files, directory is used. This directory could contain all kinds of stuff, but make sure that it has fastq or fastq.gz", 
                   nargs='?',
                   type = str)
parser.add_argument("-o", "--OutPath", help = "Directory to write report in. This pipeline makes another subdirectory inside it and all stuff goes there.",
                   nargs='?',
                   type = str)
parser.add_argument("-overlap", "--Overlap", help = "Number of nucleotides to use while merging, default to 7",
                   nargs='?',
                   type = str, 
                   default=7)
parser.add_argument("-d", "--Depth", help = "Taxonomy depth to use for analysis, possible arguments are R,D, P, C, O, F, G, default to G",
                   nargs='?',
                   type = str, 
                   default="G")
parser.add_argument("-c", "--Cutoff", help = "Minimum percent to use as cutoff, default to 1",
                   nargs='?',
                   type = str, 
                   default=1)
parser.add_argument('--pe', action='store_true', help = 'Use this flag to merge reads (e.g. when using illumina pe reads)')
args = parser.parse_args()

curDir = os.path.dirname(os.path.realpath(__file__))

dbPath = curDir + '/kraken2/silva/'
pearBin = curDir + '/pear-0.9.11-linux-x86_64/bin/pear'

def getNow():
    return str(datetime.now()).replace(' ', '_').replace(':', '').split('.')[0]

def merge(f, inpath, outpath, overlap):
    os.system(f'{pearBin} -v {overlap} -f {inpath}{f[0]} -r {inpath}{f[1]} -o {outpath}{f[0].split("_S")[0]}')

def releaseKraken(fil, inpath, outpath):
    print(f'kraken2 --db {dbPath} --report {outpath}{fil.split(".")[0]}.kreport --output {outpath}{fil.split(".")[0]} {inpath}{fil}')
    os.system(f'kraken2 --db {dbPath} --report {outpath}/{fil.split(".")[0]}.kreport --output {outpath}/{fil.split(".")[0]} {inpath}{fil}')

c = getNow()
root = args.OutPath+c
if not os.path.exists(root):
    os.system(f'mkdir {root}')
mergedPath = root + '/merged'
if not os.path.exists(mergedPath):
    os.system(f'mkdir {mergedPath}')
resultsPath = root + '/results'
if not os.path.exists(resultsPath):
    os.system(f'mkdir {resultsPath}')

if args.pe:
    fils = sorted([fil for fil in os.listdir(args.InPath) 
                   if fil.endswith('.fq.gz') 
                   or fil.endswith('.fastq.gz') 
                   or fil.endswith('.fq') 
                   or fil.endswith('.fastq')])
    filsDoubled = [(fils[e], fils[e+1]) for e, fil in enumerate(fils) if e%2 == 0]
    for fil in tqdm(filsDoubled):
        merge(fil, args.InPath, mergedPath, args.Overlap)

if args.pe:
    fils = [fil for fil in os.listdir(args.InPath) if fil.endswith('.assembled.fastq')]
    
    for fil in tqdm(fils):
        releaseKraken(fil, mergedPath, resultsPath)
else:
    fils = sorted([fil for fil in os.listdir(args.InPath) 
                   if fil.endswith('.fq.gz') 
                   or fil.endswith('.fastq.gz') 
                   or fil.endswith('.fq') 
                   or fil.endswith('.fastq')])
    for fil in tqdm(fils):
        releaseKraken(fil, args.InPath, resultsPath)

results = []
stat = []
general = {}
fils = [fil for fil in os.listdir(resultsPath) if fil.endswith('report')]
print(fils)
for fil in sorted(fils):
    with open(resultsPath+'/'+fil, 'r') as handle:
        sub = {}
        name = fil.split('.')[0]
        results.append([name, '', '', '', '', '', ''])
        for lin in handle:
            lins = [l.strip() for l in lin.strip().split('\t')]
            if float(lins[0]) >= args.Cutoff:
                if str(lins[3]) == args.Depth:
                    sub[lins[5]] = lins[0]
                res = ['', *lins[0:6]]
                results.append(res)
        general[fil.split('.')[0]] = sub

results_df = pd.DataFrame(results, columns = ['sample', 
                                              'fraction', 
                                              'reads for this TXID',
                                              'reads for this TXID with no further classification', 
                                              'tax depth', 
                                              'TXID',
                                              'name'])

results_df.to_excel(os.path.join(root, 'results.xlsx'))
