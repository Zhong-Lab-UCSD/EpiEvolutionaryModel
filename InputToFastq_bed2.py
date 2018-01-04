import sys, argparse
import pyBigWig
from xplib.Annotation import Bed
from xplib import DBI
import os, string
from multiprocessing import *

rev_table=string.maketrans('ACGTacgtN', 'TGCATGCAN')

def ParseArg():
  p=argparse.ArgumentParser( description = "Generate input file in fastq format for the Epi-evolutionary model." )
  p.add_argument("q_region",type=str, nargs='+', help="Files containing coordinates of query regions in bed format. Note that the order should be the same as species.")
  p.add_argument("-s", "--species", nargs=2, type=str, default=["hg38","rheMac8"], help="A list of genome assembly. Note that the list should have the same order as the peak files")
  p.add_argument("-f","--fasta_path", type=str, default="/home/jialu/evm/Sequence", help="Path to the genome sequence files (.fa files).")
  p.add_argument("--histone", nargs="+",type=str, default=["H3K4me3"], help="Name of histone modifications. Note that the list should have the same order as histone.")
  p.add_argument("--bg", type=str, help="File name. The file contains paths to ChIP-Seq peak calling files (.bed files).")
#  p.add_argument("-t","--bg_cutoff", type=int, default=4, help="Cutoff for coverage")
  p.add_argument("--s_path",type=str, default="samtools", help="Path to samtools")
  p.add_argument("-p","--p_num",type=int, default=5,help="Number of processes to be used.")
  p.add_argument("-o","--output",type=str,help="Output file name. ")
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args() 

def overlap(bed1,bed2):
  '''
  Compute the overlap between two regions. 
  If the two regions overlap, return the length of the overlap, otherwise return False.
  '''
  if (bed1.stop>bed2.start) and (bed1.start<bed2.stop):
    return min(bed1.stop,bed2.stop)-max(bed1.start,bed2.start)+1
  else:
    return False

def ReadHistones(fhist_name):
  '''
  Read ChIP-Seq peak files.
  '''
  sp1=[]
  sp2=[]
  with open(fhist_name,"r") as fhist:
    fhist.readline()
    while True:
      line=fhist.readline().strip()
      if "@" in line:
        break
      line=line.split("\t")
      sp1.append(DBI.init(line[0],"bed"))
    while True:
      line=fhist.readline().strip()
      if line=="":
        break
      line=line.split("\t")
      sp2.append(DBI.init(line[0],"bed"))
 
  return sp1,sp2


def revcomp(seq):
  return seq.translate(rev_table)[::-1]

def fetchSeq(qr,fasta,s_path):
  ''' Fetch sequence from the genome sequence file.'''
 
  region = '%s:%d-%d'%(qr.chr,qr.start,qr.stop-1)
  seq="".join(os.popen(s_path+" faidx "+fasta+" "+region).read().split('\n')[1:])
  if qr.strand=="-":
    seq = revcomp(seq)
  return seq.upper()

def fetchHistModSeq(qr, dbi):
  ''' Convert histone modification measurement to sequences of 1's and 0's based on the peak calling results.'''
  #[a,b)

  pointer=qr.start
  hseq=""
  interval_list=sorted(list(dbi.query(qr)))
  for line in interval_list:
    if not overlap(line,qr): 
      continue

    bstart=line.start
    bend=line.stop
    if bstart<pointer:
      hseq+="1"*(min(bend,qr.stop)-pointer)
    elif bstart>pointer:
      hseq+="0"*(bstart-pointer)+"1"*(min(bend,qr.stop)-bstart)
    else:
      hseq+="1"*(min(bend,qr.stop)-bstart)
    pointer=min(bend,qr.stop)

  if pointer!=qr.stop:
    hseq+="0"*(qr.stop-pointer)
  if qr.strand=="-":
    hseq=hseq[::-1]

  return hseq

def Generate_output_str(bed_pair):
  '''
  Generate region pairs. Note that the program will discard regions with N's in their sequences.
  '''

  bed1=bed_pair[0]
  bed2=bed_pair[1]

  output_list=[]
  for s,sp,qrs in zip(args.species, [sp1,sp2], [bed1,bed2]):
    output_list.append("@"+s+"_"+"_".join(args.histone)+"_"+bed1.id)

    seq=fetchSeq(qrs,args.fasta_path.rstrip('/')+'/'+s+'.fa',args.s_path)
    if "N" in seq:
      output_list=[]
      return
    output_list.append(seq)
    for fsp_name in sp:
      output_list.append("+")
      hseq=fetchHistModSeq(qrs, fsp_name)
      if len(hseq)!=len(seq):
        print bed1.id
        print len(seq)
        print len(hseq)
        print >>sys.stderr,"Wrong hseq length!"
        exit(0)
      output_list.append(str(hseq))

  if len(output_list)!=0:
    return "\n".join(output_list)


def Main():
  global args
  args=ParseArg()

  print >>sys.stderr, "Indexing..."
  global sp1, sp2
  sp1, sp2=ReadHistones(args.bg)
  print >>sys.stderr, "Done indexing"
  if len(sp1)!=len(sp2) or len(sp1)!=len(args.histone):
    print >>sys.stderr, "Check the number of histone modifications!"
    exit(0)

  fout=open(args.output,"w")


  input_list=[]
  with open(args.q_region[0],"r") as fbed1, open(args.q_region[1],"r") as fbed2:
    while True:
      line1=fbed1.readline().strip().split("\t")
      line2=fbed2.readline().strip().split("\t")
      if line1[0]=="":
        break
      if line2[0]=="":
        break

      bed1=Bed(line1)
      bed2=Bed(line2)
      if bed1.chr!="chrX" and bed1.chr!="chrY":
        try:
          int(bed1.chr.lstrip("chr"))
        except:
          continue
      if bed2.chr!="chrX" and bed2.chr!="chrY":
        try:
          int(bed2.chr.lstrip("chr"))
        except:
          continue

      input_list.append((bed1, bed2))

  p=Pool(args.p_num)
  out_queue=p.map(Generate_output_str,input_list)

  with open(args.output,"w") as fout:
    for item in out_queue:
      if item:
        print >>fout, item

Main()





  
