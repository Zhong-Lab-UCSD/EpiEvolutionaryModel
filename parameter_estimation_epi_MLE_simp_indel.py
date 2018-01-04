import argparse
import sys,copy
import math
from math import floor, log10, log, exp
from time import time
import scipy.optimize
from multiprocessing import *
from functools import partial

def ParseArg():
  p=argparse.ArgumentParser( description = "Parameter estimation for the evolutionary model using MLE. Model N and Model M" )
  p.add_argument("Input",type=str,help="Input file name")
  p.add_argument("-e", "--equil_file", type=str, help="The parameter file containing initial guesses for s, mu and k, and equilibrium probabilities estimated from the input data.")
  p.add_argument("-p","--process_num", type=int, default=6, help="Number of processes.")
  p.add_argument("--hypN",type=int, default=1, help="Use which hypothesis for parameter estimation. hypN=1: Model N; hypN=2: Model M.")
  p.add_argument("-n","--iter_non",action='store_true', help="If specified, the program will only report the first likelihood. It can be used to calculate the likelihood when parameters are known.")
  p.add_argument("-o","--output",type=str,help="Likelihood file name. If -n is specified, the program will output region name and the likelihood. Otherwise it will only output the likelihood in each iteration.")
  p.add_argument("-O","--out_allvec",type=str,help="Parameter file name. The program will output the tested parameters in each iteration.")
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args()  

class HomoRegion:
  '''
  class for homologous regions.
  '''
  def __init__(self):
    self.S1=[]
    self.S2=[]
    self.L=0
    self.name=""


def ReadInput(fin_name):
  Slist=[]
  s1_count=0
  s2_count=0
  flag=0
  with open(fin_name,"r") as fin:
    S=[]
    line=fin.readline().strip()
    if "@" not in line:
      print >> sys.stderr,"The input file is not FASTQ"
      exit()
    flag=1
    while True:
      line=fin.readline().strip()
      if len(line)==0:
        if s1_count>s2_count:
          Sobj.S2=S
          Slist.append(Sobj)
          s2_count+=1
        else:
          print >>sys.stderr, "The number of sequences are different!"
          exit()
        break
      if line=="+":
        i=0
        flag=0
        continue
      if "@" in line:
        if s1_count>s2_count:
          Sobj.S2=S
          Slist.append(Sobj)
          s2_count+=1
        else:
          Sobj=HomoRegion()
          Sobj.name=line[1:]
          Sobj.S1=S
          s1_count+=1   
        flag=1
        S=[]
        continue
      if flag==1:
        line=line.upper()
        S+=[(x,) for x in line]
      else:
        S=S[0:i]+[a+(b,) for a, b in zip(S[i:(i+len(line))],line)]+S[(i+len(line)):]
        i+=len(line)
#    if len(S1)!=len(S2):
#      print >>sys.stderr, "The number of sequences are different!"
#      exit()
  return Slist


def ReadParameters(f_name):
  '''
  Build the dictionaries of equilibrium probabilities on the linear and log scales. 
  Sample return: 
  {'A': 0.25, 1: [0.9, 0.1], 'C': 0.25, 'T': 0.25, 'G': 0.25} 
  {'A': -1.386, 1: [-0.105, -2.303], 'C': -1.386, 'T': -1.386, 'G': -1.386}
  '''
  n_epi=0
  x=[]
  equil_dict={}
  log_equil_dict={}
  with open(f_name,"r") as fin:
    x.append(float(fin.readline().strip()))
    x.append(float(fin.readline().strip()))
    for line in fin:
      line=line.strip().split("\t")
      if len(line)>1:
        if n_epi==0:
          for item in line:
            equil_dict[item.split(":")[0]]=float(item.split(":")[1])
            log_equil_dict[item.split(":")[0]]=log(float(item.split(":")[1]))
          n_epi+=1
        elif n_epi>0:
          equil_dict[n_epi]=[0.0,0.0]
          log_equil_dict[n_epi]=[0.0,0.0]
          for item in line:
            equil_dict[n_epi][int(item.split(":")[0])]=float(item.split(":")[1])
            log_equil_dict[n_epi][int(item.split(":")[0])]=log(float(item.split(":")[1]))
          n_epi+=1
      else:
        x.append(float(line[0]))
  return x, equil_dict, log_equil_dict

def Epi_equilibrium(n_epi):
  '''
  Products of the equilibrium probabilities of epigenetic marks
  The keys are combinations of 1 and 0.
  '''
  S_epi={}
  log_S_epi={}
  for i in xrange(pow(2,n_epi)):
    k=bin(i)[2:].zfill(n_epi)
    v=1.0
    lv=0.0
    for j in xrange(0,n_epi):
      v=v*equil_dict[j+1][int(k[j])]
      lv+=log_equil_dict[j+1][int(k[j])]
    S_epi[k]=v
    log_S_epi[k]=lv

  return S_epi, log_S_epi

def PA_equilibrium(S,lamb,mu):
  '''
  Compute the equilibrium of all ancestrol sequence.
  '''

  PA=log(Gamma(len(S), lamb, mu))

  for a in S:
    PA+=log_equil_dict[a[0]]
  return PA


def Gamma(n, lam, mu):
  '''
  Equilibrium probability of a sequence of length n
  '''
  return (1-lam/mu)*math.pow(lam/mu,n)

def Link_prob(prime, n, b, lam, mu):
  """
  Compute p', p'', p'''
  """
  if prime==1 and n==0:
    return mu*b
  elif prime==0 and n==1:
    return math.exp(-1*mu)*(1-lam*b)
  elif prime==1 and n==1:
    return (1-math.exp(-1*mu)-mu*b)*(1-lam*b)
  elif prime==2 and n==1:
    return (1-lam*b)

def Transition_f(s, base1, base2):
  e=math.exp(-s)
  if base1==base2:
    return e+equil_dict[base2]*(1-e)
  else:
    return equil_dict[base2]*(1-e)

def Transition_g(i, e1, e2, k):
  e=exp(-k)
  if e1==e2:
    return e+equil_dict[i][int(e2)]*(1-e)
  else:
    return equil_dict[i][int(e2)]*(1-e)


def Trans_matrix(n_epi, x):
  '''
  Transition probabilities of epigenetic marks. Keys: the i-th epigenetic mark. Dict of dict of dict. 
  {i:{0:{0:a,1:b},1:{0:c,1:d}}}
  '''
  trans_dic={}
  base=["A","C","G","T"]
  for b in base:
    trans_dic[b]={}
    for b1 in base:
      trans_dic[b][b1]=Transition_f(x[0], b, b1)
  for i in xrange(1,n_epi+1):
    trans_dic[i]={}
    for e1 in ['0','1']:
      trans_dic[i][e1]={}
      for e2 in ['0','1']:
        trans_dic[i][e1][e2]=Transition_g(i, e1, e2, x[1+i])
  return trans_dic

def Transition_g_prod(tuple1, tuple2, trans_dic):
  g=1.0
  if hypN==1 or tuple1[0]==tuple2[0]:
    for i in xrange(1,len(tuple1)):
      g=g*trans_dic[i][tuple1[i]][tuple2[i]]
    return g
  else:
    return S_epi["".join(tuple2[1:])]


def Log_sum(tmp):
  '''
  tmp=[logA, logB, logC, logD]
  return log(A+B+C-D)
  '''
 
  tm=max(tmp)
  sm=sum([exp(t-tm) for t in tmp[0:3]])
  sm=sm-exp(tmp[3]-tm)
  if sm>0:
    return tm+log(sm)
  else:
    return float("-Inf")

def Log_sum2(tmp):
  '''
  tmp=[logA, logB, logC]
  return log(A+B+C)
  '''
 
  tm=max(tmp)
  sm=sum([exp(t-tm) for t in tmp])
  if sm>0:
    return tm+log(sm)
  else:
    return float("-Inf")



def LocalPi(S1, trans_dic):
 # print >>sys.stderr, trans_dic
  epi_n=len(S1[0][1:])
  SS_epi={}
  log_SS_epi={}
  phi_list=[]
  m=len(S1)*1.0

  for i in xrange(1, epi_n+1):
    s=sum([1 for b in S1 if b[i]=='1'])
    phi_list.append([1-s/m, s/m])

  for i in xrange(pow(2,epi_n)):
    k=bin(i)[2:].zfill(epi_n)
    v=1.0
    for j in xrange(0,epi_n):
      v=v*phi_list[j][0]*trans_dic[j+1]['0'][k[j]] + phi_list[j][1]*trans_dic[j+1]['1'][k[j]]
    SS_epi[k]=v
    log_SS_epi[k]=log(v)

  return SS_epi, log_SS_epi


def Manhattan(S, X, trans_dic):
  ##put the argument to be mapped at first
  Na=float('-Inf')
  S1=S.S1
  S2=S.S2
  m=len(S1)
  n=len(S2)
  s=X[0]
  mu=X[1]
  kappa=X[2:]
  lamb=mu*(n+m)/(n+m+2)


  ##Initial guess
  #Sequence
  beta=(1-math.exp(lamb-mu))/(mu-lamb*math.exp(lamb-mu))
  link_p=[Link_prob(1,0,beta,lamb,mu),Link_prob(0,1,beta,lamb,mu),Link_prob(1,1,beta,lamb,mu),Link_prob(2,1,beta,lamb,mu)]
  if (1-exp(-mu)-mu*beta)<0:
    return float('Inf')

  log_link_p=[log(lp) for lp in link_p]
  log_lamb_mu=log(lamb/mu)
  lamb_beta=lamb*beta
  log_lamb_beta=log(lamb*beta)


  PA=PA_equilibrium(S1,lamb, mu)

  S1_epi, log_S1_epi=LocalPi(S1, trans_dic)
  S2_epi, log_S2_epi=LocalPi(S2, trans_dic)

  manh=[[Na]*(n+1) for i in range(2)]
  manh[0][0]=log_link_p[3]  

  manh[1][0]=log_link_p[3]+log_link_p[0]+log_S2_epi["".join(S1[0][1:])]   
  manh[0][1]=log_link_p[3]+log_lamb_beta+log_equil_dict[S2[0][0]]+log_S1_epi["".join(S2[0][1:])]   
  for i in xrange(2,n+1):
    manh[0][i]=manh[0][i-1]+log_lamb_beta+log_equil_dict[S2[i-1][0]]+log_S1_epi["".join(S2[i-1][1:])]
 
  for i in xrange(1,(m+1)):
    for j in xrange(1,(n+1)):

      g1=log( S_epi["".join(S1[i-1][1:])] * trans_dic[S1[i-1][0]][S2[j-1][0]] * Transition_g_prod(S1[i-1],S2[j-1], trans_dic) * link_p[1] + equil_dict[S2[j-1][0]] * S1_epi["".join(S2[j-1][1:])] * S2_epi["".join(S1[i-1][1:])] * link_p[2])

      g2=log_equil_dict[S2[j-1][0]]+log_S1_epi["".join(S2[j-1][1:])]+log_lamb_beta

      manh[1][j]=Log_sum([log_link_p[0]+manh[0][j]+log_S2_epi["".join(S1[i-1][1:])], g2+manh[1][j-1], g1+manh[0][j-1], g2+log_link_p[0]+manh[0][j-1]+log_S2_epi["".join(S1[i-1][1:])] ])


    manh[0]=manh[1]
    manh[1]=[Na]*(n+1)
    if i<m:
      manh[1][0]=manh[0][0]+log_link_p[0]+log_S2_epi["".join(S1[i][1:])] 


  S.L=manh[0][n]+PA
  return (manh[0][n]+PA, S)


def Manhattan_obj(x,S,fout):
 
  global iter_i
  iter_i+=1
  if iter_i%50==0:
    print >>sys.stderr, "Iteration:%d"%(iter_i)

  if min(x)<0:
    print >>sys.stderr, "minus p"
    return float("inf")
  Transition_dic=Trans_matrix(len(S[0].S1[0])-1, x)

  p=Pool(p_num)
  mp_queue=p.map(partial(Manhattan, X=x, trans_dic=Transition_dic),S)

  ll=sum([l[0] for l in mp_queue])
  S[:]=[l[1] for l in mp_queue]
  p.close()
  p.join()
  if iter_non:
    return
  if ll>0:
    print >>sys.stderr, "minus p"
    return float("inf")  
  print "\t".join([str(xx) for xx in x])
  print >>fout, str(-ll)
  return -ll

def Main():
  ##All of the global vairables cannot be changed 
  args=ParseArg()
  S=ReadInput(args.Input)

  global hypN
  hypN=args.hypN

  global p_num
  p_num=args.process_num

  ##equil_dic, log_equil_dict: 
  global equil_dict, log_equil_dict
  x, equil_dict, log_equil_dict=ReadParameters(args.equil_file)

  global S_epi, log_S_epi
  S_epi,log_S_epi=Epi_equilibrium(len(S[0].S1[0])-1)

  global iter_i
  iter_i=0

  global iter_non
  iter_non=args.iter_non

  t0=time()
  fout=open(args.output,"w")
  if args.iter_non:
    Manhattan_obj(x,S,fout)
    for s in S:
      print >>fout, "\t".join([s.name, str(s.L)])
  else:
    xopt, fopt, itera, funcalls, warnflag, allvecs = scipy.optimize.fmin(Manhattan_obj,x,(S,fout),disp=False, full_output=True, retall=True)
    t1=time()
    print >>sys.stderr, "Time:%f"%((t1-t0)/60)
    print >>sys.stderr, "Iterations: %d"%(itera)
    fout2=open(args.out_allvec,"w")
    print >>fout2, "\n".join(["\t".join([str(f) for f in x]) for x in allvecs])
  fout.close()



Main()







