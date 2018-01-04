import sys, argparse
from math import exp, log
import numpy
import numpy.random

class Parameters():
  def __init__(self):
    self.s=0.0
    self.mu=0.0
    self.lamb=0.0
    self.kappa=[]
    self.equil_dict={}


def ParseArg():
  p=argparse.ArgumentParser( description = "Generate simulation data under the indel-independent hypothesis." )
  p.add_argument("input",type=str,help="Parameter file.")
  p.add_argument("--seq_len",type=int, default=500, help="Length of the ancestoral sequence")
  p.add_argument("--seq_num",type=int, default=1000,help="Number of the region pairs")
  p.add_argument("--hypN",type=int, default=1, help="Use which hypothesis to generate the data. hypN=1: mutation-independent; hypN=2: mutation-dependent.")
  p.add_argument("-r","--ran_seed",type=int,help="Set random seed. The program will generate a random seed if not specified.")
  p.add_argument("-o","--output",type=str, help="Output file name. ")
  p.add_argument("-O","--refout",type=str,help="If specified, the program will output the reference alignment (correct answer) at the same time")
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args() 

def ReadParameters(f_name, par,n):
  '''
  Read the parameter file. 
  Store s, m, k and global pi's. 
  '''
  n_epi=0
  with open(f_name,"r") as fin:
    par.s=float(fin.readline().strip())
    par.mu=float(fin.readline().strip())
    par.lamb=par.mu*(n)/(n+1)
    for line in fin:
      line=line.strip().split("\t")
      if len(line)>1:
        if n_epi==0:
          for item in line:
            par.equil_dict[item.split(":")[0]]=float(item.split(":")[1])
          n_epi+=1
        elif n_epi>0:
          par.equil_dict[str(n_epi)]=[0.0,0.0]
          for item in line:
            par.equil_dict[str(n_epi)][int(item.split(":")[0])]=float(item.split(":")[1])
          n_epi+=1
      else:
        par.kappa.append(float(line[0]))
  return


def Transition_f(base1, base2, par):
  e=exp(-par.s)
  if base1==base2:
    return e+par.equil_dict[base2]*(1-e)
  else:
    return par.equil_dict[base2]*(1-e)

def Transition_g(i, e1, e2, par):
  e=exp(-par.kappa[int(i)-1])
  if e1==e2:
    return e+par.equil_dict[i][int(e2)]*(1-e)
  else:
    return par.equil_dict[i][int(e2)]*(1-e)

def Immortal_p(lamb, beta):
  t=(1-lamb*beta)
  imm_p=[t]
  while True:
    t=t*lamb*beta
    if t<1e-323:
      break
    else:
      imm_p.append(t)
  return imm_p

def Normal_p(mu, lamb, beta):
  t=(exp(-mu)*(1-lamb*beta))
  norm_p=[t]
  while True:
    t=t*lamb*beta
    if t<1e-323:
      break
    else:
      norm_p.append(t)
  n=len(norm_p)
  norm_p.append(mu*beta)
  norm_p.append((1-exp(-mu)-mu*beta)*(1-lamb*beta))
  t=(1-exp(-mu)-mu*beta)*(1-lamb*beta)
  while True:
    t=t*lamb*beta
    if t<1e-323:
      break
    else:
      norm_p.append(t)
  return norm_p, n

class Evolved_seq():
  def __init__(self, epi_n):
    self.S=""
    self.S_epi=[""]*epi_n
    self.Sa=""
    self.Sa_epi=[""]*epi_n
    self.Sd=""
    self.Sd_epi=[""]*epi_n
    self.struct=""
    self.flag=1


def LocalPhi(b_list, trans_dic, epi_n):
  '''
  b_list: local phi list
  '''
 # print >>sys.stderr, trans_dic
  SS_epi={}
  log_SS_epi={}

  for epi in xrange(0,epi_n):
    phi_0 = (1-b_list[epi])*trans_dic[str(epi+1)+'0'][0] + b_list[epi]*trans_dic[str(epi+1)+'1'][0]
    phi_1 = (1-b_list[epi])*trans_dic[str(epi+1)+'0'][1] + b_list[epi]*trans_dic[str(epi+1)+'1'][1]
  SS_epi[str(epi+1)]=[phi_0, phi_1]
 # log_SS_epi[str(epi+1)]=[log(f) for f in SS_epi[str(epi+1)]]

#  print SS_epi

  return SS_epi

def Evolution(n, par, base, base_p, trans_dic, imm_p, norm_p, norm_cut):
  '''
  Evole the ancestral sequence and histone modification signals in silico. 
  n: length of the ancestral sequence
  '''
  epi_n=len(par.kappa)
  #Simulate the ancestoral sequence S
  seq_pair=Evolved_seq(epi_n)

  phi_list1=[]
  for epi in xrange(1,epi_n+1):
    phi1=numpy.random.beta(1.5,(1-par.equil_dict[str(epi)][1])/par.equil_dict[str(epi)][1]*1.5, 1)[0] 
    phi_list1.append(phi1)

  S1_epi=LocalPhi(phi_list1, trans_dic, epi_n)



  for i in xrange(n):
    seq_pair.S+=numpy.random.choice(base, p=base_p)

  #immortal
  j=numpy.random.choice(numpy.arange(len(imm_p)), p=imm_p)

  if j>0:
    for i in xrange(j):
      seq_pair.Sa+="-"
      seq_pair.Sd+=numpy.random.choice(base, p=base_p)
      seq_pair.struct+=" "

      if seq_pair.flag:
        for epi in xrange(1,epi_n+1):
          seq_pair.Sd_epi[epi-1]+=numpy.random.choice(['0','1'], p=S1_epi[str(epi)])
          seq_pair.Sa_epi[epi-1]+="-"

  #normal
  for si, sa in enumerate(seq_pair.S):
    j=numpy.random.choice(numpy.arange(len(norm_p)), p=norm_p)

    if j<norm_cut:
      #the original link is alive
      seq_pair.Sa+=sa
      tmp_base=numpy.random.choice(base, p=trans_dic[sa])
      seq_pair.Sd+=tmp_base
      equ_f=0
      if sa==tmp_base:
        seq_pair.struct+="|"
        equ_f=1
      else:
        seq_pair.struct+=" "

      if seq_pair.flag:
        for epi in xrange(1,epi_n+1):
          tmp_epi=numpy.random.choice(['0','1'], p=[1-phi_list1[epi-1], phi_list1[epi-1]] )
          seq_pair.Sa_epi[epi-1]+=tmp_epi
          if hypN==1 or equ_f==1:
            seq_pair.Sd_epi[epi-1]+=numpy.random.choice(['0','1'], p=trans_dic[str(epi)+tmp_epi])
          else:
            seq_pair.Sd_epi[epi-1]+=numpy.random.choice(['0','1'], p=par.equil_dict[str(epi)])
          

      for i in xrange(j):
        seq_pair.Sd+=numpy.random.choice(base, p=base_p)
        seq_pair.Sa+="-"
        seq_pair.struct+=" "
        if seq_pair.flag:
          for epi in xrange(1,epi_n+1):
            seq_pair.Sd_epi[epi-1]+=numpy.random.choice(['0','1'], p=S1_epi[str(epi)])
            seq_pair.Sa_epi[epi-1]+="-"

    else:
      seq_pair.Sd+="-"
      seq_pair.Sa+=sa
      seq_pair.struct+=" "
      if seq_pair.flag:
        for epi in xrange(1,epi_n+1):
          seq_pair.Sd_epi[epi-1]+="-"
          seq_pair.Sa_epi[epi-1]+="*"   #seq_pair.S_epi[epi-1][si]

      for i in xrange(j-norm_cut):
        seq_pair.Sd+=numpy.random.choice(base, p=base_p)
        seq_pair.Sa+="-"
        seq_pair.struct+=" "
        if seq_pair.flag:
          for epi in xrange(1,epi_n+1):
            seq_pair.Sd_epi[epi-1]+=numpy.random.choice(['0','1'], p=S1_epi[str(epi)])
            seq_pair.Sa_epi[epi-1]+="-"

  m1=sum([1 for i in seq_pair.Sa_epi[0] if i !="-" and i!="*"])
  n=sum([1 for i in seq_pair.Sd if i !="-"])

  phi_list2=[]
  for epi in xrange(1, epi_n+1):
    phi_list2.append(seq_pair.Sd_epi[epi-1].count("1")*1.0/n)


  S2_epi=LocalPhi(phi_list2, trans_dic, epi_n)


  for epi in xrange(1, epi_n+1):
    seq_pair.Sa_epi[epi-1]=list(seq_pair.Sa_epi[epi-1])
    for i in xrange(len(seq_pair.Sa_epi[epi-1])):
      if seq_pair.Sa_epi[epi-1][i]=="*":
        seq_pair.Sa_epi[epi-1][i]=numpy.random.choice(['0','1'],  p=S2_epi[str(epi)])

    seq_pair.Sa_epi[epi-1]="".join(seq_pair.Sa_epi[epi-1])
 
  m=sum([1 for i in seq_pair.Sa if i !="-"])

  phi_list1_flip=[]
  for epi in xrange(1,epi_n+1):  
    phi_list1_flip.append(seq_pair.Sa_epi[epi-1].count("1")*1.0/m)
 

  S1_epi_flip=LocalPhi(phi_list1_flip, trans_dic, epi_n)

  for epi in xrange(1, epi_n+1):
    seq_pair.Sd_epi[epi-1]=list(seq_pair.Sd_epi[epi-1])
    for i in xrange(len(seq_pair.Sd_epi[epi-1])):
      if seq_pair.Sa_epi[epi-1][i]=="-":
        seq_pair.Sd_epi[epi-1][i]=numpy.random.choice(['0','1'],  p=S1_epi_flip[str(epi)])

    seq_pair.Sd_epi[epi-1]="".join(seq_pair.Sd_epi[epi-1])

  

  return seq_pair


def Main():
  args=ParseArg()
  if args.ran_seed:
    ran_seed=args.ran_seed
    numpy.random.seed(ran_seed)
  else:
    ran_seed=numpy.random.random_integers(0,100000000)
    numpy.random.seed(ran_seed)
  print >>sys.stderr, "The random seed is: %d\n"%(ran_seed)
  
  par=Parameters()
  n=args.seq_len
  global hypN
  hypN=args.hypN

  ReadParameters(args.input, par, n)

  beta=(1-exp(par.lamb-par.mu))/(par.mu-par.lamb*exp(par.lamb-par.mu))

  base=["A","C","G","T"]
  base_p=[par.equil_dict[x] for x in base]

  #transition probabilities
  Transition_dic={}
  for b in base:
    Transition_dic[b]=[]
    for b1 in base:
      Transition_dic[b].append(Transition_f(b,b1, par))
  for i in xrange(1,len(par.kappa)+1):
    for e1 in ['0','1']:
      Transition_dic[str(i)+e1]=[]
      for e2 in ['0','1']:
        Transition_dic[str(i)+e1].append(Transition_g(str(i),e1,e2,par))

  p_immortal=Immortal_p(par.lamb, beta)

  #normal link
  p_normal, norm_cut=Normal_p(par.mu, par.lamb, beta)


  with open(args.output,"w") as fout, open(args.refout,"w") as fout2:
    for i in xrange(args.seq_num):
      j=0
      Seq_pair=Evolution(n, par, base, base_p, Transition_dic, p_immortal, p_normal, norm_cut)
 
      print >>fout, "@species1"
      print >>fout, Seq_pair.S
      if Seq_pair.flag:
        for seq in Seq_pair.Sa_epi:
          print >>fout, "+"
          print >>fout, seq.replace("-","")        
      print >>fout, "@species2"
      print >>fout, Seq_pair.Sd.replace("-","")
      if Seq_pair.flag:
        for seq in Seq_pair.Sd_epi:
          print >>fout, "+"
          print >>fout, seq.replace("-","")

      if args.refout:
        while j+100<len(Seq_pair.Sa):
          print >>fout2, "@Sequence No. %d"%(i)
          if Seq_pair.flag:
            for seq in Seq_pair.Sa_epi[::-1]:
              print >>fout2, seq[j:j+100]
          print >>fout2, Seq_pair.Sa[j:j+100]
          print >>fout2, Seq_pair.struct[j:j+100]
          print >>fout2, Seq_pair.Sd[j:j+100]
          if Seq_pair.flag:
            for seq in Seq_pair.Sd_epi:
              print >>fout2, seq[j:j+100]
          print >>fout2, ""
          j+=100

        if Seq_pair.flag:
          for seq in Seq_pair.Sa_epi[::-1]:
            print >>fout2, seq[j:j+100]
        print >>fout2, Seq_pair.Sa[j:j+100]
        print >>fout2, Seq_pair.struct[j:j+100]
        print >>fout2, Seq_pair.Sd[j:j+100]
        if Seq_pair.flag:
          for seq in Seq_pair.Sd_epi:
            print >>fout2, seq[j:j+100]
        print >>fout2, ""

      if i%50==0:
        print >>sys.stderr, "Sequence pairs:%d\r"%(i),


Main()






