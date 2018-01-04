import sys

if len(sys.argv)==1:
  print >>sys.stderr, "Usage: python Extend_lift.py original.bed lift.bed lift_back.bed output1 output2"
  sys.exit()

def Overlap(bed1,bed2):
  if bed1[0]!=bed2[0]:
    return False
  if bed1[5]!=bed2[5]:
    return False
  if int(bed2[2])<int(bed1[1]) or int(bed1[2])<int(bed2[1]):
    return False
  return True

Origin={}
with open(sys.argv[1],"r") as fin:
  for line in fin:
    line=line.strip().split("\t")
    Origin[line[3]]=line

Lift={}
with open(sys.argv[2],"r") as fin:
  for line in fin:
    line=line.strip().split("\t")
    Lift[line[3]]=line

with open(sys.argv[3],"r") as fin, open(sys.argv[4], "w") as sp1, open(sys.argv[5],"w") as sp2:
  for line in fin:
    line=line.strip().split("\t")
    if line[3] in Lift and line[3] in Origin:
      if Overlap(line, Origin[line[3]]):
        o_len=(float(Origin[line[3]][2])-float(Origin[line[3]][1]))
        l_len=(float(line[2])-float(line[1]))
        print >>sp1, "\t".join(Origin[line[3]])
        if l_len/o_len<0.9:
          eleft_len=max(0, int(line[1])-int(Origin[line[3]][1]))
          eright_len=max(0, int(Origin[line[3]][2])-int(line[2]))
          if Lift[line[3]][5]==Origin[line[3]][5]:
            Lift[line[3]][1]=str( int(Lift[line[3]][1])-eleft_len )
            Lift[line[3]][2]=str( int(Lift[line[3]][2])+eright_len+1 )
          else:
            Lift[line[3]][1]=str( int(Lift[line[3]][1])-eright_len )
            Lift[line[3]][2]=str( int(Lift[line[3]][2])+eleft_len+1 )

        print >> sp2, "\t".join(Lift[line[3]])







