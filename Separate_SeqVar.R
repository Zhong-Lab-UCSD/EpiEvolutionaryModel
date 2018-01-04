library(argparse)

logsum <- function(x, a){
	A <- x[1]
	B <- x[2]
	b <- 1-a
	if( (a+A)>(b+B) ){
		part <- exp(b+B-a-A)
		return( a+A+log(1+part) )
	}
	else{
		part <- exp(A+a-b-B)
		return( b+B+log(1+part) )
	}
}

# argparse
parser <- ArgumentParser(description="Mixture model")
# specify options
parser$add_argument("Input",type="character", help="Input file name. Input file columns should be: regionid, (likelihood of) model I, model B, model N, model M. A header is required.")
#regionid, modelI, modelB, modelN, modelM
parser$add_argument("-o","--output",help="Output file name")

args <- parser$parse_args()


fin <- read.table(args$Input, header=TRUE, sep="\t", stringsAsFactors = F)

##a: log-ind=0; b: log-mut=0
Type_count <- as.data.frame(table(fin$group))
#print(Type_count)
Total <-  sum(Type_count[,2])
la <- log( (Type_count[Type_count=="MN",2] + Type_count[Type_count=="MM",2]) / Total )
lb <- log( (Type_count[Type_count=="MI",2] + Type_count[Type_count=="MN",2]) / Total )

print(la)
print(lb)
fin$MI <- apply(fin[,c(4,2)], 1, logsum, a=la)
fin$MD <- apply(fin[,c(5,3)], 1, logsum, a=la)
fin$II <- apply(fin[,c(4,5)], 1, logsum, a=lb)
fin$ID <- apply(fin[,c(2,3)], 1, logsum, a=lb)
fin$M_diff <- (fin$MI - fin$MD)/abs(fin$MI + fin$MD)
fin$I_diff <- (fin$II - fin$ID)/abs(fin$II + fin$ID)
fin <- fin[order(-fin$diff),]

fout <- fin[,c("peak_name", "MI", "MD", "II", "ID", "group", "diff", "M_diff", "I_diff")]
write.table(fout, file=args$output, sep="\t", quote=F, col.names=T,row.names=F)