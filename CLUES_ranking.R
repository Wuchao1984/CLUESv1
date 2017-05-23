#!/usr/local/bin/R Rscript
args <- commandArgs(TRUE)

InputERs= args[1];
InputERs = as.character(InputERs);

type= args[2];
type = as.character(type);

RankedERs= args[3];
RankedERs = as.character(RankedERs);

ERs=read.table(InputERs, header=FALSE);
print (type);

RankERs <- function (ERs, type, RankedERs)
{
	if (type=="Bymotif")
	{
		print ("Bymotif")
		Tmp=ERs[order(-as.numeric(ERs[,10])),];
		Tmp=cbind(Tmp,1:length(ERs[,10]));
		Tmp=Tmp[order(-as.numeric(Tmp[,7])),];
		Tmp=cbind(Tmp,1:length(ERs[,10]));
		TmpRank=as.numeric(Tmp[,11])+as.numeric(Tmp[,12]);
		Tmp=cbind(Tmp,TmpRank);
		Tmp=Tmp[order(as.numeric(Tmp[,13])),1:10];
		write.table(Tmp, RankedERs, row.names=FALSE, col.names=FALSE);	
	}

	if (type=="Bygene")
	{
		Tmp=ERs[order(-as.numeric(ERs[,4])),];
		write.table(Tmp, RankedERs, row.names=FALSE, col.names=FALSE);	
	}

	if (type=="Bydomain")
	{
		Tmp=as.numeric(ERs[,3])-as.numeric(ERs[,2])+1;
		ERs=ERs[order(-Tmp),];
		Tmp=cbind(ERs,1:length(ERs[,8]));
		Tmp=Tmp[order(-as.numeric(Tmp[,8])),];
		Tmp=cbind(Tmp,1:length(ERs[,8]));		
		TmpRank=as.numeric(Tmp[,10])+as.numeric(Tmp[,11]);		
		Tmp=cbind(Tmp,TmpRank);
		Tmp=Tmp[order(as.numeric(Tmp[,12])),1:9];
		write.table(Tmp, RankedERs, row.names=FALSE, col.names=FALSE);	
	}

}

RankERs (ERs, type, RankedERs)
