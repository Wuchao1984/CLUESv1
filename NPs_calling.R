#!/usr/local/bin/R Rscript
args <- commandArgs(TRUE)

#Case_sorted_bed="Mikkelsen_H3K27me3_dupRemove.bed";
Case_sorted_bed= args[1];
Case_sorted_bed = as.character(Case_sorted_bed);

#Ctrl_sorted_bed="Mikkelsen_WCE_dupRemove.bed";
Ctrl_sorted_bed= args[2];
Ctrl_sorted_bed = as.character(Ctrl_sorted_bed);

#Species="mm9";
Species= args[3];
Species = as.character(Species);

#Mappability=0.8;
Mappability= args[4];
Mappability = as.numeric(Mappability);

#OutputPath="testNode1_H3K27me3";
OutputPath= args[5];
OutputPath = as.character(OutputPath);

#FragmentSize=200;
FragmentSize= args[6];
FragmentSize = as.numeric(FragmentSize);

#SigThreshold=0.05;
SigThreshold= args[7];
SigThreshold = as.numeric(SigThreshold);

#ncpu=10;
ncpu= args[8];
ncpu = as.numeric(ncpu);

#NPs_file="NPs_H3K27me3.txt";
NPs_file= args[9];
NPs_file = as.character(NPs_file);
###################################################################
# This function will call NPs from the input data.
# Last modified by Chao Wu on 2017-05-21

		source("CLUES_functions.R");

		# Get genome information
		genome_info=addGenomeInfo (Species, Mappability);
		genome_size=genome_info$size;
		Chr=genome_info$chr_info;

		OutPath=paste("mkdir ", OutputPath, sep="");
		system(OutPath);

		case_name <- paste (Chr[1], Case_sorted_bed, sep="");
		tmp=paste("cp ", Case_sorted_bed, " ./", OutputPath, "/",case_name , sep="");
		print (tmp);
		system(tmp);


		control_name <- paste (Chr[1], Ctrl_sorted_bed, sep="");
		tmp=paste("cp ", Ctrl_sorted_bed, " ./", OutputPath,"/",control_name , sep="");
		system(tmp);

		tmp=paste("./", OutputPath, sep="");
		setwd(tmp);

		#Count reads number in case and control samples
		len=paste (case_name,"_length.txt",sep="");
		xtmp=paste ("wc -l ",case_name," >",len, sep="");
		system(xtmp);
		ytmp=read.table(len, header=FALSE);
		case_len=ytmp[1];
		rm_len=paste("rm ",len,sep="");
		system(rm_len);

		len=paste (control_name,"_length.txt",sep="");
		xtmp=paste ("wc -l ",control_name," >",len, sep="");
		system(xtmp);
		ytmp=read.table(len, header=FALSE);
		control_len=ytmp[1];
		rm_len=paste("rm ",len,sep="");
		system(rm_len);

		Scale_factor=as.numeric(case_len/control_len);
		print ("Scale_factor");
		print (Scale_factor);
		Genome_density=as.numeric(case_len/genome_size);
		print ("Genome_density");
		print (Genome_density);

		

		#####################################################################################################################################
		# Reads shifiting module
		# Tune for shift parameter
		nlines=1234;
		case_name_nlines=paste(case_name, nlines,sep="");
		case_select <- paste ("cat ",case_name, " >", case_name_nlines, sep="");
		print (case_select);
		system(case_select);

		
		
		# separate reads in positive and negative strand
		for (i in 1:length(Chr[1:5]))
		{
			case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name_nlines, " > ",Chr[i],".bed",sep="");
			case_pos <- paste("awk '{OFS==\"\\t\"} {if($6==\"+\") print ($2) }' ", Chr[i],".bed > ",Chr[i],"_seq_case_pos.txt",sep=""); 
			case_neg <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3) }' ", Chr[i],".bed > ",Chr[i],"_seq_case_neg.txt",sep=""); 	
			system(case_tmp);
			system(case_pos);
			system(case_neg);
		
			rm_file <- paste ("rm ",Chr[i],".bed",sep="");
			system(rm_file);
		}
		rm_file <- paste ("rm ",case_name_nlines,sep="");
		system(rm_file);
	
		shift_x=as.numeric(FragmentSize); # initial shift paramter	
		shift_set=seq(round(shift_x*0.2), shift_x, by=10);

		shift_res=c();
		for (mmm in shift_set)
		{
			shift=mmm;
			print(mmm);
			all_bins=c();
			for (i in 1:length(Chr[1:5]))
			{
				case_pos<-paste (Chr[i],"_seq_case_pos.txt",sep="");
				case_neg<-paste (Chr[i],"_seq_case_neg.txt",sep="");	
				if(file.info(case_pos)[1]>0)
				{
					case_reads_pos=read.delim(case_pos,header=FALSE);
					case_reads_pos=data.matrix(case_reads_pos);		
					case_reads_pos=case_reads_pos+shift;
					case_reads_neg=read.delim(case_neg,header=FALSE);
					case_reads_neg=data.matrix(case_reads_neg);				
					case_reads_neg=case_reads_neg-shift;
					case_reads=c(case_reads_pos, case_reads_neg);	
					case_reads=sort(case_reads);
			
					temp=case_reads[2:length(case_reads)]-case_reads[1:(length(case_reads)-1)];
					all_bins=rbind(all_bins, temp);	
				}
			}

			all_bins=sort(all_bins);
			bin_res=sum(all_bins<=1);
			tmp=c(shift, bin_res);
			shift_set=rbind(shift_set, tmp);
		}

		a=max(which(shift_set[,2]==max(shift_set[,2])));
		Shift=shift_set[a,1];
		Shift=as.character(Shift);
		print (Shift);

		rm_file <- paste ("rm *_seq_case_pos.txt",sep=""); 
		system(rm_file);
		rm_file <- paste ("rm *_seq_case_neg.txt",sep=""); 
		system(rm_file);
		####################################################################################################################################################


		####################################################################################################################################################
		# NP calling module
		
		# Tune for BinSize_limit parameter
		# initial parameter set for parallel computing
		# The training chromosomes
		chr_test=Chr[1];
		for (i in Chr[2:5])
		{
			chr_test=paste(chr_test,"M",i,sep="");
		}		
		# Get the candidates of binsize_limit
		para=c(0.1,0.2,0.3, 0.4, 0.5, 0.6,0.7,0.8, 0.9, 1);
		para_set_all=c();
		for (i in para)
		{
			x <- paste(chr_test,case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, i, sep=",");
			para_set_all=c(para_set_all, x);
		}
		
		# parallel computing to get length distribution of NPs under different binsize_limit paramters
		library(parallel)
		cl <- makeCluster(getOption("cl.cores", ncpu))
		NPsDist <- parLapply(cl, para_set_all, multicore_BinSizeLimit)
		
		NPsDist_res=c();
		for (i in 1:length(NPsDist))
		{
			tmp=NPsDist[i];
			tmp=tmp[[1]]$NPs_lenDistr;
			NPsDist_res=rbind(NPsDist_res, tmp);
		}


		#write.table(NPsDist_res,"NPsLenDistr.txt",row.names=FALSE, col.names=FALSE);
		write.table(Shift,"ShiftParameter.txt",row.names=FALSE, col.names=FALSE);	
		
		# Identify the default BinSize_limit	
		tmp=NPsDist_res[,3];
		tmp1=NPsDist_res[as.numeric(tmp)>(2*as.numeric(Shift)),1];
		if (length(tmp1)>0) {BinSize_Limit=min(as.numeric(tmp1));} else BinSize_Limit=2*as.numeric(Shift);
		print ("BinSize_Limit");
		print (BinSize_Limit);
		
		# Call NPs, multicore running version
		# initial parameter set for parallel computing
		para_set_all=c();
		for (i in 1:length(Chr))
		{
			x <- paste(Chr[i],case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, BinSize_Limit, sep=",");
			para_set_all=c(para_set_all, x);
		}

		# Call NPs by chromosome
		# library(parallel)
		cl <- makeCluster(getOption("cl.cores", ncpu))
		NPs_chr <- parLapply(cl, para_set_all, multicore_NPs)
		NPs_all=c();
		for (i in 1:length(NPs_chr))
		{
			tmp=NPs_chr[i];
			tmp=tmp[[1]]$NPs;
			NPs_all=rbind(NPs_all, tmp);
		}
		qfdr=p.adjust(10^(-as.numeric(NPs_all[,8])),"BY");
		qfdr=-log10(qfdr);
		NPs_all[,10]=qfdr;
		NPs_all=NPs_all[as.numeric(NPs_all[,10])>-log10(SigThreshold),];

		write.table(NPs_all, file=NPs_file,sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);
		rm_case_name=paste("rm ",case_name,sep="");
		system(rm_case_name);
		rm_control_name=paste("rm ",control_name,sep="");
		system(rm_control_name);



