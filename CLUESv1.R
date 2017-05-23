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

#FR_D_threshold=0.01;
FR_D_threshold= args[9];
FR_D_threshold = as.numeric(FR_D_threshold);

#FC_threshold=1.5
FC_threshold= args[10];
FC_threshold = as.numeric(FC_threshold);

#NPs_file="NPs_H3K27me3.txt";
NPs_file= args[11];
NPs_file = as.character(NPs_file);

#SERs_file="SERs_H3K27me3.txt";
SERs_file= args[12];
SERs_file = as.character(SERs_file);

#LERs_file="LERs_H3K27me3.txt";
LERs_file= args[13];
LERs_file = as.character(LERs_file);


# This function will estimate the shift parameter from the input data.
# Last modified by Chao Wu on 2017-05-17

## If you want to source() a bunch of files, something like
## the following may be useful:
 sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
       if(trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }
#path="/public/home/chaowu/Ongoing_Project/clus_20151117/CLUES_package20170515/CLUES_v1";
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
		
		
		###################################################################################################
		# Deci_SER module
		# Calculate FR-D of NPs
		
		# Estimate fragment rate by distance(FR-D) of the NPs
		NPs_len=as.numeric(NPs_all[,3])-as.numeric(NPs_all[,2])+1;
		LeftNPs=NPs_len[1:(length(NPs_len)-1)];
		RightNPs=NPs_len[2:length(NPs_len)];
		max_Len=cbind(LeftNPs,RightNPs);
		max_Len=1/2*apply(max_Len,1,max);
		
		Between_NPs=as.numeric(NPs_all[2:length(NPs_all[,2]),2])-as.numeric(NPs_all[1:(length(NPs_all[,2])-1),3])+1;
		Between_NPs[Between_NPs<0]=max(max_Len);
		tmp=Between_NPs<max_Len;
		FR_NPs=which(tmp==1);
		FRD_NPs=length(FR_NPs)/(length(max_Len)+1);

		print ("FRD_NPs");
		print (FRD_NPs);
		write.table(FRD_NPs,"FRD_NPs.txt",row.names=FALSE, col.names=FALSE);

		if (FRD_NPs>FR_D_threshold) 
		{
			SER_label=1;
		} else {
		  SER_label=0;
		}
		
		#######################################################################################################
		# SER calling module
		if (SER_label==1)
		{

			# Tune StepWindowSize_limit parameter for SER calling
			# Get the StepWindowSize_limit candidates
			Percentile_set=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1);
			NPs_number=length(NPs_all[,1]);
      Shuffled_NPs=sort(sample(genome_size, NPs_number));
			Shuffled_dist=Shuffled_NPs[2:length(Shuffled_NPs)]-Shuffled_NPs[1:(length(Shuffled_NPs)-1)];
			Shuffled_dist=sort(Shuffled_dist);

			# Initial parameter set for parallel computing of testing StepWindowSize_limit
			chr_test=Chr[1];
			for (i in Chr[2:5])
			{
				chr_test=paste(chr_test,"M",i,sep="");
			}
			para_set_all=c();
			for (i in Percentile_set)
			{
				StepWindowSize_limit=Shuffled_dist[max(length(Shuffled_dist)*i,1)];
				#if (StepWindowSize_limit<2*as.numeric(Shift)) next;
				#if (StepWindowSize_limit>10000) next;
				x <- paste(chr_test,case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, StepWindowSize_limit, NPs_file,i, sep=",");
				para_set_all=c(para_set_all, x);
			}
			library(parallel)
			cl <- makeCluster(getOption("cl.cores", ncpu));
			StepWindowSize_res <- parLapply(cl, para_set_all,  multicore_StepWindowSize_SERs)
			StepWindowSize_set=c();
			for (i in 1:length(StepWindowSize_res))
			{
				tmp=StepWindowSize_res[i];
				tmp=tmp[[1]]$StepWindowSize_stat;
				if (sum(tmp)>0)
				{
					StepWindowSize_set=rbind(StepWindowSize_set, tmp);
				}
			}
			tmp=min(which(StepWindowSize_set[,2]==max(StepWindowSize_set[,2])));
			if (tmp>1) StepWindowSize_set[1:(tmp-1),]=c(0,0,1);		
			write.table(StepWindowSize_set[,c(1,3)], file="StepWindowSize_limit_SERs.txt",sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);
						
			if (length(which(StepWindowSize_set[,3]<FR_D_threshold))<1|min(StepWindowSize_set[,3])==1)
			{				
				print ("no SERs are detected under the FR_D");
				InputERs=NPs_all;							
			} else {
				if (length(which(StepWindowSize_set[,3]<FR_D_threshold))>0)
				{
					tmp=	which(StepWindowSize_set[,3]<FR_D_threshold);
					StepWindowSize_limit=min(StepWindowSize_set[tmp,1]);
				}
				if (length(which(StepWindowSize_set[,3]<FR_D_threshold))<1)
				{
					tmp=	StepWindowSize_set[which(StepWindowSize_set[,3]==min(StepWindowSize_set[,3])),1];
					StepWindowSize_limit=min(tmp);				
				}
				
				
							
				# detect SERs with the StepWindowSize_limit, multicore running version				
				# initial parameter set for parallel computing
				para_set_all=c();
				for (i in 1:length(Chr))
				{
					x <- paste(Chr[i],case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, StepWindowSize_limit, NPs_file, sep=",");
					para_set_all=c(para_set_all, x);
				}

				#library(parallel)
				cl <- makeCluster(getOption("cl.cores", ncpu))
				SERs_res <- parLapply(cl, para_set_all,  multicore_SERs)

				SERs_set=c();
				for (i in 1:length(SERs_res))
				{
					tmp=SERs_res[i];
					tmp=tmp[[1]]$SERs;
					SERs_set=rbind(SERs_set, tmp);
				}

				qfdr=p.adjust(10^(-as.numeric(SERs_set[,8])),"BY");
				qfdr=-log10(qfdr);
				SERs_set[,10]=qfdr;
				SERs_set=SERs_set[as.numeric(SERs_set[,10])>-log10(SigThreshold),];
				InputERs=SERs_set;	
				write.table(SERs_set[,c(1:6,8:10)], file=SERs_file,sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);
			}
		}	else {
			print ("no SERs are detected");
			InputERs=NPs_all;
		}
		write.table(InputERs,"XXXInputERs.txt",sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);
		
		###################################################################################################
		# Deci_LER module
		# Calculate FR-RE of InputERs
		FC_set=c();
		for (i in 1:length(Chr))
		{
			print (Chr[i]);
			case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",Chr[i],".bed",sep="");
			case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", Chr[i],".bed > ",Chr[i],"_seq_case.txt",sep="");
			system(case_tmp);
			system(case_tmp1);
			rm_file <- paste ("rm ",Chr[i],".bed",sep="");
			system(rm_file);
			case <- paste (Chr[i],"_seq_case.txt",sep="");
			rm_case <- paste ("rm ",Chr[i],"_seq_case.txt",sep="");
	
			control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr[i],"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",Chr[i],".bed",sep="");
			control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", Chr[i],".bed > ",Chr[i],"_control.txt",sep=""); 
			system(control_tmp);
			system(control_tmp1);
			system(rm_file);
			rm_control <- paste ("rm ",Chr[i],"_control.txt",sep="");
			control<-paste (Chr[i],"_control.txt",sep="");
			if(file.info(case)[1]>0&&file.info(control)[1]>0)
			{
				seq_case=read.delim(case,header=FALSE);
				seq_input=read.delim(control,header=FALSE);
				seq_case=data.matrix(seq_case);
				seq_case=sort(seq_case);
				ab=which(seq_case<1);
				seq_case[ab]=1;

				seq_input=data.matrix(seq_input);
				seq_input=sort(seq_input);
				ab=which(seq_input<1);
				seq_input[ab]=1;
				x=which(InputERs[,1]==Chr[i]);
				if (length(x)<3) next;
				InputERs_chr=cbind(as.numeric(InputERs[x,2]),as.numeric(InputERs[x,3]));
				FC_res<-calculateFR_RE (seq_case, seq_input, Genome_density, Scale_factor, InputERs_chr);
			}
			system(rm_case);
			system(rm_control);
			FC_set=c(FC_set,FC_res$FC);
		}
		FR_RE=sum(FC_set>FC_threshold)/length(FC_set);
		if (FR_RE<0.01)
		{
			quit("yes");
		}
		
		#######################################################################################################
		# LER calling module				
		# Tune StepWindowSize_limit parameter for LER calling
		# Get the StepWindowSize_limit candidates and initial parameter set for parallel computing

		Percentile_set=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3, 0.4, 0.5);
		InputERs_number=length(InputERs[,1]);
    Shuffled_InputERs=sort(sample(genome_size, InputERs_number));
		Shuffled_dist=Shuffled_InputERs[2:length(Shuffled_InputERs)]-Shuffled_InputERs[1:(length(Shuffled_InputERs)-1)];
		Shuffled_dist=sort(Shuffled_dist);
		
		chr_test=Chr[1];
		for (i in Chr[2:5])
		{
			chr_test=paste(chr_test,"M",i,sep="");
		}
		para_set_all=c();
		for (i in Percentile_set)
		{
			StepWindowSize_limit=Shuffled_dist[max(length(Shuffled_dist)*i,1)];
			x <- paste(chr_test,case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, StepWindowSize_limit, "XXXInputERs.txt",i, sep=",");
			para_set_all=c(para_set_all, x);
		}

		library(parallel);
		cl <- makeCluster(getOption("cl.cores", ncpu));
		StepWindowSize_LERs_res <- parLapply(cl, para_set_all,  multicore_StepWindowSize_LERs);
		StepWindowSize_LERs_set=c();
		for (i in 1:length(StepWindowSize_LERs_res))
		{
			tmp=StepWindowSize_LERs_res[i];
			tmp=tmp[[1]];
			StepWindowSize_LERs_set=rbind(StepWindowSize_LERs_set, tmp);
		}
		tmp=min(which(StepWindowSize_LERs_set[,2]==max(StepWindowSize_LERs_set[,2])));
		if (tmp>1) StepWindowSize_LERs_set[1:(tmp-1),]=c(0,0,0,0,0,0);
		write.table(StepWindowSize_LERs_set[,c(1,4:6)], file="StepWindowSize_limit_LERs.txt",sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);
		
		if (length(which(StepWindowSize_LERs_set[,4]>FC_threshold))<1)
		{
			print ("no LERs is detected");
			rm_XXXInputERs <- paste ("rm XXXInputERs.txt");
			system(rm_XXXInputERs);
			quit("yes");
		} else {
			tmp=StepWindowSize_LERs_set[StepWindowSize_LERs_set[,4]>FC_threshold,1];
			LERs_StepWindowSize_limit=max(tmp);
		}

	
		# Calling LERs	
		para_set_all=c();	
		for (i in 1:length(Chr))
	  {
			x <- paste(Chr[i],case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, LERs_StepWindowSize_limit, "XXXInputERs.txt", sep=",");
			para_set_all=c(para_set_all, x);
		}

		library(parallel)
		cl <- makeCluster(getOption("cl.cores", ncpu));
		LERs_res <- parLapply(cl, para_set_all,  multicore_LERs)
		LERs_set=c();
		for (i in 1:length(LERs_res))
		{
			tmp=LERs_res[i];
			tmp=tmp[[1]]$LERs;
			LERs_set=rbind(LERs_set, tmp);
		}

		qfdr=p.adjust(10^(-as.numeric(LERs_set[,8])),"BY");
		qfdr=-log10(qfdr);
		LERs_set[,10]=qfdr;
		LERs_set=LERs_set[as.numeric(LERs_set[,10])>-log10(SigThreshold),];

		write.table(LERs_set[,c(1:6,8:10)], file=LERs_file,sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);		
		rm_XXXInputERs <- paste ("rm XXXInputERs.txt");
		system(rm_XXXInputERs);

		rm_case_name=paste("rm ",case_name,sep="");
		system(rm_case_name);
		rm_control_name=paste("rm ",control_name,sep="");
		system(rm_control_name);












