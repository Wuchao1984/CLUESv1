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

Shift = args[6];
Shift = as.character(Shift);

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
# This function will tune StepWindowSize_limit parameter for SER caling.
# Last modified by Chao Wu on 2017-05-21

		source("CLUES_functions.R");
		NPs_all <- read.table(NPs_file, header=FALSE);

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

		
		write.table(NPs_all,NPs_file, sep="\t",col.name=FALSE,row.name=FALSE, quote=FALSE);
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
		
		#######################################################################################################
		# SER calling module

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
		rm_NPs_file=paste("rm ",NPs_file,sep="");
		system(rm_NPs_file);
		rm_case_name=paste("rm ",case_name,sep="");
		system(rm_case_name);
		rm_control_name=paste("rm ",control_name,sep="");
		system(rm_control_name);

