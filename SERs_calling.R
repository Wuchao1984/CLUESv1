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

Shift= args[6];
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

#SERs_file="SERs_H3K27me3.txt";
SERs_file= args[10];
SERs_file = as.character(SERs_file);

StepWindowSize_limit = args[11];
StepWindowSize_limit = as.character(StepWindowSize_limit);
###################################################################
# This function will call SER.
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
					
		# detect SERs with the StepWindowSize_limit, multicore running version				
		# initial parameter set for parallel computing
		para_set_all=c();
		for (i in 1:length(Chr))
		{
			x <- paste(Chr[i],case_name, control_name, Shift, Species, Genome_density, Scale_factor, SigThreshold, StepWindowSize_limit, NPs_file, sep=",");
			para_set_all=c(para_set_all, x);
		}

		library(parallel)
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
		rm_NPs_file=paste("rm ",NPs_file,sep="");
		system(rm_NPs_file);
		rm_case_name=paste("rm ",case_name,sep="");
		system(rm_case_name);
		rm_control_name=paste("rm ",control_name,sep="");
		system(rm_control_name);

