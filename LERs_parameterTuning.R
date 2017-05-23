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

#FC_threshold=1.5
FC_threshold= args[8];
FC_threshold = as.numeric(FC_threshold);

#ncpu=10;
ncpu= args[9];
ncpu = as.numeric(ncpu);

#Input_file="NPs_H3K27me3.txt";
Input_file= args[10];
Input_file = as.character(Input_file);



#########################################################################
# This function will tune StepWindowSize_limit parameter for LER caling.
# Last modified by Chao Wu on 2017-05-21

		InputERs<- read.table(Input_file, header=FALSE);
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
		print ("FR_RE");
		print (FR_RE);
		write.table(FR_RE,"FR_RE.txt",row.names=FALSE, col.names=FALSE);
		
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
		rm_XXXInputERs <- paste ("rm XXXInputERs.txt");
		system(rm_XXXInputERs);
		rm_case_name=paste("rm ",case_name,sep="");
		system(rm_case_name);
		rm_control_name=paste("rm ",control_name,sep="");
		system(rm_control_name);


