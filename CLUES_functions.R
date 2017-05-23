scoreFinalWindow <- function (ReadsLocation_case, ReadsLocation_control, InitialWindowSize, Steps, Genome_density, Scale_factor, StepWindowSize_limit ,BinSize_limit)
# Calculate p value of each final window in the chromosome.
# The size of the step window are limited by max allowed length parameter under NP calling ...
# and are limited by the length parameter of SER and LER calling of CLUES methodlogy.
# Last modificated by Chao Wu on 2017-05-16

{

		if (length(ReadsLocation_case)==length(ReadsLocation_control))
		{
			if (sum(abs(ReadsLocation_case-ReadsLocation_control))==0) {WithCtrl=0;} else  WithCtrl=1;		
		} else WithCtrl=1;

		Input=ReadsLocation_control;

		# Estimate reads density in each postion in the control sample 
		Win_len=10; # 10bp, the reads density is estimated by every 10 bp		
		Chr_len=max(c(ReadsLocation_case, ReadsLocation_control));
		Chr_len=floor(Chr_len/Win_len)+1;
		ReadsDensity_ctrl=matrix(0, nrow=Chr_len, ncol=1);
		for (i in Input) ## estimate reads density in each position
		{
			Position=floor(i/Win_len)+1;
			ReadsDensity_ctrl[Position]=ReadsDensity_ctrl[Position]+1;
		}
		ReadsDensity_ctrl=ReadsDensity_ctrl/Win_len;

		BinSet_case=ReadsLocation_case[2:length(ReadsLocation_case)]-ReadsLocation_case[1:(length(ReadsLocation_case)-1)]; # bin chromosome with reads in case samples

		# initiate p value matrix of the final windows of the chromosome
		if (Steps>0) 
		{
			xx=seq(1,(length(ReadsLocation_case)-1),by=Steps);		
		}
		if (Steps==0)
		{
			xx=1:(length(ReadsLocation_case)-1);
		}
		FinalWindow_p=matrix(0,ncol=3,nrow=length(xx));
		FinalWindow_p_tmp=matrix(10,ncol=1,nrow=length(xx)*InitialWindowSize);

		kkk=1; # index for FinalWindow_p matrix
		kkktmp=1; # index for FinalWindow_p_tmp matrix
		for (ii in xx)
		{
			if (ii%%10000<1&&ii>10000) {print(ii);}
			if (ii>((length(ReadsLocation_case)-1)-InitialWindowSize)) {kkk=kkk+1;next;}
			
			Initial_size=1+InitialWindowSize;
			Initial_window=matrix(0,nrow=Initial_size, ncol=1);
			Len=0:(InitialWindowSize);
			Initial_window=ii+Len;
		
			# Filter step windows with BinSize_limit parameter for NP calling
			StepWindow=c(); StepWindow=BinSet_case[Initial_window];
			StepWindow_filter=which(StepWindow>=BinSize_limit);
			if (length(StepWindow_filter)==0) 
			{
				StepWindow=Initial_window;
			}	else {
				StepWindow_filter=min(StepWindow_filter);
				if (StepWindow_filter<4) {kkk=kkk+1; next;} 
				StepWindow=Initial_window[1:(StepWindow_filter-1)];
			}

			# Set the sequence of step windows
			len_neighbor=length(StepWindow)-1;
			if (len_neighbor<1000&&len_neighbor>=100) { sa=seq(11,20,1);sb=seq(21,100,5);sc=seq(101, len_neighbor+1, 20);bin_steps=c(sa,sb,sc);}
			if (len_neighbor<100&&len_neighbor>=20) { sa=seq(3,20,1);sb=seq(21,len_neighbor+1,5);bin_steps=c(sa,sb);}
			if (len_neighbor<20) {bin_steps=seq(3, len_neighbor+1, by=1);}

			bin_steps=bin_steps[order(-bin_steps)];
			bin_steps=c((len_neighbor+1),bin_steps);
			bin_steps=unique(bin_steps);
			p_list=matrix(10,ncol=2, nrow=len_neighbor+1);
		
			# calculate p value for each final window
			for (j in bin_steps)
			{
				OS=c();
				OS=StepWindow[1:j];
				reads=j+1;
				start=min(OS);
				start=ReadsLocation_case[start];
				end=max(OS);
				end=ReadsLocation_case[end+1];
				distance=end-start+1;
				
				if (distance>=StepWindowSize_limit) next; # Filter step windows with StepWindowSize_limit parameter for SER and LER calling
	
				local_density=0;
				if (WithCtrl>0)
				{			
					# Estimate reads density in control sample in the region
					down=floor(start/Win_len)+1;
					up=floor(end/Win_len)+1;
					if (down<1) down=1;
					if (up<1) up=1;	
					local_density=sum(ReadsDensity_ctrl[down:up])/(up-down+1);
				}
				
				# get p value for each step window
				p=ppois((reads-1),max(Genome_density,local_density*Scale_factor)*distance,lower=FALSE);
				p_list[j,]=c(j,p);
				
				FinalWindow_p_tmp[kkktmp]=p;
				kkktmp=kkktmp+1;
				if (p<0.001/length(ReadsLocation_case)) break;
			}
		
			p_min=min(p_list[,2]);
			p_min_index=p_list[max(which(p_list[,2]==p_min)),1];
			FinalWindow_p[kkk,]<- c(ii, p_min, p_min_index);
			kkk=kkk+1;	
		}
		FinalWindow_p=FinalWindow_p[FinalWindow_p[,1]>0,];
		FinalWindow_p=FinalWindow_p[FinalWindow_p[,2]<2,];
		FinalWindow_p_tmp=FinalWindow_p_tmp[FinalWindow_p_tmp<2];
		FinalWindow_p<- list(FinalWindow_p=FinalWindow_p, FinalWindow_p_tmp=FinalWindow_p_tmp);
		return (FinalWindow_p);
}



getSERs <- function (SigWindow, ReadsLocation_case, ReadsLocation_control, Genome_density, Scale_factor,Shift, NPs)
# getSERs will merge significant overlapped final window to get enriched regions (ERs) ...
# then it keeps the shortest ERs covering NPs and calculate p value of the ERs
# Last modified by Chao Wu on 2017-05-18

{

		Shift=as.numeric(Shift);
		if (length(ReadsLocation_case)==length(ReadsLocation_control))
		{
			if (sum(abs(ReadsLocation_case-ReadsLocation_control))==0) {WithCtrl=0;} else  WithCtrl=1;		
		} else WithCtrl=1;		
		
		if (sum(SigWindow)==0)
		{
			SERs_res <- list (SERs=c(0,0,0,0,0,0,0,0,0));
			return(SERs_res);
		}		
		
		a=length(SigWindow)/4;
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}

		SigWindow <- SigWindow[order(SigWindow[,1],SigWindow[,2]),];

		a=length(SigWindow)/4;
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}
		
		Tags=ReadsLocation_case;
		Input=ReadsLocation_control;

		# Build reads density background in genome in case samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_case_10=matrix(0, nrow=Len, ncol=1);
		for (i in Tags)
		{
			mn=floor(i/10)+1;
			RD_case_10[mn]=RD_case_10[mn]+1;
		}
		RD_case_10=RD_case_10/10;

		# Build reads density background in genome in control samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_ctrl_10=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10)+1;
			RD_ctrl_10[mn]=RD_ctrl_10[mn]+1;
		}
		RD_ctrl_10=RD_ctrl_10/10;
		
		# 1000 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/1000)+1;
		RD_ctrl_1000=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/1000)+1;
			RD_ctrl_1000[mn]=RD_ctrl_1000[mn]+1;
		}
		RD_ctrl_1000=RD_ctrl_1000/1000;

		# 10000bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10000)+1;
		RD_ctrl_10000=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10000)+1;
			RD_ctrl_10000[mn]=RD_ctrl_10000[mn]+1;
		}
		RD_ctrl_10000=RD_ctrl_10000/10000;

		# merge significant final windows to get ERs
		mm=sum(SigWindow[,2]-SigWindow[,1]+1);
		SigWindow_index=matrix(0, nrow=mm, ncol=1);
		mm=length(SigWindow[,1]);
		m=1;
		for (i in 1:mm)
		{
			x=SigWindow[i,2]-SigWindow[i,1]+1;
			SigWindow_index[m:(m+x-1)]=SigWindow[i,1]:SigWindow[i,2];
			m=m+x;
		}		
		SigWindow_index=unique(SigWindow_index);
		SigWindow_index=sort(SigWindow_index);
		SigWindow=matrix(0,nrow=length(SigWindow_index), ncol=4);
		m=SigWindow_index[1];
		for (i in 2:length(SigWindow_index))
		{
			if (SigWindow_index[i]-SigWindow_index[i-1]>1) {tmp=c(m,SigWindow_index[i-1],0,0); SigWindow[(i-1),]=tmp; m=SigWindow_index[i];}
		}
		tmp=c(m,SigWindow_index[i],0,0); SigWindow[i,]=tmp; m=SigWindow_index[i];
		SigWindow=SigWindow[which(SigWindow[,1]>0),];

		a=length(SigWindow)/4;
		if (length(a)==0)
		{
			SERs_res <- list (SERs=c(0,0,0,0,0,0,0,0,0));
			return(SERs_res);
		}				
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}
		SigWindow[,3]=SigWindow[,2]-SigWindow[,1]+2;
		ERs=SigWindow;
		ERs[,1]=Tags[ERs[,1]];
		ERs[,2]=Tags[ERs[,2]+1];

		# count NPs in ERs and filter the ERs without NPs
		m=length(ERs)/4;
		ERs_copy=matrix(0,nrow=m,ncol=4);
		for (i in 1:m)
		{
			start=ERs[i,1];
			end=ERs[i,2];
			aa=NPs[,1]-start;
			bb=NPs[,2]-end;
			aa=(aa>=0);
			bb=(bb<=0);
			aa=aa+bb;
			cc=which(aa==2);
			if (length(cc)>0)
			{
				ERs_copy[i,1]=NPs[min(cc),1];
				ERs_copy[i,2]=NPs[max(cc),2];
				ERs_copy[i,4]=length(cc);
			}
		}
		ERs=ERs_copy[ERs_copy[,1]>0,];

		a=length(ERs)/4;
		if (length(a)==0)
		{
			SERs_res <- list (SERs=c(0,0,0,0,0,0,0,0,0));
			return(SERs_res);
		}				
		if (a==1)
		{
			ERs=matrix(ERs);
			ERs=t(ERs);
		}		
		
		# count reads in SERs in case and control sample
		m=length(ERs)/4;	
		SERs_tmp=matrix(0, nrow=m, ncol=6);	
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			temp=ERs[i,];
			# count reads in case sample			
			# in 10bps resolution
			down=floor(temp[1]/10)+1;
			up=floor(temp[2]/10)+1;
			if (down<1) down=1;
			if (up<1) up=1;	
			case_d10=sum(RD_case_10[down:up])/length(down:up);			
			temp[3]=case_d10*(temp[2]-temp[1]+1);
			
			if (WithCtrl>0)
			{
				# count reads in control sample			
				# in 10bps resolution
				down=floor(temp[1]/10)+1;
				up=floor(temp[2]/10)+1;
				if (down<1) down=1;
				if (up<1) up=1;	
				d10=sum(RD_ctrl_10[down:up])/length(down:up);

				# in 1000bps resolution
				if ((temp[2]-temp[1]+1)<1000)
				{
					down=floor(temp[1]/1000)+1;
					up=floor(temp[2]/1000)+1;
					d1000=sum(RD_ctrl_1000[down:up])/length(down:up);
				} else {d1000=0;}

				# in 10000bps resolution
				if ((temp[2]-temp[1]+1)<10000)
				{
					down=floor(temp[1]/10000)+1;
					up=floor(temp[2]/10000)+1;
					d10000=sum(RD_ctrl_10000[down:up])/length(down:up);
				} else {d10000=0;}

				# estimate background reads 
				c=max(d10, d1000,d10000, Genome_density/Scale_factor);
				c=max(c*(temp[2]-temp[1]+1),c*2*Shift);
			} else {
				c=Genome_density/Scale_factor;
				c=max(c*(temp[2]-temp[1]+1),c*2*Shift);						
			}
			
			temp1=c();
			temp1=c(temp,c,0);
			SERs_tmp[i,]=temp1;
		}		
			
		# estimate p value of SERs
		m=dim(SERs_tmp)[1];
		SERs_tmp1=matrix(0, nrow=m, ncol=8);
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			a=SERs_tmp[i,3];
			b=SERs_tmp[i,5];
			p=-log10(ppois((a-1),Scale_factor*b,lower=FALSE));
			temp=c();
			ratio=(a/b)/Scale_factor;
			temp=cbind(t(SERs_tmp[i,]),p, ratio);
			SERs_tmp1[i,]=temp;
		}
		
		SERs_p=matrix(0, nrow=m, ncol=9);
		SERs_p[,1:8]=SERs_tmp1;
		SERs_p[,9]=0;
		SERs_res <- list (SERs=SERs_p);
		return(SERs_res);

}

getNPs <- function(SigWindow, ReadsLocation_case,ReadsLocation_control,Genome_density, Scale_factor,Shift)
# getNPs will merge significant overlapped final window to get narrow peaks (NPs) ...
# and calculate p value of the NPs

{
		Shift=as.numeric(Shift);
		if (length(ReadsLocation_case)==length(ReadsLocation_control))
		{
			if (sum(abs(ReadsLocation_case-ReadsLocation_control))==0) {WithCtrl=0;} else  WithCtrl=1;		
		} else WithCtrl=1;

		if (sum(SigWindow)==0)
		{
			NPs_res <- list (NPs=c(0,0,0,0,0,0,0,0,0));
			return(NPs_res);
		}
				
		a=length(SigWindow)/4;
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}
		
		SigWindow <- SigWindow[order(SigWindow[,1],SigWindow[,2]),];

		Tags=ReadsLocation_case;
		Input=ReadsLocation_control;

		# Build reads pipe up model in case samples for each bp in genome for summit calling
		Len=max(Tags);
		Len=Len+floor(4*Shift); 
		Pipeup_mat=matrix(0, nrow=Len, ncol=1);
		for (i in Tags)
		{
			a=i-Shift; # shift reads
			if (a<1) {a=1;}
			b=i+Shift; # shift reads
			Pipeup_mat[a:b]=Pipeup_mat[a:b]+1;
		}

		# Build reads density background in genome in control samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_ctrl_10=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10)+1;
			RD_ctrl_10[mn]=RD_ctrl_10[mn]+1;
		}
		RD_ctrl_10=RD_ctrl_10/10;
		
		# 1000 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/1000)+1;
		RD_ctrl_1000=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/1000)+1;
			RD_ctrl_1000[mn]=RD_ctrl_1000[mn]+1;
		}
		RD_ctrl_1000=RD_ctrl_1000/1000;

		# 10000bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10000)+1;
		RD_ctrl_10000=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10000)+1;
			RD_ctrl_10000[mn]=RD_ctrl_10000[mn]+1;
		}
		RD_ctrl_10000=RD_ctrl_10000/10000;

		# merge significant final windows to get NPs
		mm=sum(SigWindow[,2]-SigWindow[,1]+1);
		SigWindow_index=matrix(0, nrow=mm, ncol=1);
		mm=length(SigWindow[,1]);
		m=1;
		for (i in 1:mm)
		{
			x=SigWindow[i,2]-SigWindow[i,1]+1;
			SigWindow_index[m:(m+x-1)]=SigWindow[i,1]:SigWindow[i,2];
			m=m+x;
		}		
		SigWindow_index=unique(SigWindow_index);
		SigWindow_index=sort(SigWindow_index);
		SigWindow=matrix(0,nrow=length(SigWindow_index), ncol=4);
		m=SigWindow_index[1];
		for (i in 2:length(SigWindow_index))
		{
			if (SigWindow_index[i]-SigWindow_index[i-1]>1) {tmp=c(m,SigWindow_index[i-1],0,0); SigWindow[(i-1),]=tmp; m=SigWindow_index[i];}
		}
		tmp=c(m,SigWindow_index[i],0,0); SigWindow[i,]=tmp; m=SigWindow_index[i];
		SigWindow=SigWindow[which(SigWindow[,1]>0),];
		
		a=length(SigWindow)/4;
		if (length(a)==0)
		{
			NPs_res <- list (NPs=c(0,0,0,0,0,0,0,0,0));
			return(NPs_res);
		}				
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}
		SigWindow[,3]=SigWindow[,2]-SigWindow[,1]+2;
		NPs=SigWindow;
		NPs[,1]=Tags[NPs[,1]];
		NPs[,2]=Tags[NPs[,2]+1];
		
		m=length(NPs)/4;	
		NPs_tmp=matrix(0, nrow=m, ncol=6);	
		# count reads in NPs regions in control sample
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			temp=NPs[i,];
			if (WithCtrl>0)
			{
				# count reads in control sample			
				# in 10bps resolution
				down=floor(temp[1]/10)+1;
				up=floor(temp[2]/10)+1;
				if (down<1) down=1;
				if (up<1) up=1;	
				d10=sum(RD_ctrl_10[down:up])/length(down:up);

				# in 1000bps resolution
				if ((temp[2]-temp[1]+1)<1000)
				{
					down=floor(temp[1]/1000)+1;
					up=floor(temp[2]/1000)+1;
					d1000=sum(RD_ctrl_1000[down:up])/length(down:up);
				} else {d1000=0;}

				# in 10000bps resolution
				if ((temp[2]-temp[1]+1)<10000)
				{
					down=floor(temp[1]/10000)+1;
					up=floor(temp[2]/10000)+1;
					d10000=sum(RD_ctrl_10000[down:up])/length(down:up);
				} else {d10000=0;}

				# estimate background reads 
				c=max(d10, d1000,d10000, Genome_density/Scale_factor);
				c=max(c*(temp[2]-temp[1]+1),c*2*Shift);
			} else {
				c=Genome_density/Scale_factor;
				c=max(c*(temp[2]-temp[1]+1),c*2*Shift);						
			}
				
			# identify summit of NPs
			mat=Pipeup_mat[temp[1]:temp[2]];
			x=which(mat==max(mat));
			cc=max(mat);
			x=floor(median(x)+temp[1]-1);

			temp1=c();
			temp1=c(NPs[i,1:3], x, c,cc);
			NPs_tmp[i,]=temp1;
		}		
		
		# estimate p value of NPs
		m=dim(NPs_tmp)[1];
		NPs_tmp1=matrix(0, nrow=m, ncol=8);
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			a=NPs_tmp[i,3];
			b=NPs_tmp[i,5];
			p=-log10(ppois((a-1),Scale_factor*b,lower=FALSE));
			temp=c();
			ratio=(a/b)/Scale_factor;
			temp=cbind(t(NPs_tmp[i,]),p, ratio);
			NPs_tmp1[i,]=temp;
		}
		
		NPs_p=matrix(0, nrow=m, ncol=9);
		NPs_p[,1:8]=NPs_tmp1;
		NPs_p[,9]=0;
		NPs_res <- list (NPs=NPs_p);
		return(NPs_res);
		rm(Pipeup_mat);
}



getLERs <- function (SigWindow, ReadsLocation_case, ReadsLocation_control, Genome_density, Scale_factor, Shift, InputERs)
# getLERs will merge significant overlapped final window to get enriched regions (ERs) then it keeps the shortest ERs covering InputERs ...
# and calculate FR-RE and p value of the ERs.
# Last modified by Chao Wu on 2017-05-19

{

		Shift=as.numeric(Shift);
		if (length(ReadsLocation_case)==length(ReadsLocation_control))
		{
			if (sum(abs(ReadsLocation_case-ReadsLocation_control))==0) {WithCtrl=0;} else  WithCtrl=1;		
		} else WithCtrl=1;	
		
		if (sum(SigWindow)==0)
		{
			LERs_res <- list (LERs=c(0,0,0,0,0,0,0,0,0),FC=list(0));
			return(LERs_res);
		}		
		
		a=length(SigWindow)/4;
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}

		SigWindow <- SigWindow[order(SigWindow[,1],SigWindow[,2]),];

		a=length(SigWindow)/4;
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}
		
		Tags=ReadsLocation_case;
		Input=ReadsLocation_control;

		# Build reads density background in genome in case samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_case_10=matrix(0, nrow=Len, ncol=1);
		for (i in Tags)
		{
			mn=floor(i/10)+1;
			RD_case_10[mn]=RD_case_10[mn]+1;
		}
		RD_case_10=RD_case_10/10;

		# Build reads density background in genome in control samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_ctrl_10=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10)+1;
			RD_ctrl_10[mn]=RD_ctrl_10[mn]+1;
		}
		RD_ctrl_10=RD_ctrl_10/10;
		
		# 1000 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/1000)+1;
		RD_ctrl_1000=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/1000)+1;
			RD_ctrl_1000[mn]=RD_ctrl_1000[mn]+1;
		}
		RD_ctrl_1000=RD_ctrl_1000/1000;

		# 10000 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10000)+1;
		RD_ctrl_10000=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10000)+1;
			RD_ctrl_10000[mn]=RD_ctrl_10000[mn]+1;
		}
		RD_ctrl_10000=RD_ctrl_10000/10000;

		# merge significant final windows to get ERs
		mm=sum(SigWindow[,2]-SigWindow[,1]+1);
		SigWindow_index=matrix(0, nrow=mm, ncol=1);
		mm=length(SigWindow[,1]);
		m=1;
		for (i in 1:mm)
		{
			x=SigWindow[i,2]-SigWindow[i,1]+1;
			SigWindow_index[m:(m+x-1)]=SigWindow[i,1]:SigWindow[i,2];
			m=m+x;
		}		
		SigWindow_index=unique(SigWindow_index);
		SigWindow_index=sort(SigWindow_index);
		SigWindow=matrix(0,nrow=length(SigWindow_index), ncol=4);
		m=SigWindow_index[1];
		for (i in 2:length(SigWindow_index))
		{
			if (SigWindow_index[i]-SigWindow_index[i-1]>1) {tmp=c(m,SigWindow_index[i-1],0,0); SigWindow[(i-1),]=tmp; m=SigWindow_index[i];}
		}
		tmp=c(m,SigWindow_index[i],0,0); SigWindow[i,]=tmp; m=SigWindow_index[i];
		SigWindow=SigWindow[which(SigWindow[,1]>0),];

		a=length(SigWindow)/4;
		if (length(a)==0)
		{
			LERs_res <- list (LERs=c(0,0,0,0,0,0,0,0,0),FC=list(0));
			return(LERs_res);
		}				
		if (a==1)
		{
			SigWindow=matrix(SigWindow);
			SigWindow=t(SigWindow);
		}
		SigWindow[,3]=SigWindow[,2]-SigWindow[,1]+2;
		ERs=SigWindow;
		ERs[,1]=Tags[ERs[,1]];
		ERs[,2]=Tags[ERs[,2]+1];

		# calculate fold enrichment of neighboring InputERs
		ERs_filter=c();
		FC_count <- list();
		mn=0;
		m=length(ERs)/4;
		for (i in 1:m)
		{
			start=ERs[i,1];
			end=ERs[i,2];
			aa=InputERs[,1]-start;
			bb=InputERs[,2]-end;
			aa=(aa>=0);
			bb=(bb<=0);
			aa=aa+bb;
			cc=which(aa==2);
			if (length(cc)>0)
			{
				mn=mn+1;
				n=length(cc);
				ERs_filter=rbind(ERs_filter,c(InputERs[cc[1],1],InputERs[cc[n],2],0,length(cc)));
				if (n>1)
				{
					FoldEnrich_set=c(); 
					for (j in 1:(n-1))
					{
						temp=c(InputERs[cc[j],1],InputERs[cc[j+1],2]);
						down=floor(temp[1]/10)+1;
						up=floor(temp[2]/10)+1;
						if (down<1) down=1;
						if (up<1) up=1;
						d10_case=sum(RD_case_10[down:up])/length(down:up);
						d10_ctrl=sum(RD_ctrl_10[down:up])/length(down:up);

						Reads_case=d10_case*(temp[2]-temp[1]+1);
						Reads_ctrl=d10_ctrl*(temp[2]-temp[1]+1);
						if (WithCtrl==0) 
						{
							FC=d10_case/Genome_density;
						} else {
							FC=min(Reads_case/(Reads_ctrl*Scale_factor),d10_case/Genome_density);						
						}
						FoldEnrich_set=c(FoldEnrich_set,FC);
					}
					FC_count[[mn]] <- FoldEnrich_set;
				} else FC_count[[mn]] <- -1;
			}
		}

		a=length(ERs_filter)/4;
		if (length(a)==0)
		{
			LERs_res <- list (LERs=c(0,0,0,0,0,0,0,0,0),FC=list(0));
			return(LERs_res);
		}	
		if (a==1)
		{
			ERs_filter=matrix(ERs_filter);
			ERs_filter=t(ERs_filter);
		}	

		# count reads in LERs in case and control sample
		m=length(ERs_filter)/4;	
		LERs_tmp=matrix(0, nrow=m, ncol=6);	
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			temp=ERs_filter[i,];
			# count reads in case sample			
			# in 10bps resolution
			down=floor(temp[1]/10)+1;
			up=floor(temp[2]/10)+1;
			if (down<1) down=1;
			if (up<1) up=1;	
			case_d10=sum(RD_case_10[down:up])/length(down:up);			
			temp[3]=case_d10*(temp[2]-temp[1]+1);
			
			if (WithCtrl>0)
			{
				# count reads in control sample			
				# in 10bps resolution
				down=floor(temp[1]/10)+1;
				up=floor(temp[2]/10)+1;
				if (down<1) down=1;
				if (up<1) up=1;	
				d10=sum(RD_ctrl_10[down:up])/length(down:up);

				# in 1000bps resolution
				if ((temp[2]-temp[1]+1)<1000)
				{
					down=floor(temp[1]/1000)+1;
					up=floor(temp[2]/1000)+1;
					d1000=sum(RD_ctrl_1000[down:up])/length(down:up);
				} else {d1000=0;}

				# in 10000bps resolution
				if ((temp[2]-temp[1]+1)<10000)
				{
					down=floor(temp[1]/10000)+1;
					up=floor(temp[2]/10000)+1;
					d10000=sum(RD_ctrl_10000[down:up])/length(down:up);
				} else {d10000=0;}

				# estimate background reads 
				c=max(d10, d1000,d10000, Genome_density/Scale_factor);
				c=max(c*(temp[2]-temp[1]+1),c*2*Shift);
			} else {
				c=Genome_density/Scale_factor;
				c=max(c*(temp[2]-temp[1]+1),c*2*Shift);						
			}
			
			temp1=c();
			temp1=c(temp,c,0);
			LERs_tmp[i,]=temp1;
		}		
			
		# estimate p value of LERs
		m=dim(LERs_tmp)[1];
		LERs_tmp1=matrix(0, nrow=m, ncol=8);
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			a=LERs_tmp[i,3];
			b=LERs_tmp[i,5];
			p=-log10(ppois((a-1),Scale_factor*b,lower=FALSE));
			temp=c();
			ratio=(a/b)/Scale_factor;
			temp=cbind(t(LERs_tmp[i,]),p, ratio);
			LERs_tmp1[i,]=temp;
		}
		
		LERs_p=matrix(0, nrow=m, ncol=9);
		LERs_p[,1:8]=LERs_tmp1;
		LERs_p[,9]=0;
		LERs_res <- list (LERs=LERs_p,FC=FC_count);
		return(LERs_res);
}


callSigFinalWindow <- function (FinalWindow_p, SigThreshold)
# Correct the p value of all the final windows by bonferroni corrected test ... 
# and call the significant final windows.
# Last modified by Chao Wu on 2017-05-16.

{
	# Correct multiple test p value by bonferroni test
	cor_p=p.adjust(FinalWindow_p$FinalWindow_p_tmp,"bonferroni");
	
	# Identify the p threshold for significant final windows
	index_p=FinalWindow_p$FinalWindow_p_tmp[cor_p<SigThreshold];
	cor_p_threshold=max(index_p);

	a=which(FinalWindow_p$FinalWindow_p[,2]<=cor_p_threshold);
	if (length(a)==0) 
	{
		print ("no significant final windows is detected");
		SigWindow=t(c(0,0,0,0)); SigWindow=as.matrix(SigWindow);
		SigWindow_res <- list (SigWindow=SigWindow);
		return (SigWindow_res);
	} else {
		SigWindow=FinalWindow_p$FinalWindow_p[a,];
		m=SigWindow[,3];
		start=SigWindow[,1];
		end=SigWindow[,1]+m-1;
		SigWindow=cbind(start, end, 0, 0);
		SigWindow_res <- list (SigWindow=SigWindow);	
		return (SigWindow_res);	
	}
}




callSERs_Main <- function (ReadsLocation_case, ReadsLocation_control, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit,NPs)
# This function is the main function to call SERs.
# ReadsLocation_case and ReadsLocation_control are the sorted mapping locations of reads from case and control samples.
# In the ChIP-Seq data without control, please input ReadsLocation_case and ReadsLocation_control both with the sorted mapping location of reads from case sample.
# Species identifies the input speicies.
# InitialWindowSize sets the number of bins in an initial window.
# Genome_density is the average mapped reads density in the genome in case sample.
# Scale_factor is the ratio of the mapped reads depth between case and control samles.
# SigThreshold is the parameter of bonforroni corrected p value to call Enrichment signals.
# Shift is the estimated shift parameter
# StepWindowSize_limit is the length parameter to set the max length of an step window.
# BinSize_limit is the size parameter to set the max size of bins in an step window.
# NPs is the detected narrow peaks
# Last modified by Chao Wu on 2017-05-19.

{

	finalWins_p=c();
	Steps=3;

	# Calculate p value of each final widnow
	finalWins_p<- scoreFinalWindow (ReadsLocation_case, ReadsLocation_control, InitialWindowSize, Steps, Genome_density, Scale_factor, StepWindowSize_limit ,BinSize_limit)

	# Identify significant final windows
	sig_finalWins<- callSigFinalWindow (finalWins_p, SigThreshold)

	
	# Identify SERs
	SERs <-getSERs (sig_finalWins$SigWindow, ReadsLocation_case, ReadsLocation_control, Genome_density, Scale_factor,Shift, NPs)
	return(SERs);
}


callNPs_Main <- function (ReadsLocation_case, ReadsLocation_control, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit)
# This function is the main function to call NPs.
# ReadsLocation_case and ReadsLocation_control are the sorted mapping locations of reads from case and control samples.
# In the ChIP-Seq data without control, please input ReadsLocation_case and ReadsLocation_control both with the sorted mapping location of reads from case sample.
# Species identifies the input speicies.
# InitialWindowSize sets the number of bins in an initial window.
# Genome_density is the average mapped reads density in the genome in case sample.
# Scale_factor is the ratio of the mapped reads depth between case and control samles.
# SigThreshold is the parameter of bonforroni corrected p value to call Enrichment signals.
# Shift is the estimated shift parameter
# StepWindowSize_limit is the length parameter to set the max length of an step window.
# BinSize_limit is the size parameter to set the max size of bins in an step window.
# Last modified by Chao Wu on 2017-05-17.

{

	finalWins_p=c();
	Steps=3;

	# Calculate p value of each final widnow
	finalWins_p<- scoreFinalWindow (ReadsLocation_case, ReadsLocation_control, InitialWindowSize, Steps, Genome_density, Scale_factor, StepWindowSize_limit ,BinSize_limit)

	# Identify significant final windows
	sig_finalWins<- callSigFinalWindow (finalWins_p, SigThreshold)

	# Identify NPs
	NPs<-getNPs (sig_finalWins$SigWindow, ReadsLocation_case,ReadsLocation_control,Genome_density, Scale_factor,Shift)
	
	return(NPs);
}


callLERs_Main <- function (ReadsLocation_case, ReadsLocation_control, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit,InputERs)
# This function is the main function to call LERs.
# ReadsLocation_case and ReadsLocation_control are the sorted mapping locations of reads from case and control samples.
# In the ChIP-Seq data without control, please input ReadsLocation_case and ReadsLocation_control both with the sorted mapping location of reads from case sample.
# Species identifies the input speicies.
# InitialWindowSize sets the number of bins in an initial window.
# Genome_density is the average mapped reads density in the genome in case sample.
# Scale_factor is the ratio of the mapped reads depth between case and control samles.
# SigThreshold is the parameter of bonforroni corrected p value to call Enrichment signals.
# Shift is the estimated shift parameter
# StepWindowSize_limit is the length parameter to set the max length of an step window.
# BinSize_limit is the size parameter to set the max size of bins in an step window.
# InputERs is the input ERs.
# Last modified by Chao Wu on 2017-05-19.


{
	finalWins_p=c();
	Steps=3;

	# Calculate p value of each final widnow
	finalWins_p<- scoreFinalWindow (ReadsLocation_case, ReadsLocation_control, InitialWindowSize, Steps, Genome_density, Scale_factor, StepWindowSize_limit ,BinSize_limit)

	# Identify significant final windows
	sig_finalWins<- callSigFinalWindow (finalWins_p, SigThreshold)
	
	# Identify LERs
	LERs <-getLERs (sig_finalWins$SigWindow, ReadsLocation_case, ReadsLocation_control, Genome_density, Scale_factor,Shift, InputERs)
	return(LERs);
}

addGenomeInfo <- function (species, ratio)
# The genome annotation information
# Please add the genome in the function if the genome information of your ChIP-Seq is not listed.  
# Last modificated by Chao Wu on 2017-05-17

{
		x=0; # genome size
		y=c(); # chr information
		if (species == "hg18")
		{
			x=3022646526*ratio;		
			y=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX");
		}
		if (species == "hg19")
		{
			x=3036303846*ratio;
			y=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX");
		}
		if (species == "mm8")
		{
			x=2628048285*ratio;
			y=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chrX");
		}
		if (species == "mm9")
		{
			x=2638992663*ratio;
			y=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chrX");			
		}
		
		genome_info <- list(size=x, chr_info=y);
		return(genome_info);
}

calculateFR_RE <- function (ReadsLocation_case, ReadsLocation_control, Genome_density, Scale_factor, Input_chr)
# calculateFR_RE will calculate the fold enrichment between neighboring Input ERs (including the neighbor ERs themseleves).
# Last modified by Chao Wu on 2017-05-19

{
		if (length(ReadsLocation_case)==length(ReadsLocation_control))
		{
			if (sum(abs(ReadsLocation_case-ReadsLocation_control))==0) {WithCtrl=0;} else  WithCtrl=1;		
		} else WithCtrl=1;	
		
		Tags=ReadsLocation_case;
		Input=ReadsLocation_control;

		# Build reads density background in genome in case samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_case_10=matrix(0, nrow=Len, ncol=1);
		for (i in Tags)
		{
			mn=floor(i/10)+1;
			RD_case_10[mn]=RD_case_10[mn]+1;
		}
		RD_case_10=RD_case_10/10;

		# Build reads density background in genome in control samples 
		# 10 bps resolution reads density background
		Len=max(c(Tags,Input));
		Len=floor(Len/10)+1;
		RD_ctrl_10=matrix(0, nrow=Len, ncol=1);
		for (i in Input)
		{
			mn=floor(i/10)+1;
			RD_ctrl_10[mn]=RD_ctrl_10[mn]+1;
		}
		RD_ctrl_10=RD_ctrl_10/10;

		m=length(Input_chr)/2;
		Neighbor_start=Input_chr[1:(m-1),1];
		Neighbor_end=Input_chr[2:m,2];
		Nieghbor=cbind(Neighbor_start,Neighbor_end);

		m=length(Nieghbor)/2;
		FC=c();
		for (i in 1:m)
		{
			if (i%%1000==0) {print (i);}
			temp=Nieghbor[i,];
			# count reads in case sample			
			# in 10bps resolution
			down=floor(temp[1]/10)+1;
			up=floor(temp[2]/10)+1;
			if (down<1) down=1;
			if (up<1) up=1;	
			case_d10=sum(RD_case_10[down:up])/length(down:up);
		
			if (WithCtrl>0)
			{
				# count reads in control sample			
				# in 10bps resolution
				down=floor(temp[1]/10)+1;
				up=floor(temp[2]/10)+1;
				if (down<1) down=1;
				if (up<1) up=1;	
				d10=sum(RD_ctrl_10[down:up])/length(down:up);

				# estimate background reads 
				c=max(d10*Scale_factor, Genome_density);
			} else {
				c=Genome_density;				
			}		
			FC=c(FC,(case_d10/c))
		}
		FC_res <- list (FC=FC);
}


multicore_BinSizeLimit <- function(input_para)
# This function will calculate the length distribution of NPs under different BinSize_limit candidates
# Last modified by Chao Wu on 2017-05-18

{
	source("../CLUES_functions.R"); ## need functions in the script to run multi-copy examples
	para_set=unlist(strsplit(input_para, ","));
	
	Chr_test=unlist(strsplit(para_set[1], "M"));
	case_name=para_set[2];
	control_name=para_set[3];
	Shift=para_set[4];
	Species=para_set[5];
	Genome_density=as.numeric(para_set[6]);
	Scale_factor=as.numeric(para_set[7]);
	SigThreshold=as.numeric(para_set[8]);
	xm=as.numeric(para_set[9]);
	
	BinSize_limit=floor(2*as.numeric(Shift)*xm);
	#print ("BinSize_limit");
	#print (BinSize_limit);

	# Estimate InitialWindowSize
	InitialWindowSize=max(floor(2*as.numeric(Shift)*Genome_density*5),20);
	InitialWindowSize=min(InitialWindowSize, 100);
	NPs_pvalue=c();
	NPs_set=c();
	for (i in 1:length(Chr_test))
	{
		print (Chr_test[i]);
		case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",Chr_test[i],xm,".bed",sep="");
		case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", Chr_test[i],xm,".bed > ",Chr_test[i],xm,"_seq_case.txt",sep="");
		system(case_tmp);
		system(case_tmp1);
		rm_file <- paste ("rm ",Chr_test[i],xm,".bed",sep="");
		system(rm_file);
		case <- paste (Chr_test[i],xm,"_seq_case.txt",sep="");
		rm_case <- paste ("rm ",Chr_test[i],xm,"_seq_case.txt",sep="");
		control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",Chr_test[i],xm,".bed",sep="");
		control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", Chr_test[i],xm,".bed > ",Chr_test[i],xm,"_control.txt",sep=""); 
		system(control_tmp);
		system(control_tmp1);
		system(rm_file);
		rm_control <- paste ("rm ",Chr_test[i],xm,"_control.txt",sep="");
		control<-paste (Chr_test[i],xm,"_control.txt",sep="");
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
			StepWindowSize_limit=10000000000;
			
			
			NPs_res <-callNPs_Main (seq_case, seq_input, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit)
			if (sum(NPs_res$NPs)!=0)
			{
				tmp=matrix(0,nrow=length(NPs_res$NPs)/9,ncol=10);
				tmp[,2:10]=NPs_res$NPs;
				NPs_pvalue=c(NPs_pvalue, NPs_res$NPs[,7]);
				NPs_set=rbind(NPs_set,tmp);
			}
		}
		system(rm_case);
		system(rm_control);
	}

	qfdr=p.adjust(10^(-NPs_pvalue),"BY"); # multiple pvalue correction
	qfdr=-log10(qfdr);
	NPs_set[,10]=qfdr;
	NPs_set=NPs_set[which(as.numeric(NPs_set[,10])>-log10(SigThreshold)),];
				
	# get the length distribution of NPs
	if (length(NPs_set)<1000) 			
	{
		NPs_lenDistr=c(BinSize_limit,rep(0,10));					
	} else {
		NPs_len=NPs_set[,3]-NPs_set[,2]+1;
		NPs_lenDistr=c(BinSize_limit,quantile(NPs_len,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)));
	}
	
	NPs_Distr <- list (NPs_lenDistr=NPs_lenDistr);
	return(NPs_Distr);
}	

multicore_LERs <- function (input_para)
# This function will call LERs.
# Last modified by Chao Wu on 2017-05-19

{
	
	source("../CLUES_functions.R"); ## need functions in the script to run multi-copy examples

	para_set=unlist(strsplit(input_para, ","));
	chr_test=para_set[1];
	case_name=para_set[2];
	control_name=para_set[3];
	Shift=para_set[4];
	Species=para_set[5];
	Genome_density=as.numeric(para_set[6]);
	Scale_factor=as.numeric(para_set[7]);
	SigThreshold=as.numeric(para_set[8]);
	LERs_StepWindowSize_limit=as.numeric(para_set[9]);
	Input=as.character(para_set[10]);	

	InputERs=read.table(Input, header=FALSE);
	# Estimate InitialWindowSize
	InitialWindowSize=max(floor(LERs_StepWindowSize_limit*Genome_density*5),20);
	InitialWindowSize=min(500,InitialWindowSize);

	LERs_set=c();
	for (i in 1:length(chr_test))
	{
		print (chr_test[i]);
		case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",chr_test[i],".bed",sep="");
		case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", chr_test[i],".bed > ",chr_test[i],"_seq_case.txt",sep="");
		system(case_tmp);
		system(case_tmp1);
		rm_file <- paste ("rm ",chr_test[i],".bed",sep="");
		system(rm_file);
		case <- paste (chr_test[i],"_seq_case.txt",sep="");
		rm_case <- paste ("rm ",chr_test[i],"_seq_case.txt",sep="");
	
		control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",chr_test[i],".bed",sep="");
		control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", chr_test[i],".bed > ",chr_test[i],"_control.txt",sep=""); 
		system(control_tmp);
		system(control_tmp1);
		system(rm_file);
		rm_control <- paste ("rm ",chr_test[i],"_control.txt",sep="");
		control<-paste (chr_test[i],"_control.txt",sep="");
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

			x=which(InputERs[,1]==chr_test[i]);
			if (length(x)<3) next;		
			InputERs_chr=cbind(as.numeric(InputERs[x,2]),as.numeric(InputERs[x,3]),as.numeric(InputERs[x,4]));			
			BinSize_limit=10000000000000;
      LERs_res <- callLERs_Main (seq_case, seq_input, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, LERs_StepWindowSize_limit ,BinSize_limit,InputERs_chr)

			if (sum(LERs_res$LERs)!=0)
			{
				tmp=matrix(0,nrow=length(LERs_res$LERs)/9,ncol=10);
				tmp[,2:10]=LERs_res$LERs;
				tmp[,1]=chr_test[i];
				LERs_set=rbind(LERs_set,tmp);				
			}
		}	
		system(rm_case);
		system(rm_control);
	}
	LERs_set<-list(LERs=LERs_set);
	return(LERs_set);
}

multicore_NPs <- function(input_para)
# This function will call NPs.
# Last modified by Chao Wu on 2017-05-18

{
		source("../CLUES_functions.R"); ## need functions in the script to run multi-copy examples

		para_set=unlist(strsplit(input_para, ","));
		Chr=para_set[1];
		case_name=para_set[2];
		control_name=para_set[3];
		Shift=para_set[4];
		Species=para_set[5];
		Genome_density=as.numeric(para_set[6]);
		Scale_factor=as.numeric(para_set[7]);
		SigThreshold=as.numeric(para_set[8]);
		BinSize_limit=as.numeric(para_set[9]);
		
		# Estimate InitialWindowSize
		InitialWindowSize=max(floor(2*as.numeric(Shift)*Genome_density*5),20);
		InitialWindowSize=min(InitialWindowSize, 100);

 		NPs_set=c();
		case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr,"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",Chr,".bed",sep="");
		case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", Chr,".bed > ",Chr,"_seq_case.txt",sep="");
		system(case_tmp);
		system(case_tmp1);
		rm_file <- paste ("rm ",Chr,".bed",sep="");
		system(rm_file);
		case <- paste (Chr,"_seq_case.txt",sep="");
		rm_case <- paste ("rm ",Chr,"_seq_case.txt",sep="");

		control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",Chr,"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",Chr,".bed",sep="");
		control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", Chr,".bed > ",Chr,"_control.txt",sep=""); 
		system(control_tmp);
		system(control_tmp1);
		system(rm_file);
		rm_control <- paste ("rm ",Chr,"_control.txt",sep="");
		control<-paste (Chr,"_control.txt",sep="");
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
			
			StepWindowSize_limit=10000000000;			
			NPs_res <-callNPs_Main (seq_case, seq_input, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit)
			if (sum(NPs_res$NPs)!=0)
			{
				tmp=matrix(0,nrow=length(NPs_res$NPs)/9,ncol=10);
				tmp[,2:10]=NPs_res$NPs;
				tmp[,1]=Chr;
				NPs_set=rbind(NPs_set,tmp);
			}
		}
		system(rm_case);
		system(rm_control);
		NPs_res<-list(NPs=NPs_set);
		return(NPs_res);
}



multicore_SERs <- function(input_para)
# This function will call SERs.
# Last modified by Chao Wu on 2017-05-19

{

	source("../CLUES_functions.R"); ## need functions in the script to run multi-copy examples

	para_set=unlist(strsplit(input_para, ","));
	chr=para_set[1];
	case_name=para_set[2];
	control_name=para_set[3];
	Shift=para_set[4];
	Species=para_set[5];
		

		
	Genome_density=as.numeric(para_set[6]);
	Scale_factor=as.numeric(para_set[7]);
	SigThreshold=as.numeric(para_set[8]);
	StepWindowSize_limit=as.numeric(para_set[9]);
	NPs_file=as.character(para_set[10]);
	# Estimate InitialWindowSize
	InitialWindowSize=max(floor(StepWindowSize_limit*Genome_density*5),20);
	InitialWindowSize=min(500,InitialWindowSize);		
	NPs=read.table(NPs_file, header=FALSE);
	SERs_set=c();
	for (i in 1:length(chr))
	{
		print (chr[i]);
		case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",chr[i],".bed",sep="");
		case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", chr[i],".bed > ",chr[i],"_seq_case.txt",sep="");
		system(case_tmp);
		system(case_tmp1);
		rm_file <- paste ("rm ",chr[i],".bed",sep="");
		system(rm_file);
		case <- paste (chr[i],"_seq_case.txt",sep="");
		rm_case <- paste ("rm ",chr[i],"_seq_case.txt",sep="");
	
		control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr[i],"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",chr[i],".bed",sep="");
		control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", chr[i],".bed > ",chr[i],"_control.txt",sep=""); 
		system(control_tmp);
		system(control_tmp1);
		system(rm_file);
		rm_control <- paste ("rm ",chr[i],"_control.txt",sep="");
		control<-paste (chr[i],"_control.txt",sep="");
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
			x=which(NPs[,1]==chr[i]);
			if (length(x)<3) next;
			NPs_chr=cbind(as.numeric(NPs[x,2]),as.numeric(NPs[x,3]),as.numeric(NPs[x,4]));
			BinSize_limit=10000000000000;		
				
			SERs_res <-callSERs_Main (seq_case, seq_input, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit,NPs_chr);
			if (sum(SERs_res$SERs)!=0)
			{
				tmp=matrix(0,nrow=length(SERs_res$SERs)/9,ncol=10);
				tmp[,2:10]=SERs_res$SERs;
				tmp[,1]=chr[i];
				SERs_set=rbind(SERs_set,tmp);				
			}
		}
		system(rm_case);
		system(rm_control);
	}	
	SERs_set<-list(SERs=SERs_set);
	return(SERs_set);
}





multicore_StepWindowSize_LERs <- function (input_para)
# This function will calculate the FR-RE of input ERs under different StepWindowSize_limit candidates
# Last modified by Chao Wu on 2017-05-19


{
	
	source("../CLUES_functions.R"); ## need functions in the script to run multi-copy examples

	para_set=unlist(strsplit(input_para, ","));
	chr_test=unlist(strsplit(para_set[1], "M"))
	case_name=para_set[2];
	control_name=para_set[3];
	Shift=para_set[4];
	Species=para_set[5];	
	Genome_density=as.numeric(para_set[6]);
	Scale_factor=as.numeric(para_set[7]);
	SigThreshold=as.numeric(para_set[8]);
	StepWindowSize_limit=as.numeric(para_set[9]);
	Input=as.character(para_set[10]);	
	xm=as.numeric(para_set[11]);
		
	InputERs=read.table(Input, header=FALSE);
	if (StepWindowSize_limit>99999) return (c(0,0,0,0,0,0));
	# Estimate InitialWindowSize
	InitialWindowSize=max(floor(StepWindowSize_limit*Genome_density*5),20);
	InitialWindowSize=min(500,InitialWindowSize);
	FC_set=c();
	LERs_pvalue=c();
	LERs_set=c();
	Input_Count=0;
	for (i in 1:length(chr_test))
	{
		print (chr_test[i]);
		case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",chr_test[i],xm,".bed",sep="");
		case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", chr_test[i],xm,".bed > ",chr_test[i],xm,"_seq_case.txt",sep="");
		system(case_tmp);
		system(case_tmp1);
		rm_file <- paste ("rm ",chr_test[i],xm,".bed",sep="");
		system(rm_file);
		case <- paste (chr_test[i],xm,"_seq_case.txt",sep="");
		rm_case <- paste ("rm ",chr_test[i],xm,"_seq_case.txt",sep="");
	
		control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",chr_test[i],xm,".bed",sep="");
		control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", chr_test[i],xm,".bed > ",chr_test[i],xm,"_control.txt",sep=""); 
		system(control_tmp);
		system(control_tmp1);
		system(rm_file);
		rm_control <- paste ("rm ",chr_test[i],xm,"_control.txt",sep="");
		control<-paste (chr_test[i],xm,"_control.txt",sep="");
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

			x=which(InputERs[,1]==chr_test[i]);
			if (length(x)<3) next;
			Input_Count=Input_Count+length(x);			
			InputERs_chr=cbind(as.numeric(InputERs[x,2]),as.numeric(InputERs[x,3]),as.numeric(InputERs[x,4]));			
			BinSize_limit=10000000000000;
      LERs_res <- callLERs_Main (seq_case, seq_input, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit,InputERs_chr)

			if (sum(LERs_res$LERs)!=0)
			{
				mn=length(LERs_res$FC);
				FC_res <- list();
				for (mmm in 1:mn)
				{
					FC_res[[mmm]]=LERs_res$FC[[mmm]];
				}
				FC_set=c(FC_set, FC_res);
				tmp=matrix(0,nrow=length(LERs_res$LERs)/9,ncol=10);
				tmp[,2:10]=LERs_res$LERs;
				LERs_pvalue=c(LERs_pvalue, LERs_res$LERs[,7]);
				LERs_set=rbind(LERs_set,tmp);				
			}
		}	
		system(rm_case);
		system(rm_control);
	}
	LERs_Inputcount=sum(as.numeric(LERs_set[,5]));
	qfdr=p.adjust(10^(-LERs_pvalue),"BY");
	qfdr=-log10(qfdr);
	LERs_set[,10]=qfdr;
	#FC_set=FC_set[as.numeric(LERs_set[,10])>-log10(SigThreshold)];
	#LERs_set=LERs_set[as.numeric(LERs_set[,10])>-log10(SigThreshold),];
	

	#if ((LERs_Inputcount/Input_Count)<0.95|length(LERs_set)<200)
	if (length(LERs_set)<200)
	{
		return (c(0,0,0,0,0,0));				
	}
	FC_count=c();
	for (mmm in 1:length(FC_set))
	{
		tmp=FC_set[[mmm]];
		for (nnn in 1:length(tmp))
		{
			FC_count=c(FC_count,tmp[[nnn]])
		}
	}
	FC_count=FC_count[FC_count>0];
	if (length(FC_count)==0) 
	{
		return (c(0,0,0,0,0,0));
	}	
	FC_count=sort(FC_count);
	StepWindowSize_LERs=c(StepWindowSize_limit,length(LERs_set[,10]),length(FC_count),FC_count[max(length(FC_count)*0.01,1)],FC_count[max(length(FC_count)*0.02,1)],FC_count[max(length(FC_count)*0.05,1)]);
	return(StepWindowSize_LERs);
}

multicore_StepWindowSize_SERs <- function(input_para)
# This function will calculate the FR-D of SERs under different StepWindowSize_limit candidates
# Last modified by Chao Wu on 2017-05-19

{

	source("../CLUES_functions.R"); ## need functions in the script to run multi-copy examples

	para_set=unlist(strsplit(input_para, ","));
	chr_test=unlist(strsplit(para_set[1], "M"));
	case_name=para_set[2];
	control_name=para_set[3];
	Shift=para_set[4];
	Species=para_set[5];
	Genome_density=as.numeric(para_set[6]);
	Scale_factor=as.numeric(para_set[7]);
	SigThreshold=as.numeric(para_set[8]);
	StepWindowSize_limit=as.numeric(para_set[9]);
	NPs_file=as.character(para_set[10]);
	xm=as.numeric(para_set[11]);
	
	NPs=read.table(NPs_file, header=FALSE);
	# Estimate InitialWindowSize
	InitialWindowSize=max(floor(StepWindowSize_limit*Genome_density*5),20);
	InitialWindowSize=min(500,InitialWindowSize);
	SERs_pvalue=c();
	SERs_set=c();
	NPschr_Count=0;
	for (i in 1:length(chr_test))
	{
		print (chr_test[i]);
		case_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", case_name, " > ",chr_test[i],xm,".bed",sep="");
		case_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+",Shift,") }' ", chr_test[i],xm,".bed > ",chr_test[i],xm,"_seq_case.txt",sep="");
		system(case_tmp);
		system(case_tmp1);
		rm_file <- paste ("rm ",chr_test[i],xm,".bed",sep="");
		system(rm_file);
		case <- paste (chr_test[i],xm,"_seq_case.txt",sep="");
		rm_case <- paste ("rm ",chr_test[i],xm,"_seq_case.txt",sep="");
	
		control_tmp <- paste("awk '{OFS==\"\\t\"} {if($1==\"",chr_test[i],"\") print $1, $2, $3, $4, $5, $6}' ", control_name, " > ",chr_test[i],xm,".bed",sep="");
		control_tmp1 <- paste("awk '{OFS==\"\\t\"} {if($6==\"-\") print ($3-", Shift, "); else print ($2+", Shift, ") }' ", chr_test[i],xm,".bed > ",chr_test[i],xm,"_control.txt",sep=""); 
		system(control_tmp);
		system(control_tmp1);
		system(rm_file);
		rm_control <- paste ("rm ",chr_test[i],xm,"_control.txt",sep="");
		control<-paste (chr_test[i],xm,"_control.txt",sep="");
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

			x=which(NPs[,1]==chr_test[i]);
			if (length(x)<3) next;
			NPschr_Count=NPschr_Count+length(x);
			NPs_chr=cbind(as.numeric(NPs[x,2]),as.numeric(NPs[x,3]),as.numeric(NPs[x,4]));
			BinSize_limit=10000000000000;
			SERs_res <-callSERs_Main (seq_case, seq_input, Species, InitialWindowSize, Genome_density, Scale_factor,SigThreshold, Shift, StepWindowSize_limit ,BinSize_limit,NPs_chr);
			if (sum(SERs_res$SERs)!=0)
			{
				tmp=matrix(0,nrow=length(SERs_res$SERs)/9,ncol=10);
				tmp[,2:10]=SERs_res$SERs;
				SERs_pvalue=c(SERs_pvalue, SERs_res$SERs[,7]);
				SERs_set=rbind(SERs_set,tmp);				
			}
		}
		system(rm_case);
		system(rm_control);
	}
	SERs_NPscount=sum(as.numeric(SERs_set[,5]));
	qfdr=p.adjust(10^(-SERs_pvalue),"BY");
	qfdr=-log10(qfdr);
	SERs_set[,10]=qfdr;
	#SERs_set=SERs_set[as.numeric(SERs_set[,10])>-log10(SigThreshold),];


	#if ((SERs_NPscount/NPschr_Count)<0.95|length(SERs_set)<200)
	if (length(SERs_set)<200)
	{
		StepWindowSize_res<-list(StepWindowSize_stat=c(0,0,1));
		return(StepWindowSize_res);					
	}

	# Estimate fragment rate by distance(FR-D) of the SERs
	SERs_len=as.numeric(SERs_set[,3])-as.numeric(SERs_set[,2])+1;
	LeftSERs=SERs_len[1:(length(SERs_len)-1)];
	RightSERs=SERs_len[2:length(SERs_len)];
	max_Len=cbind(LeftSERs,RightSERs);
	max_Len=1/2*apply(max_Len,1,max);
	
	Between_SERs=as.numeric(SERs_set[2:length(SERs_set[,2]),2])-as.numeric(SERs_set[1:(length(SERs_set[,2])-1),3])+1;
	Between_SERs[Between_SERs<0]=max(max_Len);
	tmp=Between_SERs<max_Len;
	FR_SERs=which(tmp==1);
	FRD_SERs=length(FR_SERs)/(length(max_Len)+1);

	tmp=c(StepWindowSize_limit, length(SERs_set[,2]), FRD_SERs);
	StepWindowSize_res<-list(StepWindowSize_stat=tmp);
	return(StepWindowSize_res);
}
