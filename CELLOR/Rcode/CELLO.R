# CELLOR
#frame_files <- lapply(sys.frames(), function(x) x$ofile)
#frame_files <- Filter(Negate(is.null), frame_files)
#TOPDIR <- dirname(frame_files[[length(frame_files)]])
#setwd(TOPDIR)

library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
library("igraph")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}


#function 1: read in savi report & preprocessing the data
mutRead <- function(savi_report_file,refdepth_Blood_cutoff,altdepth_Blood_cutoff,freq_cutoff){
	savi <- read.delim(savi_report_file)
	savi <- savi[which( savi$refdepth_Blood >= refdepth_Blood_cutoff & savi$altdepth_Blood <= altdepth_Blood_cutoff & savi$Sgt1_max_frequency >= freq_cutoff ),]
	savi[is.na(savi)] <- 0
	#savi <- savi[which(!(savi$Effect_Impact %in% 'LOW')),]

 	results <- savi
}


#function 2: analysis the number of mutations in primary & recurrent samples; and analysis the mutation information of selected genes
mutStats <- function(savi_table, selected_geneList, freq_cutoff,remove_LOW=TRUE){

	savi <- savi_table
	cutoff <- freq_cutoff
	if(remove_LOW == TRUE){
		savi<-savi[which(!(savi$Effect_Impact %in% 'LOW')),]  ## exclude the synonymous_variant
	}

	case <- unique(savi$CaseID)


	#counting mutations
	mut.primary <- rep(0,length(case))
	mut.recurrence <- rep(0,length(case))
	mut.common <- rep(0,length(case))

	for(i in 1:length(case)){
		mut.common[i] <- length(which(savi$CaseID == case[i] & savi$Primary_freq >= cutoff & savi$Recurrent_freq >= cutoff))
		mut.primary[i] <- length(which(savi$CaseID == case[i] & savi$Primary_freq >= cutoff & savi$Recurrent_freq < cutoff))
		mut.recurrence[i] <- length(which(savi$CaseID == case[i] & savi$Primary_freq < cutoff & savi$Recurrent_freq >= cutoff))
	}

	mutation_num_table<- data.frame( case, mut.primary, mut.common, mut.recurrence )
	colnames(mutation_num_table)<-c('Patients','Primary','Common','Recurrent')


	#analysis for selected genes
	selGene <- selected_geneList

	gene.Matrix <- rep('N',length(case))
	for( i in 2:length(selGene)){
		gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
	}

	savi.Sel <- savi[which(savi$Gene_Name %in% selGene),]

	for(i in 1:length(case)){
		for(j in 1:length(selGene)){
			temp.P <- savi.Sel$Primary_freq[which(savi.Sel$CaseID == case[i] & savi.Sel$Gene_Name == selGene[j])]
			temp.R <- savi.Sel$Recurrent_freq[which(savi.Sel$CaseID == case[i] & savi.Sel$Gene_Name == selGene[j])]
			if( any(temp.P >= cutoff) & any(temp.R >= cutoff) ){
				gene.Matrix[i,j] <- 'C'
			}
			else if( any(temp.P	>= cutoff) ){
				gene.Matrix[i,j] <- 'P'
			}
			else if( any(temp.R	>= cutoff) ){
				gene.Matrix[i,j] <- 'R'
			}
			else{}
		}
	}

	colnames(gene.Matrix) <- selGene
	rownames(gene.Matrix) <- case

	returnList <- list("mutNum.table" = mutation_num_table, "mutGenes.table" = gene.Matrix)
	return(returnList)

}


#function 3: plot the mutation landscape
mutLandscape <- function(mutation_num_table, mutation_gene_table){

	plot_1a.data <- data.frame(rep(mutation_num_table$Patients,3), c( c(mutation_num_table$Recurrent), c(mutation_num_table$Common), c(mutation_num_table$Primary) ),
	                            c( rep('a_Recurrent',length(mutation_num_table[,1])), rep('b_Common',length(mutation_num_table[,1])), rep('c_Primary',length(mutation_num_table[,1])) ) )


	colnames(plot_1a.data)<-c('Patients','Mutations','Groups')
	plot_1a.data$Patients<-as.character(plot_1a.data$Patients)
	plot_1a.data$Mutations<-as.numeric(plot_1a.data$Mutations)
	plot_1a.data$Groups<-as.factor(plot_1a.data$Groups)

	x.scale <<- as.character(mutation_num_table$Patients)
	j=0
	k=0
	for(i in 1:length(mutation_num_table$Patients)){
	  if((i-1)%%2 != 0){
	    k=i-j
	    x.scale <<- x.scale[-k]
	    j<-j+1
	  }
	}

	orderID <<- c(1:nrow(plot_1a.data))


	F1a.plot<-ggplot()+theme_classic()
	F1a.plot<-F1a.plot+geom_bar(data=plot_1a.data,aes(x=reorder(Patients,orderID),y=Mutations,fill=Groups),width=0.7,color="black",stat='identity',position=position_stack())
	F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.4,2,0.4,2),'lines'),
	                      plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
						  legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.8,'cm'),legend.position=c(0.92,0.8),legend.text=element_text(size=16,face='bold.italic'),
						  axis.text.x=element_text(size=16,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold',color='black'),
						  axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.9,vjust=2,face='bold',color='black'))
	F1a.plot<-F1a.plot+ggtitle(NULL)+xlab(NULL)+ylab('Number of somatic mutations')+scale_fill_manual(name=NULL,values=c('black',gg_color_hue(2)[2],gg_color_hue(2)[1]),labels=c('Recurrence','Common','Initial'))
	F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,15000),breaks=seq(0,1500,100))+scale_x_discrete(breaks=x.scale)

	yaxis <- mutation_num_table[,2] + mutation_num_table[,3] + mutation_num_table[,4]
	maxy <- max(yaxis)

	if(maxy > 1000){
		split1 <- F1a.plot + coord_cartesian(ylim = c(0, 390)) +
	 	theme(legend.position='none',plot.margin=unit(c(0.2,0.6,0.1,2),'lines')) +ylab(NULL)
		split2 <- F1a.plot + coord_cartesian(ylim = c(1000, 1200)) + ggtitle(NULL)+
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',plot.margin=unit(c(1,0.6,0.2,2),'lines'))+ylab(NULL)
	}
	else if(maxy > 800){
		split1 <- F1a.plot + coord_cartesian(ylim = c(0, 390)) +
	 	theme(legend.position='none',plot.margin=unit(c(0.2,0.6,0.1,2),'lines')) +ylab(NULL)
		split2 <- F1a.plot + coord_cartesian(ylim = c(800, 1000)) + ggtitle(NULL)+
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',plot.margin=unit(c(1,0.6,0.2,2),'lines'))+ylab(NULL)
	}
	else if(maxy > 600){
		split1 <- F1a.plot + coord_cartesian(ylim = c(0, 390)) +
	 	theme(legend.position='none',plot.margin=unit(c(0.2,0.6,0.1,2),'lines')) +ylab(NULL)
		split2 <- F1a.plot + coord_cartesian(ylim = c(600, 800)) + ggtitle(NULL)+
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',plot.margin=unit(c(1,0.6,0.2,2),'lines'))+ylab(NULL)
	}
	else if(maxy > 400){
		split1 <- F1a.plot + coord_cartesian(ylim = c(0, 600)) +
	 	theme(legend.position='none',plot.margin=unit(c(1,0.6,0.1,2),'lines')) +ylab(NULL)
	}
	else if(maxy > 200){
		split1 <- F1a.plot + coord_cartesian(ylim = c(0, 400)) +
	 	theme(legend.position='none',plot.margin=unit(c(1,0.6,0.1,2),'lines')) +ylab(NULL)
	}
	else{
		split1 <- F1a.plot + coord_cartesian(ylim = c(0, 200)) +
	 	theme(legend.position='none',plot.margin=unit(c(1,0.6,0.1,2),'lines')) +ylab(NULL)
	}


	plot_1b.data <- melt(mutation_gene_table)


	F1b.plot<-ggplot()+theme_classic()
	F1b.plot<-F1b.plot+geom_tile(data = plot_1b.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
	F1b.plot<-F1b.plot+scale_fill_manual(name=NULL,values=c(N='white',C=gg_color_hue(2)[2],P=gg_color_hue(2)[1],R='black'),labels=c(N='None',C='Common',P='Primary',R='Recurrence'))
	F1b.plot<-F1b.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
	                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='bottom',
							              legend.direction="horizontal",legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
							            axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
	F1b.plot<-F1b.plot+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)


	gene.scale = 1 + (ncol(mutation_gene_table) - 12)/12

	if(maxy > 600){
		figure_1<-rbind(ggplotGrob(split2),ggplotGrob(split1),ggplotGrob(F1b.plot),size="last")
		panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
		figure_1$heights[panels][1] <- unit(0.513,'null')
		figure_1$heights[panels][3] <- unit(gene.scale,'null')
	}
	else{
		figure_1<-rbind(ggplotGrob(split1),ggplotGrob(F1b.plot),size="last")
		panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
		figure_1$heights[panels][3] <- unit(gene.scale,'null')
	}


	grid.draw(figure_1)
	ggsave(file="Figure1_mutation_landscape.pdf", plot=figure_1,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

	x.scale <<- NULL
	orderID <<- NULL
}


#function 4: Co-mutation analysis
mutCorrelation <- function(mutation_gene_table){
	mutGene.P <- mutation_gene_table
	mutGene.P[mutGene.P == 'C' | mutGene.P == 'P'] <- 1
	mutGene.P[mutGene.P == 'R' | mutGene.P == 'N'] <- 0

	mutGene.R <- mutation_gene_table
	mutGene.R[mutGene.R == 'C' | mutGene.R == 'R'] <- 1
	mutGene.R[mutGene.R == 'P' | mutGene.R == 'N'] <- 0

	plot_2.data <- data.frame( rep(0,ncol(mutation_gene_table)^2),rep(0,ncol(mutation_gene_table)^2),rep(0,ncol(mutation_gene_table)^2),rep(0,ncol(mutation_gene_table)^2))
	n=1
	for(i in 1:ncol(mutation_gene_table)){
		for(j in 1:ncol(mutation_gene_table)){
			plot_2.data[n,1] <- colnames(mutation_gene_table)[i]
			plot_2.data[n,2] <- colnames(mutation_gene_table)[j]

			Mut.FEtest <- cbind(c(0,0),c(0,0))
			if(i > j){ #Primary Tumor
				Mut.FEtest[1,1] <- length(which(mutGene.P[,i] == 1 & mutGene.P[,j] == 1))
				Mut.FEtest[1,2] <- length(which(mutGene.P[,i] == 1 & mutGene.P[,j] == 0))
				Mut.FEtest[2,1] <- length(which(mutGene.P[,i] == 0 & mutGene.P[,j] == 1))
				Mut.FEtest[2,2] <- length(which(mutGene.P[,i] == 0 & mutGene.P[,j] == 0))
				coMutation <- ((Mut.FEtest[1,1]+1)*(Mut.FEtest[2,2]+1)) / ((Mut.FEtest[1,2]+1)*(Mut.FEtest[2,1]+1))

				#fisher exact test
				pValue <- fisher.test(Mut.FEtest,alternative ="two.sided")$p.value

				if(pValue < 0.1){
					plot_2.data[n,3] <- -log10(pValue) + 1 #for dot size, but if pValue > 0.1 , filling with an instinct number
					if(coMutation > 1){
						plot_2.data[n,4] <- 'D_red'
					}
					else{
						plot_2.data[n,4] <- 'A_blue'
					}
				}
				else{
					plot_2.data[n,3] <- 1 #for dot size, but if pValue > 0.1 , filling with an instinct number
					plot_2.data[n,4] <- 'C_grey'
				}
			}

			else if(i < j){
				Mut.FEtest[1,1] <- length(which(mutGene.R[,i] == 1 & mutGene.R[,j] == 1))
				Mut.FEtest[1,2] <- length(which(mutGene.R[,i] == 1 & mutGene.R[,j] == 0))
				Mut.FEtest[2,1] <- length(which(mutGene.R[,i] == 0 & mutGene.R[,j] == 1))
				Mut.FEtest[2,2] <- length(which(mutGene.R[,i] == 0 & mutGene.R[,j] == 0))
				coMutation <- ((Mut.FEtest[1,1]+1)*(Mut.FEtest[2,2]+1)) / ((Mut.FEtest[1,2]+1)*(Mut.FEtest[2,1]+1))

				#fisher exact test
				pValue <- fisher.test(Mut.FEtest,alternative ="two.sided")$p.value

				if(pValue < 0.1){
					plot_2.data[n,3] <- -log10(pValue) + 1#for dot size, but if pValue > 0.1 , filling with an instinct number
					if(coMutation > 1){
						plot_2.data[n,4] <- 'E_black'
					}
					else{
						plot_2.data[n,4] <- 'B_green'
					}
				}
				else{
					plot_2.data[n,3] <- 1 #for dot size, but if pValue > 0.1 , filling with an instinct number
					plot_2.data[n,4] <- 'C_grey'
				}
			}

			else{  # i==j, filled with 'NA'
				plot_2.data[n,3] <- NA #for dot size, but if pValue > 0.1 , filling with an instinct number
				plot_2.data[n,4] <- NA #red for Co-mutation in P; green stand for exclusive in P,black for Co-mutation in P; blue stand for exclusive in R; grey for not significant in both P and R
			}
		n <- n+1
		}
	}

	#Figuring
	colnames(plot_2.data) <- c('listA','listB','dotSize','dotColor')

	orderID <<- c(1:nrow(plot_2.data))

	coMut.plot<-ggplot()+theme_classic()
	coMut.plot<-coMut.plot+geom_point(data = plot_2.data,aes(x=reorder(listA,orderID),y=reorder(listB,orderID),size=dotSize,color=as.factor(dotColor)))+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
	coMut.plot<-coMut.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.7,0.7,0.7,0.7),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
	                text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='bottom',legend.direction="horizontal",
							   legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
							   axis.text.x=element_text(size=14,angle=90,vjust=0.5,hjust=1,face='bold.italic',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),
							   axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))

	coMut.plot<-coMut.plot+scale_colour_manual(name=NULL,values=c(E_black='black',A_blue='blue',C_grey='grey',B_green='green',D_red='red'),
                                         labels=c(E_black='Co-occurrence in recurrent',A_blue='Exclusive in primary',C_grey='Not significant',B_green='Exclusive in recurrent',D_red='Co-occurrence in primary'),
                                         guide = guide_legend(override.aes=list(size=4),nrow=3),na.translate = F)+guides(size=FALSE)


	figure_2<-rbind(ggplotGrob(coMut.plot),size="first")
	grid.draw(figure_2)
	ggsave(file="Figure2_CoMutation.pdf", plot=figure_2,bg = 'white', width = 20, height = 22, units = 'cm', dpi = 600)

	x.scale <<- NULL
	orderID <<- NULL
}


#function 5: 3D-bubble plot of mutation frequency
mutFrequency<- function(savi_table, selected_geneList, mutation_gene_table, freq_cutoff){

	savi <- savi_table
	cutoff <- freq_cutoff
	savi<-savi[which(!(savi$Effect_Impact %in% 'LOW')),]  ## exclude the synonymous_variant

	selGene <- selected_geneList
	case <- unique(savi$CaseID)

	Mut.freq <- cbind(rep(0,length(selGene)),rep(0,length(selGene)),rep(0,length(selGene)))
	for(i in 1:length(selGene)){
	  #savi.Sel <- savi[which(savi$Gene_Name == selGene[i]),]
	  savi.Sel <- savi[which(savi$Gene_Name == selGene[i]),]
	  temp.P <- savi.Sel$Primary_freq
	  temp.R <- savi.Sel$Recurrent_freq
	  Mut.freq[i,1] <- length(unique(savi.Sel$CaseID[which( (savi.Sel$Primary_freq >= cutoff) & (savi.Sel$Recurrent_freq < cutoff))]) )
	  Mut.freq[i,2] <- length(unique(savi.Sel$CaseID[which( (savi.Sel$Primary_freq < cutoff) & (savi.Sel$Recurrent_freq >= cutoff))]) )
	  Mut.freq[i,3] <- length( which(mutation_gene_table[,i] ==  'C' ))
	  #Mut.freq[i,2] <- length(which( (temp.P < cutoff) & (temp.R >= cutoff) ))
	  #Mut.freq[i,3] <- length(which( (temp.P >= cutoff) & (temp.R >= cutoff) ))
	}

	colnames(Mut.freq) <- c('Primary','Recurent','Common')
	rownames(Mut.freq) <- selGene

	#build the fake 3D axis
	freq.plot<-ggplot()+theme_classic()
	freq.plot<-freq.plot+geom_segment(aes(x = 0, xend = 0, y = 0, yend = 45), arrow = arrow(), colour = gg_color_hue(2)[2], size = 1)
	freq.plot<-freq.plot+geom_segment(aes(x = 0, xend = cospi(45/180)*25, y = 0, yend = -(sinpi(45/180)*25)), arrow = arrow(), colour = 'black', size = 1)
	freq.plot<-freq.plot+geom_segment(aes(x = 0, xend = -cospi(45/180)*25, y = 0, yend = -(sinpi(45/180)*25)), arrow = arrow(), colour = gg_color_hue(2)[1], size = 1)
	freq.plot<-freq.plot+geom_text(aes(x = 0, y = 45,label = 'Common'), colour = gg_color_hue(2)[2], size = 5,vjust=-0.7,hjust=0.5)
	freq.plot<-freq.plot+geom_text(aes(x = cospi(45/180)*25, y = -(sinpi(45/180)*25),label = 'Recurrence'), colour = 'black', size = 5,vjust=1.2,hjust=0.5)
	freq.plot<-freq.plot+geom_text(aes(x = -cospi(45/180)*25, y = -(sinpi(45/180)*25),label = 'Primary'), colour = gg_color_hue(2)[1], size = 5,vjust=1.1,hjust=0.5)

	ng <-length(selGene)
	color1 <- c(255/255,109/255,96/255)
	color2 <- c(0/255,0/255,0/255)
	color3 <- c(0/255,197/255,204/255)
	x.position <- rep(0,ng)
	y.position <- rep(0,ng)
	sample.size <- rep(0,ng)
	total <- rep(0,ng)
	temp <- cbind(rep(0,ng),rep(0,ng),rep(0,ng))
	color.value <- rep(0,ng)
	for(i in 1:length(selGene)){
		x.position[i] <- Mut.freq[i,2] - Mut.freq[i,1]
		y.position[i] <- 2*Mut.freq[i,3] - Mut.freq[i,2] - Mut.freq[i,1]
		total[i] <- Mut.freq[i,3] + Mut.freq[i,2] + Mut.freq[i,1]
		sample.size[i] <- (Mut.freq[i,1] + Mut.freq[i,2] + Mut.freq[i,3])/3 +5
		temp[i,] <- c(Mut.freq[i,1]*1.1/total[i] * color1 + Mut.freq[i,2]*1.1/total[i] * color2 + Mut.freq[i,3]*1.1/total[i] * color3)
		color.value[i] <- rgb(temp[i,1],temp[i,2],temp[i,3])
	}

	plot_3.data<-data.frame(x.position,y.position,selGene,color.value,sample.size)
	colnames(plot_3.data) <- c('x.position','y.position','Gene_Name','color.value','sample.size')

	freq.plot<-freq.plot+geom_point(data=plot_3.data,aes(x =x.position, y = y.position), size = sample.size/2,colour = color.value,alpha=0.7)
	for(i in 1:length(x.position)){
		if( (x.position[i]^2 + y.position[i]^2 > 50 && y.position[i] > 1) || (selGene[i] == 'LTBP4') || (sample.size[i] > 10) ){
		  if(x.position[i] >= 0 ){
			  freq.plot<-freq.plot+geom_text(data=plot_3.data,x =x.position[i], y = y.position[i],label = plot_3.data$Gene_Name[i], size = sample.size[i]/3, colour = 'black',vjust = 0.1, hjust=-0.25,fontface = "italic")
		  }
		  else{
		    freq.plot<-freq.plot+geom_text(data=plot_3.data,x =x.position[i], y = y.position[i],label = plot_3.data$Gene_Name[i], size = sample.size[i]/3, colour = 'black',vjust = 0.1, hjust=1.25,fontface = "italic")
		  }
		}
	}
	freq.plot<-freq.plot+theme(panel.background=element_blank(),panel.border= element_blank(),axis.line = element_blank(),title=element_blank(),text=element_text(size=16,face='bold.italic'),
	                           legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position=c(0.2,0.87),legend.direction="vertical",legend.text=element_text(size=16,face='bold.italic'),
							               axis.ticks = element_blank(), axis.text = element_blank())
	freq.plot<-freq.plot+guides(size=FALSE,color=FALSE)+ylab(NULL)+xlab(NULL)+scale_y_continuous(expand=c(0,1),limits=c(-25,50))+scale_x_continuous(expand=c(0,1),limits=c(-30,30))


	figure_3 <- rbind(ggplotGrob(freq.plot),size="first")
	grid.draw(figure_3)

	ggsave(file="Figure3_3D_bubble_plot.pdf", plot=figure_3,bg = 'white', width = 28, height = 20, units = 'cm', dpi = 600)

	return(Mut.freq)
}


#function 6: hypermutation analysis
mutSignature <- function(savi_table, mut_freq_cutoff, mut_num_cutoff, HM_score_cutoff){

	varPresentcut <- mut_freq_cutoff #Please be aware that the threshold here is different from that of figure 1 (to filter noise)
	cutMutLoad <<- mut_num_cutoff
	cutHMscore <<- HM_score_cutoff
	savi <- savi_table
	rownames(savi) = c(1:nrow(savi))
	case <- unique(savi$CaseID)

	Mut.Load.P <- rep(0,length(case))
	Mut.Load.R <- rep(0,length(case))
	HM.Score.P <- rep(0,length(case))
	HM.Score.R <- rep(0,length(case))
	savi.P <- rep(0,length(savi[,1]))
	savi.R <- rep(0,length(savi[,1]))
	savi.isHM <- rep(0,length(savi[,1]))
	
	savi$isC2T <- ((savi$ref == 'C' & savi$alt == 'T') | (savi$ref == 'G' & savi$alt == 'A'))
	savi$isCC2TC <- ((savi$ref == 'C' & savi$alt == 'T' & savi$varSurffix == 'C') | (savi$ref == 'G' & savi$alt == 'A' & savi$varPrefix == 'G'))
	savi$isCT2TT <- ((savi$ref == 'C' & savi$alt == 'T' & savi$varSurffix == 'T') | (savi$ref == 'G' & savi$alt == 'A' & savi$varPrefix == 'A'))
	#savi <- cbind(savi,isC2T,isCC2TC,isCT2TT)
	
	p0 = 1
		
	for(i in 1:length(case)){
		#for each patient
		saviPatient <- savi[which(savi$CaseID == case[i]),]

		#for primary
		Pvar <- saviPatient[which((saviPatient$Blood_freq == 0) & (saviPatient$Primary_freq > varPresentcut)),]
		nPvar <- length(Pvar[,1])
		savi.P[as.integer(rownames(Pvar))] <- 1
		Mut.Load.P[i] <- nPvar

		#nPctga <- length( which((Pvar$ref == 'C' & Pvar$alt == 'T') | (Pvar$ref == 'G' & Pvar$alt == 'A') ) )
		#nPccgg <- length( which((Pvar$ref == 'C' & Pvar$varSurffix == 'C') | (Pvar$ref == 'G' & Pvar$varPrefix == 'G') ) )
		#nPcggc <- length( which((Pvar$ref == 'C' & Pvar$varSurffix == 'G') | (Pvar$ref == 'G' & Pvar$varPrefix == 'C') ) )
    
		pC2T <- sum(Pvar$isC2T) / (Mut.Load.P[i] + p0 )
		pCC2TC <- sum(Pvar$isCC2TC) / (sum(Pvar$isC2T) + p0)
		pCT2TT <- sum(Pvar$isCT2TT) / (sum(Pvar$isC2T) + p0)
		
		#cat(paste0(sum(Pvar$isC2T[which(savi.P==1)]),"\t\t\t\t", (Mut.Load.P[i] + p0 )," \t\t\t"," ", pCC2TC,"\t","3", pCT2TT,"\n"))
		
		#HM.Score.P[i] = (nPccgg - nPcggc)/(nPctga + 1) + nPctga/(nPvar + 1)
		HM.Score.P[i] = pC2T + pCC2TC + sign(pCC2TC - pCT2TT) * pCT2TT 

		#for recurrence
		Rvar <- saviPatient[which((saviPatient$Blood_freq == 0) & (saviPatient$Recurrent_freq > varPresentcut)),]
		nRvar <- length(Rvar[,1])
		savi.R[as.integer(rownames(Rvar))] <- 1
		Mut.Load.R[i] <- nRvar

		#nRctga <- length( which((Rvar$ref == 'C' & Rvar$alt == 'T') | (Rvar$ref == 'G' & Rvar$alt == 'A') ) )
		#nRccgg <- length( which((Rvar$ref == 'C' & Rvar$varSurffix == 'C') | (Rvar$ref == 'G' & Rvar$varPrefix == 'G') ) )
		#nRcggc <- length( which((Rvar$ref == 'C' & Rvar$varSurffix == 'G') | (Rvar$ref == 'G' & Rvar$varPrefix == 'C') ) )
		
		rC2T <- sum(Rvar$isC2T) / (Mut.Load.R[i] + p0 )
		rCC2TC <- sum(Rvar$isCC2TC) / (sum(Rvar$isC2T) + p0)
		rCT2TT <- sum(Rvar$isCT2TT) / (sum(Rvar$isC2T) + p0)

		#HM.Score.R[i] = (nRccgg - nRcggc)/(nRctga + 1) + nRctga/(nRvar + 1)
		HM.Score.R[i] = rC2T + rCC2TC + sign(rCC2TC - rCT2TT) * rCT2TT
	}

	HM.mark.P <- rep(0,length(case))
	HM.mark.R <- rep(0,length(case))
	for(i in 1:length(case)){
		if(Mut.Load.P[i] > cutMutLoad & HM.Score.P[i] > cutHMscore){
			savi.isHM[which(savi$CaseID == case[i])] <- 1
			HM.mark.P[i] <- 1
		}

		if(Mut.Load.R[i] > cutMutLoad & HM.Score.R[i] > cutHMscore){
			savi.isHM[which(savi$CaseID == case[i])] <- 1
			HM.mark.R[i] <- 1
		}
	}

	plot_4a.data <- data.frame(c(as.character(case),as.character(case)),c(Mut.Load.P+1,Mut.Load.R+1),c(HM.Score.P,HM.Score.R),c(HM.mark.P,HM.mark.R),c(rep('P',length(case)),rep('R',length(case))))
	colnames(plot_4a.data) <- c('caseID','mutNumber','HMscore','HMmark','PRmark')

	#figure 4a (hypermuation detection)
	HM.plot<-ggplot()+theme_classic()
	HM.plot<-HM.plot+geom_point(data = plot_4a.data,aes(x=mutNumber,y=HMscore,fill=as.factor(PRmark),shape=as.factor(HMmark)),color="black",size=3,alpha=0.7,stroke=0.5)+ylab('HM score')+xlab('Number of mutations')+ggtitle(NULL)
	HM.plot<-HM.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
						   text=element_text(size=2,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='top',legend.direction="horizontal",
						   legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,face='bold',color='black'),
						   axis.text.x=element_text(size=14,face='bold',color='black'),axis.title.x=element_text(size=18,face='bold',color='black'),axis.title.y=element_text(size=18,hjust=0.5,vjust=2,face='bold',color='black'))
	HM.plot<-HM.plot+geom_hline(aes(yintercept=cutHMscore),color='black',size=.5,linetype='dashed')
	HM.plot<-HM.plot+geom_vline(aes(xintercept=cutMutLoad),color='black',size=.5,linetype='dashed')
	HM.plot<-HM.plot+scale_x_continuous(trans="log10",expand=c(0,0.1))
	HM.plot<-HM.plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
	#HM.plot<-HM.plot+scale_color_manual(name=NULL,values=c(R='black',P=gg_color_hue(2)[1]),labels=c(R='Recurrence',P='Primary'),guide = guide_legend(override.aes=list(size=4),nrow=2))+guides(size=FALSE)
	HM.plot<-HM.plot+scale_fill_manual(name=NULL,values=c(R='black',P=gg_color_hue(2)[1]),labels=c(R='Recurrence',P='Primary'),guide = guide_legend(override.aes=list(size=4,shape=21),nrow=2))+guides(size=FALSE)
	HM.plot<-HM.plot+scale_shape_manual(name=NULL,values=c(22,24),labels=c('1'='Hypermutation','0'='Non-hypermuation'),guide = guide_legend(override.aes=list(size=4),keywidth=0.1,keyheight=0.2,nrow=2))
	#HM.plot<-HM.plot+scale_y_continuous(expand=c(0,0.1),limits=c(-0.5,1.8))
	HM.plot<-HM.plot+scale_y_continuous(expand=c(0,0.1))
	#figure_4a<-rbind(ggplotGrob(HM.plot),size="first")
	#grid.draw(figure_4a)
	#ggsave(file="hypermutation_detection.pdf", plot=figure_4a,bg = 'white', width = 18, height = 18, units = 'cm', dpi = 600)


	#nucleotides changes in all mutations
	subTypes.P <- savi[which(savi.P == 1),]
	subTypes.R <- savi[which(savi.isHM == 0 & savi.R == 1),]
	subTypes.HR <- savi[which(savi.isHM == 1 & savi.R == 1),]

	#for primary
	Mut.CTGA.P <- length( which((subTypes.P$ref == 'C' & subTypes.P$alt == 'T') | (subTypes.P$ref == 'G' & subTypes.P$alt == 'A') ) )
	Mut.CGGC.P <- length( which((subTypes.P$ref == 'C' & subTypes.P$alt == 'G') | (subTypes.P$ref == 'G' & subTypes.P$alt == 'C') ) )
	Mut.CAGT.P <- length( which((subTypes.P$ref == 'C' & subTypes.P$alt == 'A') | (subTypes.P$ref == 'G' & subTypes.P$alt == 'T') ) )
	Mut.ATTA.P <- length( which((subTypes.P$ref == 'A' & subTypes.P$alt == 'T') | (subTypes.P$ref == 'T' & subTypes.P$alt == 'A') ) )
	Mut.AGTC.P <- length( which((subTypes.P$ref == 'A' & subTypes.P$alt == 'G') | (subTypes.P$ref == 'T' & subTypes.P$alt == 'C') ) )
	Mut.ACTG.P <- length( which((subTypes.P$ref == 'A' & subTypes.P$alt == 'C') | (subTypes.P$ref == 'T' & subTypes.P$alt == 'G') ) )
	sumMut.P <- Mut.CTGA.P+Mut.CGGC.P+Mut.CAGT.P+Mut.ATTA.P+Mut.AGTC.P+Mut.ACTG.P

	#for recurrence without hypermutation
	Mut.CTGA.R <- length( which((subTypes.R$ref == 'C' & subTypes.R$alt == 'T') | (subTypes.R$ref == 'G' & subTypes.R$alt == 'A') ) )
	Mut.CGGC.R <- length( which((subTypes.R$ref == 'C' & subTypes.R$alt == 'G') | (subTypes.R$ref == 'G' & subTypes.R$alt == 'C') ) )
	Mut.CAGT.R <- length( which((subTypes.R$ref == 'C' & subTypes.R$alt == 'A') | (subTypes.R$ref == 'G' & subTypes.R$alt == 'T') ) )
	Mut.ATTA.R <- length( which((subTypes.R$ref == 'A' & subTypes.R$alt == 'T') | (subTypes.R$ref == 'T' & subTypes.R$alt == 'A') ) )
	Mut.AGTC.R <- length( which((subTypes.R$ref == 'A' & subTypes.R$alt == 'G') | (subTypes.R$ref == 'T' & subTypes.R$alt == 'C') ) )
	Mut.ACTG.R <- length( which((subTypes.R$ref == 'A' & subTypes.R$alt == 'C') | (subTypes.R$ref == 'T' & subTypes.R$alt == 'G') ) )
	sumMut.R <- Mut.CTGA.R+Mut.CGGC.R+Mut.CAGT.R+Mut.ATTA.R+Mut.AGTC.R+Mut.ACTG.R

	#for recurrence with hypermutation
	Mut.CTGA.HR <- length( which((subTypes.HR$ref == 'C' & subTypes.HR$alt == 'T') | (subTypes.HR$ref == 'G' & subTypes.HR$alt == 'A') ) )
	Mut.CGGC.HR <- length( which((subTypes.HR$ref == 'C' & subTypes.HR$alt == 'G') | (subTypes.HR$ref == 'G' & subTypes.HR$alt == 'C') ) )
	Mut.CAGT.HR <- length( which((subTypes.HR$ref == 'C' & subTypes.HR$alt == 'A') | (subTypes.HR$ref == 'G' & subTypes.HR$alt == 'T') ) )
	Mut.ATTA.HR <- length( which((subTypes.HR$ref == 'A' & subTypes.HR$alt == 'T') | (subTypes.HR$ref == 'T' & subTypes.HR$alt == 'A') ) )
	Mut.AGTC.HR <- length( which((subTypes.HR$ref == 'A' & subTypes.HR$alt == 'G') | (subTypes.HR$ref == 'T' & subTypes.HR$alt == 'C') ) )
	Mut.ACTG.HR <- length( which((subTypes.HR$ref == 'A' & subTypes.HR$alt == 'C') | (subTypes.HR$ref == 'T' & subTypes.HR$alt == 'G') ) )
	sumMut.HR <- Mut.CTGA.HR+Mut.CGGC.HR+Mut.CAGT.HR+Mut.ATTA.HR+Mut.AGTC.HR+Mut.ACTG.HR

	bases <- c('CTGA', 'CGGC', 'CAGT', 'ATTA', 'AGTC', 'ACTG')
	Fraction <- rep(0,18)
	i<-1
	for(x in bases){
	  filenameP<-paste0('Mut.',x,'.P')
	  filenameR<-paste0('Mut.',x,'.R')
	  filenameHR<-paste0('Mut.',x,'.HR')
	  Fraction[i] <- get(filenameP)/sumMut.P
	  Fraction[i+6] <- get(filenameR)/sumMut.R
	  Fraction[i+12] <- get(filenameHR)/sumMut.HR
	  i <- i+1
	}

	plot_4b.data <- data.frame(c(rep('Primary',6),rep('Rec_noHM',6),rep('Rec_HM',6)), Fraction,rep(c('a','b','c','d','e','f'),3))
	colnames(plot_4b.data) <- c('sampleType','Fraction','MuType')

	#figure 4b
	orderID_4b<<-c(1:nrow(plot_4b.data))
	Fra.plot<-ggplot()+theme_classic()
	Fra.plot<-Fra.plot+geom_bar(data=plot_4b.data,aes(x=reorder(sampleType,orderID_4b),y=Fraction*100,fill=MuType),width=0.4,color="black",stat='identity',alpha=0.8,position=position_stack())+coord_cartesian(ylim = c(0, 100))
	Fra.plot<-Fra.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
							 legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
							 legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
							 axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
							 axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=18,face='bold',color='black'))

	Fra.plot<-Fra.plot+ggtitle(NULL)+xlab(NULL)+ylab('Mutation fraction (%)')+scale_fill_manual(name=NULL,values=c('black',gg_color_hue(5)),
							labels=c(a='C>T/G>A ',b='C>G/G>C ',c='C>A/G>T ',d='A>T/T>A ',e='A>G/T>C ',f='A>C/T>G '),guide = guide_legend(override.aes=list(size=0.5),nrow=2))
	Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Primary='Primary',Rec_noHM='Non-HM\nrecurrence',Rec_HM='HM\nrecurrence'))

	#figure_4b<-cbind(ggplotGrob(HM.plot),ggplotGrob(Fra.plot),size="first")
	#grid.draw(figure_4b)
	#ggsave(file="hypermutation_detection_fraction.pdf", plot=figure_4b,bg = 'white', width = 32, height = 17, units = 'cm', dpi = 600)

	# silent/Missense ratio
	Ratio.P <- rep(NA,length(case))
	Ratio.R <- rep(NA,length(case))
	Ratio.HR <- rep(NA,length(case))

	for(i in 1:length(case)){   #for each patient
		tmp.P <- savi[which(savi$CaseID == case[i] & savi.P == 1),]
		tmp.R <- savi[which(savi$CaseID == case[i] & savi.R == 1 & savi.isHM == 0),]
		tmp.HR <- savi[which(savi$CaseID == case[i] & savi.R == 1 & savi.isHM == 1),]

		Missense.P <- length(which(grepl(pattern = 'MISSENSE',x = tmp.P$Functional_Class,fixed=T)))
		Silent.P <- length(which(grepl(pattern = 'SILENT',x = tmp.P$Functional_Class,fixed=T)))
		Missense.R <- length(which(grepl(pattern = 'MISSENSE',x = tmp.R$Functional_Class,fixed=T)))
		Silent.R <- length(which(grepl(pattern = 'SILENT',x = tmp.R$Functional_Class,fixed=T)))
		Missense.HR <- length(which(grepl(pattern = 'MISSENSE',x = tmp.HR$Functional_Class,fixed=T)))
		Silent.HR <- length(which(grepl(pattern = 'SILENT',x = tmp.HR$Functional_Class,fixed=T)))

		#filter the samples that mutations fewer than 10
		if(length(tmp.P[,1]) > 10 & Missense.P > 0){
			Ratio.P[i] <- Silent.P / Missense.P
		}
		if(length(tmp.R[,1]) > 10 & Missense.R > 0){
			Ratio.R[i] <- Silent.R / Missense.R
		}
		if(length(tmp.HR[,1]) > 10 & Missense.HR > 0){
			Ratio.HR[i] <- Silent.HR / Missense.HR
		}
	}

	plot_4c.data <- data.frame( c(rep('Primary',length(Ratio.P)),rep('Rec_nonHM',length(Ratio.R)),rep('Rec_HM',length(Ratio.HR))),c(Ratio.P,Ratio.R,Ratio.HR),c(rep('P',length(Ratio.P)),rep('R',length(Ratio.R)),rep('R',length(Ratio.HR))) )
	plot_4c.data <- na.omit(plot_4c.data)
	colnames(plot_4c.data) <- c('Groups','Ratio','PRmark')

	#Wilcoxson rank sum testing
	if(length(which(plot_4c.data[,1] == 'Rec_HM')) > 0 & length(which(plot_4c.data[,1] == 'Primary')) > 0 & length(which(plot_4c.data[,1] == 'Rec_nonHM')) > 0 ){
		Pvalue.P_R <- wilcox.test(Ratio ~ Groups,data = plot_4c.data,subset = Groups %in% c('Primary','Rec_nonHM'),alternative = "two.sided",paired = F,conf.level = 0.05)$p.value
		Pvalue.P_HR <- wilcox.test(Ratio ~ Groups,data = plot_4c.data,subset = Groups %in% c('Primary','Rec_HM'),alternative = "two.sided",paired = F,conf.level = 0.05)$p.value
		Pvalue.R_HR <- wilcox.test(Ratio ~ Groups,data = plot_4c.data,subset = Groups %in% c('Rec_nonHM','Rec_HM'),alternative = "two.sided",paired = F,conf.level = 0.05)$p.value
		cat(paste0("P-value between Primary and NonHM Recurrence: ", Pvalue.P_R,"\n","P-value between Primary and HM Recurrence: ", Pvalue.P_HR,"\n","P-value between NonHM Recurrence and HM Recurrence: ", Pvalue.R_HR,"\n"))
	}

	#figure 2a (hypermuation detection)
	orderID<<-c(1:nrow(plot_4c.data))
	Box.plot<-ggplot()+theme_classic()
	Box.plot<-Box.plot+stat_boxplot(data = plot_4c.data,geom='errorbar',width=0.3,size=1,aes(x=reorder(Groups,orderID),y=Ratio,color=as.factor(PRmark)))
	Box.plot<-Box.plot+geom_boxplot(data = plot_4c.data,aes(x=reorder(Groups,orderID),y=Ratio,color=as.factor(PRmark)),fill='white',size=1,width=0.4,alpha=1,linetype=1)+ylab('Silent/missense ratio')+xlab(NULL)+ggtitle(NULL)
	Box.plot<-Box.plot+geom_boxplot(data = plot_4c.data,aes(x=reorder(Groups,orderID),y=Ratio,color=as.factor(PRmark),fill=as.factor(PRmark)),size=1,width=0.4,alpha=0.4,linetype=1)+ylab('Silent/missense ratio')+xlab(NULL)+ggtitle(NULL)
	Box.plot<-Box.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
	                        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
	                        legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
	                        axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
	                        axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=18,face='bold',color='black'))

	Box.plot<-Box.plot+scale_colour_manual(name=NULL,values=c(R='black',P=gg_color_hue(2)[1]),labels=c(R='Recurrence',P='Primary'),na.translate = F,guide = guide_legend(nrow=1))+guides(size=FALSE)
	Box.plot<-Box.plot+scale_fill_manual(name=NULL,values=c(R='black',P=gg_color_hue(2)[1]),na.translate = F,guide=FALSE)
	Box.plot<-Box.plot+scale_y_continuous(expand=c(0,0),limit=c(0,1.01),breaks=seq(0,1,0.2))+scale_x_discrete(labels=c(Primary='Primary',Rec_nonHM='Non-HM\nrecurrence',Rec_HM='HM\nrecurrence'))
	Box.plot<-Box.plot+stat_summary(data = plot_4c.data,fun.y=mean,geom='point',size=3,aes(x=reorder(Groups,orderID),y=Ratio),fill='white',shape=23)

	figure_4<-cbind(ggplotGrob(HM.plot),ggplotGrob(Fra.plot),ggplotGrob(Box.plot),size="first")
	grid.draw(figure_4)
	ggsave(file="Figure4_hypermutation_analysis.pdf", plot=figure_4,bg = 'white', width = 48, height = 17, units = 'cm', dpi = 600)

	cutMutLoad <<- NULL
	cutHMscore <<- NULL
	orderID_4b <<- NULL
	orderID <<- NULL
	return(plot_4c.data)
}


#function 7: plot the evoluationary clusters
mutTreeClustering <- function(mutation_num_table){

	tmp.sum <- mutation_num_table$Primary + mutation_num_table$Recurrent + mutation_num_table$Common
	plot_5.data <- data.frame( c(mutation_num_table$Primary/tmp.sum) ,c(mutation_num_table$Recurrent/tmp.sum),c(mutation_num_table$Common/tmp.sum))
	clustering <- kmeans(plot_5.data,3,iter.max=1000)
	clusters <- fitted(clustering,method = c( "classes"))
	#clusters[which(clusters == 1)] <- 'a'
	#clusters[which(clusters == 2)] <- 'b'
	#clusters[which(clusters == 3)] <- 'c'
	colnames(plot_5.data) <- c('Primary','Recurrence','Common')
	Result_cluster <- clusters
	Result_cluster[which(clusters == clusters[which(plot_5.data$Primary == max(plot_5.data$Primary))])] <- 'b_Primary'
	Result_cluster[which(clusters == clusters[which(plot_5.data$Recurrence == max(plot_5.data$Recurrence))])] <- 'a_Recurrence'
	Result_cluster[which(clusters == clusters[which(plot_5.data$Common == max(plot_5.data$Common))])] <- 'c_Common'
	plot_5.data <- cbind(plot_5.data,Result_cluster)
	colnames(plot_5.data) <- c('Primary','Recurrence','Common','Cluster')
	plot_5.data$id <- 1:nrow(plot_5.data)

  require(ggtern)
  require(ggalt)

	Ternary.plot <- ggtern(data=plot_5.data,aes(Primary,Common,Recurrence,color=Cluster,fill=Cluster))+theme_bw()
	Ternary.plot <- Ternary.plot + geom_encircle(alpha=0.2,size=1)+geom_point(size=4.5,alpha=0.6)+guides(fill=FALSE)
	Ternary.plot <- Ternary.plot + scale_colour_manual(name=NULL,values=c('black',gg_color_hue(2)[1],gg_color_hue(2)[2]),labels=c('Cluster 1','Cluster 2','Cluster 3'))
	Ternary.plot <- Ternary.plot + scale_fill_manual(name=NULL,values=c('black',gg_color_hue(2)[1],gg_color_hue(2)[2]))
	Ternary.plot <- Ternary.plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
	       text=element_text(size=13,face='bold',color='black'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position=c(0.85,0.8),legend.direction="vertical",legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
	       legend.text=element_text(size=16,face='bold.italic'))#+geom_text(aes(label = id),size = 3,color='black')
	Ternary.plot <- Ternary.plot+theme_showgrid() +theme_showarrows()+tern_limits(T = 1.02, L = 1.02, R = 1.02, verbose = F)

	figure5<-rbind(ggplotGrob(Ternary.plot),size="last")
	grid.draw(figure5)
	ggsave(file="Figure5_Case_Clustering.pdf", plot=figure5,bg = 'white', width = 24, height = 20, units = 'cm', dpi = 600)

  #detach("package:ggalt", unload=TRUE)
  #detach("package:ggtern", unload=TRUE)
  unloadNamespace("ggalt")
  unloadNamespace("ggtern")

	return(plot_5.data)
}

#function 8: Clonal mutation replacement in selected genes
mutSwitch <- function(savi_table, selected_geneList, freq_cutoff_low,freq_cutoff_high){

	savi <- savi_table[which(!(savi_table$Effect_Impact %in% 'LOW')),]
	case <- unique(savi$CaseID)
	CaseGene <- c()
	mutation.types <- c()
	x.content <- c()
	y.content <- c()
	color.factor <- c()

	for(i in 1:length(case)){   #for each patient
		tmp.case <- savi[which(savi$CaseID == case[i] & !grepl(pattern = 'LOW',x = savi$Effect_Impact,fixed=T)),]

		for(j in 1:length(selected_geneList)){
			NameCom <- paste0(case[i],'-',selected_geneList[j])
			tmp.switch <- tmp.case[which(tmp.case$Gene_Name == selected_geneList[j]),]
			if( any(tmp.switch$Primary_freq <= freq_cutoff_low & tmp.switch$Recurrent_freq >= freq_cutoff_high) ){
				if( any(tmp.switch$Primary_freq >= freq_cutoff_high & tmp.switch$Recurrent_freq <= freq_cutoff_low)){
					for( k in c(which( (tmp.switch$Primary_freq <= freq_cutoff_low & tmp.switch$Recurrent_freq >= freq_cutoff_high) | (tmp.switch$Primary_freq >= freq_cutoff_high & tmp.switch$Recurrent_freq <= freq_cutoff_low) ) )){
						CaseGene <- c(CaseGene, NameCom,NameCom)
						AAC <- gsub('^ *| *$','',as.character(tmp.switch$Amino_Acid_Change[k]))
						if( AAC %in% mutation.types){
							temp <- paste0(AAC,'NA')
							mutation.types <- c(mutation.types,as.character(temp),as.character(temp))
						}
						else{
						  mutation.types <- c(mutation.types,AAC,AAC)
						}
						x.content <- c(x.content,'Primary','Recurrence')
						y.content <- c(y.content,tmp.switch$Primary_freq[k],tmp.switch$Recurrent_freq[k])
						if(tmp.switch$Primary_freq[k] > tmp.switch$Recurrent_freq[k]){
							color.factor <- c(color.factor,'P','P')
						}
						else{
							color.factor <- c(color.factor,'R','R')
						}
					}
				}
			}
		}
	}

	plot_6.data <- data.frame(CaseGene,mutation.types,x.content,y.content,color.factor,stringsAsFactors = F)
	if(nrow(plot_6.data) == 0){
		cat(paste0("No switch event detected in the given gene list!\n"))
		return(plot_6.data)
	}

	Switch.plot<-ggplot()+theme_classic()
	Switch.plot<-Switch.plot+geom_point(data = plot_6.data,aes(x=x.content,y=y.content,color=as.factor(color.factor),shape=as.factor(color.factor)),size=3,alpha=0.9)
	Switch.plot<-Switch.plot+geom_line(data = plot_6.data,aes(x=x.content,y=y.content,color=as.factor(color.factor),group=as.factor(mutation.types)),size=1)+ylab('Mutation frequency in tumor')+xlab(NULL)+ggtitle(NULL)

	Switch.plot<-Switch.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
	                       text=element_text(size=18,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='null',legend.direction="horizontal",
	                       legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,face='bold',color='black'),
	                       axis.text.x=element_text(size=14,face='bold',color='black'),axis.title.x=element_text(size=18,face='bold',color='black'),axis.title.y=element_text(size=18,hjust=0.5,vjust=2,face='bold',color='black'),
						   strip.text = element_text(size=18,face='bold.italic'),strip.background = element_rect(colour="transparent", fill="transparent"))
	Switch.plot<-Switch.plot+scale_colour_manual(name=NULL,values=c(R='black',P=gg_color_hue(2)[1]))+scale_shape_manual(name=NULL,values=c(1,5))
	text.data.P <- plot_6.data[which(plot_6.data$y.content >= freq_cutoff_high & plot_6.data$x.content == 'Primary'),]
	text.data.R <- plot_6.data[which(plot_6.data$y.content >= freq_cutoff_high & plot_6.data$x.content == 'Recurrence'),]

	Switch.plot<-Switch.plot+geom_text(data = text.data.P,aes(x=x.content,y=y.content,color=as.factor(color.factor),label=mutation.types),size=3.5,hjust = 1, vjust = 0, nudge_x = -0.05, lineheight = 0.2,fontface = "bold")
	Switch.plot<-Switch.plot+geom_text(data = text.data.R,aes(x=x.content,y=y.content,color=as.factor(color.factor),label=mutation.types),size=3.5,hjust = 0, vjust = 0, nudge_x = 0.1, lineheight = 0.2,fontface = "bold")

	Switch.plot<-Switch.plot+scale_y_continuous(expand=c(0,1),limits=c(0,103),breaks=seq(0,100,20))


	numPlots <- length(unique(plot_6.data$CaseGene))
	numCols <- ceiling(sqrt(numPlots))
	numRows <- ceiling(numPlots/numCols)
	Switch.plot<-Switch.plot+facet_wrap( ~ CaseGene,ncol=numCols )
	figure_6<-rbind(ggplotGrob(Switch.plot), size="last")
	grid.draw(figure_6)
	ggsave(file="Figure6_mutation_Switch.pdf", plot=figure_6, bg = 'white', width = 7*numCols, height = 7*numRows-1, units = 'cm', dpi = 600)

	return(plot_6.data)
}


#function 9: Infer mutation order in tumor evolutionary history
mutDirectedGraph <- function(mutation_gene_table){

	input.table <- mutation_gene_table
	selected_geneList <- as.character(colnames(input.table))

	#calculating TEDG edge table
	temp <- rep(0,length(selected_geneList))
	edge.matrix <- temp
	for( i in 2:length(selected_geneList)){
		edge.matrix <- cbind(edge.matrix,temp)
	}

	edge.table <- c('geneA','geneB','weight','label')

	end <- length(selected_geneList)-1
	for( i in 1:end){
		start <- i+1
		for( j in start:length(selected_geneList) ){

			edge.matrix[i,j] <- length( which( (input.table[,i] == 'C') & (input.table[,j] %in% c('P','R')) ) )
			edge.matrix[j,i] <- length( which( (input.table[,j] == 'C') & (input.table[,i] %in% c('P','R')) ) )
			labelA <- paste( rownames(input.table)[which( (input.table[,i] == 'C') & (input.table[,j] %in% c('P','R')) )], collapse = ";")
			labelB <- paste( rownames(input.table)[which( (input.table[,j] == 'C') & (input.table[,i] %in% c('P','R')) )], collapse = ";")

			if(edge.matrix[i,j] < edge.matrix[j,i]){
				edge.matrix[i,j] <- 0
				edge <- c(selected_geneList[j],selected_geneList[i],edge.matrix[j,i],labelB)
				edge.table <- rbind(edge.table,edge)
			}
			else if(edge.matrix[i,j] > edge.matrix[j,i]){
				edge.matrix[j,i] <- 0
				edge <- c(selected_geneList[i],selected_geneList[j],edge.matrix[i,j],labelA)
				edge.table <- rbind(edge.table,edge)
			}
			else{
				edge.matrix[i,j] <- 0
				edge.matrix[j,i] <- 0
			}

		}
	}

	rownames(edge.matrix) <- selected_geneList
	colnames(edge.matrix) <- selected_geneList

	edge.table <- edge.table[-1,]
	colnames(edge.table) <- c('geneA','geneB','weight','label')
	rownames(edge.table) <- c(1:nrow(edge.table))

	write.table(edge.table,"TEDGedge.txt",row.names = F,quote = F,sep = '\t')

	#calculating TEDG node table
	Mut.freq <- cbind(rep(0,length(selected_geneList)),rep(0,length(selected_geneList)),rep(0,length(selected_geneList)))
	for(i in 1:length(selected_geneList)){
	  Mut.freq[i,1] <- length(which(input.table[,i] == 'P'))
	  Mut.freq[i,2] <- length(which(input.table[,i] == 'R'))
	  Mut.freq[i,3] <- length(which(input.table[,i] == 'C')) * 2
	}

	sample.size <- rep(0,length(selected_geneList))
	for(i in 1:length(selected_geneList)){
	  sample.size[i] <- Mut.freq[i,1] + Mut.freq[i,2] + Mut.freq[i,3]
	}

	ins <- rep(0, length(selected_geneList))
	outs <- rep(0,length(selected_geneList))
	for(i in 1:length(selected_geneList)){
	  ins[i] <- length(which(edge.matrix[,i]>0))
	  outs[i] <- length(which(edge.matrix[i,]>0))
	}


	pcdf <- rep(1,length(selected_geneList))
	for(i in 1:length(pcdf)){
	  if(ins[i] < outs[i]){
	    #y = binocdf(x,N,p) computes (x,y) that follow binomial dist (N,p)
	    pcdf[i] <- binom.test(ins[i], ins[i]+outs[i], 0.5)$p.value
	  }
	  else{
	    pcdf[i] <- binom.test(outs[i], ins[i]+outs[i], 0.5)$p.value
	  }
	}

	fc = log2((outs+1) / (ins+1)) # positive = early; negative = late.

	node.table <- data.frame(selected_geneList,pcdf,fc,sample.size)
	colnames(node.table) <- c('Gene','P_CDF',	'FC',	'Occurrence')
	write.table(node.table,"TEDGnode.txt",row.names = F,quote = F,sep = '\t')

	node <- data.frame(node.table)
	edge <- data.frame(edge.table)

	net <- graph_from_data_frame(d=edge[which(as.numeric(edge$weight) > 0),], vertices=node, directed=T)

	color.gradient <- function(x, colors=c("darkgreen","yellow","red"), colsteps=10000) {
	  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
	}
	V(net)$color <- color.gradient(V(net)$FC)

	V(net)$size <- V(net)$Occurrence*0.2+5
	E(net)$width <- as.numeric(E(net)$weight)/2

	E(net)$label <- NA

	E(net)$arrow.size <- 0.4
	E(net)$edge.color <- "gray90"

	graph_attr(net, "layout") <- layout_with_lgl
 	pdf(file="Figure7_TEDG.pdf", bg = 'white', width = 8, height = 6)
	plot(net)
 	dev.off()
	#figure_7<-rbind(ggplotGrob(net), size="last")
	#grid.draw(figure_7)
	#ggsave(file="Figure7_TEDG.pdf", plot=figure_7, bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)

	plot(net)

	returnList <- list("edge.table" = edge.table, "node.table" = node.table)
	return(returnList)
}




