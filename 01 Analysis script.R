#title: "Adrenals complete"
#author: "Tristan Philippe"
#date: "October 15, 2019"
#Purpose: Data analysis script for CBG KD Adrenal and Liver projects

#Data loading and wrangling
#!!!!Make sure to change the setwd.
#Make sure the "Counts.csv" and "Design table.csv" are similar to previous analyses. The row names of the design table must match the column names of the counts.
##Read design table and counts

#setwd("C:/Documents/University/PhDNSc Project/Results/7 CBG RNA seq/Adrenal complete analysis/Sex_CBG_Age")
library(limma)
library(tidyverse)

des<-read.csv("Design table.csv", header = TRUE, row.names = 1)#Read in the design
Cat.search <- read.csv("Category search.csv")#Read in the reference file

#Read the files
sample <- list.files(pattern="*.counts.genes")
dat <- data.frame(
  lapply(setNames(sample, make.names(gsub("*.counts.genes", "", sample))), 
         read.table), row.names = 1)
#remove all gene names except for rownames
selectcol <- seq(1,ncol(dat),2)
dat <- dat[,c(selectcol)]
#Check number of columns
ncol(dat)
nrow(des)
#Checking row and col names of dat and des to make sure they match
colnames(dat)
rownames(des)
#Fix the colnames to match the Design table
colnames(dat) <- rownames(des)


##Reorder the design table and data by sex then genotype

des <- des[order(factor(des$Genotype, levels = c("WT", "KO"))),]
des <- des[order(des$Sex),]
des <- des[order(des$Age),]
data <- t(dat)
data <- data[match(rownames(des), rownames(data)),]
dat <- t(data)
dm <- model.matrix(data=des, ~ Sex*Genotype*Age)#Make dm


##Remove all-zero rows (>3000) and "no_feature" etc

gdat <- data.frame(rowMeans(dat[,1:4]),rowMeans(dat[,5:8]),rowMeans(dat[,9:12]),rowMeans(dat[,13:16]),rowMeans(dat[,17:20]),rowMeans(dat[,21:24]),rowMeans(dat[,25:28]),rowMeans(dat[,29:32])) #row means per group
gdat <- data.frame(gdat > 5)*1
dat<-t(subset(t(dat), select = rowSums(gdat) > 1)) #remove genes with expression lower than 5 in all groups
dat<-dat[!rownames(dat) %in% c("no_feature", "ambiguous"),] #remove other miscelaneous tags
hist(rowSums(log2(dat + 1)))


{r setup}
knitr::opts_knit$set(root.dir = 'C:/Documents/University/PhDNSc Project/Results/Hammond CBG seq colab/Adrenal complete analysis')


##Basic heatmap to check for sample correlation and outliers

library(lattice)
heatmap <-heatmap(cor(dat), Colv=NA, Rowv=NA, symm = T)
levelplot((cor(dat)), col.regions=heat.colors)

dev.copy(pdf, 'heat sample correlation.pdf')
dev.off()
dev.copy(jpeg, 'heat sample correlation.jpeg')
dev.off()

There are few outliers, #53 is a little concerning.

#PCA

#Scale
Scdat <- t(scale(t(dat)))
#Use SVD
svd.dat <- svd(Scdat)
#Plot
barplot(svd.dat$d^2/sum(svd.dat$d^2), main = "Variance explained by Principle components", xlab = "Principal Components", ylab = "Percent Variance explained")
dev.copy(pdf, 'PCA principle components.pdf')
dev.off()
dev.copy(jpeg, 'PCA principle components.jpeg')
dev.off()
#Extracting $v
pcs <- data.frame(svd.dat$v)
#Recolour based on Sex.Genotype.Age
des.col <- data.frame(paste(des$Sex, des$Genotyp, des$Age, sep = "."))
rownames(des.col) <- rownames(des)
colnames(des.col) <- "group"
des.col$col <- c(rep("salmon", 4), rep("orchid", 4), rep ("cadetblue1", 4), rep("chartreuse", 4), rep("deeppink", 4), rep("purple", 4), rep("cadetblue4", 4), rep("chartreuse4", 4))

#Top 5 PCA dimensions
plot(pcs[ , c("X1", "X2", "X3", "X4", "X5")], pch = 21, cex = 1, bg = des.col$col)
dev.copy(pdf, 'PCA top 5.pdf')
dev.off()
dev.copy(jpeg, 'PCA top 5.jpeg')
dev.off()


#Limma pathway Main comparisons (ie three way anova)

#Factors:
SGA=factor(paste(des$Sex, des$Genotype, des$Age, sep="."))
Maindm<- model.matrix(~0 + SGA)
rownames(Maindm) <- rownames(des)
Maincm<-makeContrasts(levels=Maindm,
                      MainSex = (SGAF.KO.30+SGAF.WT.30+SGAF.KO.60+SGAF.WT.60)-(SGAM.KO.30+SGAM.WT.30+SGAM.KO.60+SGAM.WT.60),
                      MainGenotype = (SGAF.KO.30+SGAM.KO.30+SGAF.KO.60+SGAM.KO.60)-(SGAF.WT.30+SGAM.WT.30+SGAF.WT.60+SGAM.WT.60),
                      MainAge = (SGAF.KO.30+SGAM.KO.30+SGAF.WT.30+SGAM.WT.30)-(SGAF.KO.60+SGAM.KO.60+SGAF.WT.60+SGAM.WT.60),
                      MainInter = (SGAM.KO.30-SGAM.WT.30)-(SGAF.KO.30-SGAF.WT.30)-(SGAM.KO.60-SGAM.WT.60)-(SGAF.KO.60-SGAF.WT.60)
)
vdat <- voom(dat, design = Maindm)
vdatweights <- data.frame(vdat)
write.csv(vdatweights, "vData.csv")

Mainf<-lmFit(vdat, Maindm, weights=vdat$weights)
Mainf<-contrasts.fit(Mainf, Maincm)
Mainf<-eBayes(Mainf)

colnames(Mainf$coefficients)

cutoff = 1

#Small function to extract coefficients (and save as csv) and display a p value distribution
#Fill in coef = whatever the coefficient is and filename = what you want it saved as.
limma_fun1 <- function(coef, filename){
  df.1a <- topTable(Mainf, coef = coef, p.value=cutoff, n=nrow(vdat))
  df.1a <- df.1a[order(row.names(df.1a)), ]
  write.csv(df.1a, file = paste("Spreadsheet ", filename, " uncorrected.csv", sep = ""))
  assign(c(paste("Main.", filename, sep = "")), df.1a, envir = .GlobalEnv) #Save in Global environment
  
  Hist.1b <- hist(topTable(Mainf, coef = coef, number = 100000)[,"P.Value"], main = filename, ylab = "Number of genes", xlab = "P value")
  dev.copy(pdf, file = paste("Pval frequency distribution ", filename, ".pdf", sep = ""))
  dev.off()
  dev.copy(jpeg, file = paste("Pval frequency distribution ", filename, ".jpeg", sep = ""))
  dev.off() 
}

limma_fun1(coef= "MainGenotype", filename = "Genotype")

limma_fun1(coef= "MainSex", filename = "Sex")

limma_fun1(coef= "MainAge", filename = "Age")

limma_fun1(coef= "MainInter", filename = "Interaction")


#Variance explained by Genotype, Sex, Age or unexplained by model

library(variancePartition)
form <- ~ (1|Genotype) + (1|Sex) + Age
form.G.S.A <- ~ (1|Genotype:Sex:Age) + (1|Genotype:Sex) + (1|Genotype:Age) + (1|Genotype) + (1|Sex:Age) + (1|Sex) + Age
form.G.S <- ~ (1|Genotype:Sex) + (1|Genotype) + (1|Sex) + Age
form.G.A <- ~ (1|Genotype:Age) + (1|Genotype) + (1|Sex) + Age
form.S.A <- ~ (1|Sex:Age) + (1|Genotype) + (1|Sex) + Age

Varfun <- function(form, filename){
  Variance <- fitExtractVarPartModel(vdat, form, des)
  Var <- data.frame(Variance)
  write.csv(rbind(Var, colMeans(Var)), file = paste("Percent Variance ", filename,".csv", sep = ""))
  colMeans(Var)
  
  #plotVarPart(Var)
  #dev.copy(pdf, 'Percent Variance Violin.pdf')
  #dev.off()
  #dev.copy(jpeg, 'Percent Variance Violin.jpeg')
  #dev.off()
  
  Var$Residuals <- NULL
  barplot(height = colMeans(Var), ylab = "Percent of Variance", las=2)
  dev.copy(pdf, 'Percent Variance Barplot.pdf')
  dev.off()
  dev.copy(jpeg, 'Percent Variance Barplot.jpeg')
  dev.off()}

Varfun(form=form, filename = "Main")
Varfun(form=form.G.S.A, filename = "Interaction GSA")
Varfun(form=form.G.S, filename = "Interaction GS")
Varfun(form=form.G.A, filename = "Interaction GA")
Varfun(form=form.S.A, filename = "Interaction SA")


#Individual constrasts (ie Posthocs)
##Limma Pathway

Postcm <- makeContrasts(levels=Maindm,
                        M30v60     = SGAM.WT.30-SGAM.WT.60,
                        F30v60     = SGAF.WT.30-SGAF.WT.60,
                        MWTvFWT.30 = SGAM.WT.30-SGAF.WT.30,
                        MWTvFWT.60 = SGAM.WT.60-SGAF.WT.60,
                        MKOvMWT.30 = SGAM.KO.30-SGAM.WT.30,
                        FKOvFWT.30 = SGAF.KO.30-SGAF.WT.30,
                        MKOvMWT.60 = SGAM.KO.60-SGAM.WT.60,
                        FKOvFWT.60 = SGAF.KO.60-SGAF.WT.60,
                        MKOvFKO.30 = SGAM.KO.30-SGAF.KO.30,
                        MKOvFKO.60 = SGAM.KO.60-SGAF.KO.60
)

Postf<-lmFit(vdat, Maindm, weights=vdat$weights)
Postf<-contrasts.fit(Postf, Postcm)
Postf<-eBayes(Postf)


##Set up to analyse the individual comparisons (spreadsheets, frequency distribution, heatmaps, matching DEGs to KEGG categories and saving the data)

library(circlize)
library(ComplexHeatmap)
library(tidyverse)

cutoff=1
FDRcutoff <- 0.05
SGA.des <- data.frame(subset(des, select = c("Sex", "Genotype", "Age")))


Ratnames <- c(rownames(des))
colnames(vdatweights) <- Ratnames

#Big function that produces and saves (1) extracts all genes without any correction and uses this to (a) Save a csv, (b) produce a frequency distribution of all the genes per comparison (2) FDR corrects dataframes that are used to (a) dataframes with data and KEGG categories, (b) a heatmap that organizes the genes based on a group of interest, (c) a heatmap that uses the clusterign algorithm to organize the genes.

#When using specify the coef (or contrast name) as = "", filename as = "" that should be used to save everything.

limma_fun2 <- function(coef, Sx, Gn, Ag, filename){
  #Regular top table to extract DEGs
  df.1a <- topTable(Postf, coef = coef, p.value=cutoff, n=nrow(vdat))
  df.1a <- df.1a[order(row.names(df.1a)), ] #Arrange alphabetically
  write.csv(df.1a, file = paste("Spreadsheet ", filename, " uncorrected.csv", sep = ""))
  df.1b <- df.1a %>% rownames_to_column("gene")
  assign(c(paste("Post.", coef, sep = "")), df.1b, envir = .GlobalEnv) #Save in Global environment
  
  
  #P value distribution
  Hist.1 <- hist(topTable(Postf, coef = coef, number = 100000)[,"P.Value"], main = filename, xlab = "P value")
  dev.copy(pdf, file = paste("Pval frequency distribution ", filename, ".pdf", sep = ""))
  dev.off()
  dev.copy(jpeg, file = paste("Pval frequency distribution ", filename, ".jpeg", sep = ""))
  dev.off() 
  
  df.2 <- subset(df.1a, adj.P.Val<FDRcutoff)
  df.2a <- subset(vdatweights, rownames(vdatweights) %in% rownames(df.2)) #Subset data
  df.2a <- df.2a[order(rownames(df.2a)),] #Arrange alphabetically
  df.2b <- data.frame(df.2, df.2a) #merge
  
  #Adding KEGG terms if applicable to each row
  colnames (df.2b) <- c(paste("logFC for ", coef, sep = ""), paste("AveExpr for ", coef, sep = ""), paste("t for ", coef, sep = ""), paste("P.value for ", coef, sep = ""), paste("adj.p.val for ", coef, sep = ""), paste("B for ", coef, sep = ""), Ratnames) #Rename cols
  write.csv(df.2b, file = paste("Spreadsheet ",  filename, " FDR corrected.csv", sep = ""))
  df.2c <- df.2b %>% rownames_to_column("gene")
  df.2c <- merge(df.2c, Cat.search, by = "gene", all.x = TRUE) #Add categories
  assign(c(paste("FDR.dat.", coef, sep = "")), df.2c, envir = .GlobalEnv) #Save in Global environment
  write.csv(df.2c, file = paste("Spreadsheet ",  filename, " FDR corrected with Kegg.csv", sep = ""))
  
  if(nrow(df.2a)<2){print("No DEGs")}else{
    
    
    #Heatmaps
    Sc <- t(scale(t(df.2a))) #scale data
    rowMeans <- data.frame(rowMeans(Sc[,c(colnames(Sc) %in% rownames(des[des$Sex %in% Sx & des$Genotype %in% Gn & des$Age %in% Ag,]))])) #Row means match group of interest
    colnames(rowMeans) <- c("Means")
    rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
    Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means  
    
    htg <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htg)
    
    htc <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htc)
    if(nrow(Sc)>50){
      pdf(file=paste("Heatmap 'relative' ", filename, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      pdf(file=paste("Heatmap 'relative' ", filename, " order FWT60.pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, " order FWT60.jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
    }
    else{
      pdf(file=paste("Heatmap relative ", filename, ".pdf", sep = "")
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, ".jpeg", sep = "")
      )
      draw(htc)
      dev.off() 
      
      pdf(file=paste("Heatmap relative ", filename, " order FWT60.pdf", sep = "")
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, " order FWT60.jpeg", sep = "")
      )
      draw(htg)
      dev.off() 
      
    }
    pdf(file=paste("Heatmap same-size ", filename, ".pdf", sep = "")
    )
    draw(htc)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, ".jpeg", sep = "")
    )
    draw(htc)
    dev.off()
    
    pdf(file=paste("Heatmap same-size ", filename, " order FWT60.pdf", sep = "")
    )
    draw(htg)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, " order FWT60.jpeg", sep = "")
    )
    draw(htg)
    dev.off()
  }}


##Extracting individual coefficients for the KO

limma_fun2(coef = "FKOvFWT.30", Sx = "F", Gn = "WT", Ag = "30", filename = "CBG Female at 30")
#Some DEGs appear to follow the MWT at 30 days, but all of them change to a different pattern at 60 days.

limma_fun2(coef = "FKOvFWT.60", Sx = "F", Gn = "WT", Ag = "60", filename = "CBG Female at 60")
#It is more obvious at 60 days that the FKO either did not feminize, mascculinized, or that it deviates from both FWT and MWT. These points will be important to pull out of the data. Albeit the distinction between lack of feminization or masculinization is a almost philosophical.
#CBG is expressed at low levels in 30 day old females, but higher in 60 day old so CBG KO has a small effect on females at 30 days (18 genes), but a far larger one at 60 days (4624 genes). Some of these genes may not have femininized at 60 days or may diverge from either MWT or FWT and so are unique to the KO.

limma_fun2(coef = "MKOvMWT.30", Sx = "F", Gn = "WT", Ag = "30", filename = "CBG Male at 30")
#Only one gene is like the FKO, the others are not like either FWT, FKO or MWT, but at 30 days 4/5 match the MWT. While interesting there are so few, that it is not worth dissecting this at this level, but it will be interesting to know what they do.

limma_fun2(coef = "MKOvMWT.60", Sx = "F", Gn = "WT", Ag = "60", filename = "CBG Male at 60")
#4/6 genes follow the FKO and are different than FWT and MWT, 2/6 are just different then everyone else. As above, we should check what they do.
#CBG is expressed at low levels in males and Knocking out CBG has little to no effect on gene expression in males (5/6 DEGS).



Maindm
des.post <- data.frame(c(rep("A", 4), rep("B", 4), rep("A",24)), c(rep("C", 20), rep("D", 4), rep("C",8)), c(rep("E", 12), rep("G", 4), rep("E",16)), c(rep("H", 28), rep("I", 4)))
colnames(des.post) <- c("F.KO.30", "F.KO.60", "M.KO.30", "M.KO.60")
des.post$Sex <- des$Sex
des.post$Age <- des$Age
rownames(des.post) <- colnames(dat)
form.post <- ~(1|F.KO.30) + (1|F.KO.60) + (1|M.KO.30) + (1|M.KO.60) + (1|Sex) + Age

Variance.post <- fitExtractVarPartModel(vdat, form.post, des.post)

Var.post <- data.frame(Variance.post)

write.csv(rbind(Var.post, colMeans(Var.post)), "Percent Variance per group.csv")

colMeans(Var.post)

plotVarPart(Var.post)
dev.copy(pdf, 'Percent Variance per group Violin.pdf')
dev.off()
dev.copy(jpeg, 'Percent Variance per group Violin.jpeg')
dev.off()

Var.post$Residuals <- NULL
Var.post[Var.post == "NaN"] <- 0

barplot(height = colMeans(Var.post), ylab = "Total percent of Variance", las=2)
dev.copy(pdf, 'Percent Variance per group Barplot.pdf')
dev.off()
dev.copy(jpeg, 'Percent Variance per group Barplot.jpeg')
dev.off()



#Sexual dimorphsim analysis

limma_fun2(coef = "MWTvFWT.30", Sx = "F", Gn = "WT", Ag = "30", filename = "Sexual dimorphism at 30")
#Clear sexual dimorphism but the FKO follows the FWT at 30 days and it appears many but not all of these DEGs are similar at 60 days.
#Only 1/4 remain sexually dimorphic by 60 days, others have different expression levels.

limma_fun2(coef = "MWTvFWT.60", Sx = "F", Gn = "WT", Ag = "60", filename = "Sexual dimorphism at 60")
#The FKO does not clearly match the FWT at 60 days.
#Half the genes that do not appear sexually dimorphic at 30 days changed in the FWT, but not the MWT nor FKO. The other half match the FWT and change in both the MWT and FKO. Hence half the genes are feminized towards the post-pubertal FWT and half are masculinized towards the post-pubertal MWT.

limma_fun2(coef = "MKOvFKO.30", Sx = "F", Gn = "KO", Ag = "30", filename = "Sexual dimorphism at 30 in KO")

limma_fun2(coef = "MKOvFKO.60", Sx = "F", Gn = "KO", Ag = "60", filename = "Sexual dimorphism at 60 in KO")



#Final heatmaps

library(circlize)
library(ComplexHeatmap)
library(tidyverse)
Heatfun <- function(Post, filename, Sx1, Gn1, Ag1, Sx2, Gn2, Ag2){
  Heatdat <- subset(Post, adj.P.Val<FDRcutoff) %>% rownames_to_column("ffc") %>% column_to_rownames("gene")
  Heatdat <- subset(vdatweights, rownames(vdatweights) %in% rownames(Heatdat)) #Subset data
  Sc <- t(scale(t(Heatdat))) #scale data
  Sc <- data.frame(subset(Sc, select = colnames(Sc) %in% rownames(des[des$Sex %in% Sx1 & des$Genotype %in% Gn1 & des$Age %in% Ag1,])), subset(Sc, select = colnames(Sc) %in% rownames(des[des$Sex %in% Sx2 & des$Genotype %in% Gn2 & des$Age %in% Ag2,])))
  rowMeans <- data.frame(rowMeans(Sc[,1:4])) #Row means match group of interest
  colnames(rowMeans) <- c("Means")
  rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
  Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means 
  htg <- Heatmap(Sc,
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 show_column_names = FALSE, show_row_dend = FALSE, show_row_names = FALSE, use_raster = TRUE,
                 col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
                 heatmap_legend_param  = list(at = c(-4, -2, 0, 2, 4), labels = c("-4", "-2", "0", "2", "4"), title ="", legend_height=unit(1.43, "inches"), grid_width=unit(0.09, "inches"), labels_gp = gpar(fontsize = 8)) #Legend specs
  )
  draw(htg)
  
  pdf(file=paste("Final Heatmap ", filename, ".pdf", sep = ""), width=unit(1.4, "inches"), height=unit(.3, "inches"))
  draw(htg)
  dev.off()
  
}

#set height to 1.6 and legend_height to 1.43
Heatfun(Post=Post.MWTvFWT.60, filename="Sex diff 60 days", Sx1="M", Gn1="WT", Ag1="60", Sx2="F", Gn2="WT", Ag2="60")
Heatfun(Post=Post.FKOvFWT.60, filename="Female CBG 60 days", Sx1="F", Gn1="WT", Ag1="60", Sx2="F", Gn2="KO", Ag2="60")

#set height to .6 no legend
Heatfun(Post=Post.MWTvFWT.30, filename="Sex diff 30 days", Sx1="M", Gn1="WT", Ag1="30", Sx2="F", Gn2="WT", Ag2="30")

#set height to .3
Heatfun(Post=Post.MKOvMWT.60, filename="Male CBG 60 days", Sx1="M", Gn1="WT", Ag1="60", Sx2="M", Gn2="KO", Ag2="60")

#Keep KO
Heatfun <- function(Post, filename, Sx1, Gn1, Ag1, Sx2, Gn2, Ag2, Sx3, Gn3, Ag3){
  Heatdat <- subset(Post, adj.P.Val<FDRcutoff) %>% rownames_to_column("ffc") %>% column_to_rownames("gene")
  Heatdat <- subset(vdatweights, rownames(vdatweights) %in% rownames(Heatdat)) #Subset data
  Sc <- t(scale(t(Heatdat))) #scale data
  Sc <- data.frame(subset(Sc, select = colnames(Sc) %in% rownames(des[des$Sex %in% Sx1 & des$Genotype %in% Gn1 & des$Age %in% Ag1,])), subset(Sc, select = colnames(Sc) %in% rownames(des[des$Sex %in% Sx2 & des$Genotype %in% Gn2 & des$Age %in% Ag2,])),subset(Sc, select = colnames(Sc) %in% rownames(des[des$Sex %in% Sx3 & des$Genotype %in% Gn3 & des$Age %in% Ag3,])))
  rowMeans <- data.frame(rowMeans(Sc[,9:12])) #Row means match group of interest
  colnames(rowMeans) <- c("Means")
  rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
  Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means 
  htg <- Heatmap(Sc,
                 cluster_rows = TRUE, cluster_columns = FALSE,
                 show_column_names = FALSE, show_row_dend = FALSE, show_row_names = FALSE, use_raster = TRUE,
                 col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
                 heatmap_legend_param  = list(at = c(-4, -2, 0, 2, 4), labels = c("-4", "-2", "0", "2", "4"), title ="", legend_height=unit(1.43, "inches"), grid_width=unit(0.09, "inches"), labels_gp = gpar(fontsize = 8)) #Legend specs
  )
  draw(htg)
  
  pdf(file=paste("Final Heatmap ", filename, ".pdf", sep = ""), width=unit(1.4, "inches"), height=unit(1.6, "inches"))
  draw(htg)
  dev.off()
  
}
Heatfun(Post=Post.MWTvFWT.60, filename="Sex diff and CBG 60 days CLUSTERED", Sx1="M", Gn1="WT", Ag1="60", Sx2="F", Gn2="WT", Ag2="60", Sx3="F", Gn3="KO", Ag3="60")













#Template matching (aka sorting) sexually dimorphic, FKO, and age effects into the various patterns observed in venn diagrams and heatmaps.
##Set up

library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(KEGGREST) 
library(gageData)
library(gage)
library(org.Rn.eg.db)

data(kegg.sets.rn) # download the rat KEGG set
data(sigmet.idx.rn)
kegg.sets.rn = kegg.sets.rn[sigmet.idx.rn] 

#cutoff value can be changes as required
TMcutoff <- 0.05

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.

pval <- data.frame(Post.MWTvFWT.30$adj.P.Val, Post.MWTvFWT.60$adj.P.Val, Post.FKOvFWT.30$adj.P.Val, Post.FKOvFWT.60$adj.P.Val, Post.MKOvMWT.30$adj.P.Val, Post.MKOvMWT.60$adj.P.Val)
rownames(pval) <- Post.MWTvFWT.30$gene
pval <- t(pval)
write.csv(pval, file = "pval to check.csv")


#Similar function as limma_fun2 but for the data resulting from template matching and with the addition of (e) gene ontology from the KEGG datatbase (f) gene ontology from the GO database
#Where: TM = result from sorting, Pa = Post... from individual comparison A, Pb = Post... from individual comaprison B, Gn = genotype, Pc = Post... from individual comaprison C OR FALSE is optional
TM_statsfun <- function(TM, Pa, Pb, Pc, Sx, Gn, Ag, filename){
  Namea <- deparse(substitute(Pa))
  Nameb <- deparse(substitute(Pb))
  Namec <- deparse(substitute(Pc))
  NameTM <- deparse(substitute(TM))
  df.2a <- subset(vdatweights, rownames(vdatweights) %in% rownames(TM)) #Subset data
  df.2b <- df.2a %>% rownames_to_column("gene")
  TM  <- data.frame(TM) %>% rownames_to_column("gene")
  
  
  #Adding KEGG terms if applicable to each row
  if(Pc == FALSE){
    df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }else{
    df.2c <-     df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., Pc, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          paste("logFC for ", Namec, sep = ""), paste("AveExpr for ", Namec, sep = ""), paste("t for ", Namec, sep = ""), paste("P.value for ", Namec, sep = ""), paste("adj.p.val for ", Namec, sep = ""), paste("B for ", Namec, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }
  assign(c(paste("dat.", NameTM, sep = "")), df.2c, envir = .GlobalEnv) #Save in Global environment
  write.csv(df.2c, file = paste("Spreadsheet ",  filename, " FDR corrected KEGG.csv", sep = ""))
  
  
  #Heatmap
  if(nrow(df.2a)<2){print("No DEGs")}else{
    #Heatmap  
    Sc <- t(scale(t(df.2a))) #scale data
    #Arrange by row means
    rowMeans <- data.frame(rowMeans(Sc[,c(colnames(Sc) %in% rownames(des[des$Sex %in% Sx & des$Genotype %in% Gn & des$Age %in% Ag,]))])) #Row means match group of interest
    colnames(rowMeans) <- c("Means")
    rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
    Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means  
    
    htg <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htg)
    htc <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = TRUE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htc)
    if(nrow(Sc)>50){
      pdf(file=paste("Heatmap 'relative' ", filename, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      pdf(file=paste("Heatmap 'relative' ", filename, " ", Sx, Gn, Ag, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
    }
    else{
      pdf(file=paste("Heatmap relative ", filename, ".pdf", sep = "")
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, ".jpeg", sep = "")
      )
      draw(htc)
      dev.off() 
      
      pdf(file=paste("Heatmap relative ", filename, " ", Sx, Gn, Ag, ".pdf", sep = "")
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = "")
      )
      draw(htg)
      dev.off() 
      
    }
    pdf(file=paste("Heatmap same-size ", filename, ".pdf", sep = "")
    )
    draw(htc)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, ".jpeg", sep = "")
    )
    draw(htc)
    dev.off()
    
    pdf(file=paste("Heatmap same-size ", filename, " ", Sx, Gn, Ag, ".pdf", sep = "")
    )
    draw(htg)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = "")
    )
    draw(htg)
    dev.off()
    
    
    # You can feed KEGG  logFC or p values so it can rank them.
    df.2d <- df.2c[!duplicated(df.2c$gene),]
    df.2d$abs <- abs(df.2d[,8])
    df.2d <- df.2d %>%
      dplyr::select(gene, abs) %>%
      arrange(desc(abs)) %>%
      column_to_rownames("gene") 
    
    
    #Extract ontological information from KEGG database
    df.2e <- df.2d
    df.2e$entrezid = mapIds(org.Rn.eg.db,
                            keys=row.names(df.2d), 
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
    df.2e <- df.2e[complete.cases(df.2e),]
    df.2e <- df.2e %>%
      rownames_to_column("gene") %>%
      dplyr::select(entrezid, abs) %>%
      column_to_rownames("entrezid")
    
    KEGG <- gage(df.2e, gsets=kegg.sets.rn, same.dir = TRUE)
    lapply(KEGG, head, n = 30)
    KEGG.Res <- rbind(data.frame(KEGG$greater), data.frame(KEGG$less))
    KEGG.Res <- KEGG.Res[complete.cases(KEGG.Res),]
    write.csv(KEGG.Res, file=paste("Kegg log ", filename, ".csv", sep = ""))
    
  }}






#cor.test for df needed alter
cor.testdf <- function(x, Template)
{k<-cor.test(x,Template)$p.value}

#Similar to above, but tests how well the expected pattern matches the data using cor.test.
TM_statsfunB <- function(TM, Pa, Pb, Pc, Sx, Gn, Ag, filename, tTm, split){
  Namea <- deparse(substitute(Pa))
  Nameb <- deparse(substitute(Pb))
  Namec <- deparse(substitute(Pc))
  NameTM <- deparse(substitute(TM))
  df.2a <- subset(vdatweights, rownames(vdatweights) %in% rownames(TM)) #Subset data
  df.2b <- df.2a %>% rownames_to_column("gene")
  TM  <- data.frame(TM) %>% rownames_to_column("gene")
  
  #Adding KEGG terms if applicable to each row
  if(Pc == FALSE){
    df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }else{
    df.2c <-     df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., Pc, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          paste("logFC for ", Namec, sep = ""), paste("AveExpr for ", Namec, sep = ""), paste("t for ", Namec, sep = ""), paste("P.value for ", Namec, sep = ""), paste("adj.p.val for ", Namec, sep = ""), paste("B for ", Namec, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }
  assign(c(paste("dat.", NameTM, sep = "")), df.2c, envir = .GlobalEnv) #Save in Global environment
  write.csv(df.2c, file = paste("Spreadsheet ",  filename, " FDR corrected KEGG.csv", sep = ""))
  
  
  
  if(nrow(df.2a)<2){print("No DEGs")}else{
    #Heatmap  
    Sc <- t(scale(t(df.2a))) #scale data
    #Arrange by row means
    rowMeans <- data.frame(rowMeans(Sc[,c(colnames(Sc) %in% rownames(des[des$Sex %in% Sx & des$Genotype %in% Gn & des$Age %in% Ag,]))])) #Row means match group of interest
    colnames(rowMeans) <- c("Means")
    rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
    Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means  
    
    htg <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htg)
    htc <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = TRUE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htc)
    if(nrow(Sc)>50){
      pdf(file=paste("Heatmap 'relative' ", filename, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      pdf(file=paste("Heatmap 'relative' ", filename, " ", Sx, Gn, Ag, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
    }
    else{
      pdf(file=paste("Heatmap relative ", filename, ".pdf", sep = "")
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, ".jpeg", sep = "")
      )
      draw(htc)
      dev.off() 
      
      pdf(file=paste("Heatmap relative ", filename, " ", Sx, Gn, Ag, ".pdf", sep = "")
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = "")
      )
      draw(htg)
      dev.off() 
      
    }
    pdf(file=paste("Heatmap same-size ", filename, ".pdf", sep = "")
    )
    draw(htc)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, ".jpeg", sep = "")
    )
    draw(htc)
    dev.off()
    
    pdf(file=paste("Heatmap same-size ", filename, " ", Sx, Gn, Ag, ".pdf", sep = "")
    )
    draw(htg)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = "")
    )
    draw(htg)
    dev.off()
    
    # You can feed KEGG  logFC or p values so it can rank them.
    df.2d <- df.2c[!duplicated(df.2c$gene),]
    df.2d$abs <- abs(df.2d[,8])
    df.2d <- df.2d %>%
      dplyr::select(gene, abs) %>%
      arrange(desc(abs)) %>%
      column_to_rownames("gene") 
    
    
    #Extract ontological information from KEGG database
    df.2e <- df.2d
    df.2e$entrezid = mapIds(org.Rn.eg.db,
                            keys=row.names(df.2d), 
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
    df.2e <- df.2e[complete.cases(df.2e),]
    df.2e <- df.2e %>%
      rownames_to_column("gene") %>%
      dplyr::select(entrezid, abs) %>%
      column_to_rownames("entrezid")
    
    KEGG <- gage(df.2e, gsets=kegg.sets.rn, same.dir = TRUE)
    lapply(KEGG, head, n = 30)
    KEGG.Res <- rbind(data.frame(KEGG$greater), data.frame(KEGG$less))
    KEGG.Res <- KEGG.Res[complete.cases(KEGG.Res),]
    write.csv(KEGG.Res, file=paste("Kegg log ", filename, ".csv", sep = ""))
    
    #Correlation
    if(split == TRUE){
      df.2f <- data.frame(df.2a[,c(colnames(df.2a) %in% rownames(des[des$Age %in% Ag,]))])
      Correlation <- apply(df.2f, 1, cor.testdf, tTm)#Run correlation
      Correlation <- data.frame(Correlation)
      Cor.Hist <- hist(Correlation$Correlation, main = filename, xlim = c(0,1), breaks = seq(0, 1, l = 20))
      dev.copy(pdf, paste("P Value distrubiton correlation", filename, ".pdf"))
      dev.off()
      dev.copy(jpeg, paste("P Value distrubiton correlation", filename, ".jpeg"))
      dev.off()
    }else{
      df.2f <- df.2a
      Correlation <- apply(df.2f, 1, cor.testdf, tTm)#Run correlation
      Correlation <- data.frame(Correlation)
      Cor.Hist <- hist(Correlation$Correlation, main = filename, xlim = c(0,1), breaks = seq(0, 1, l = 20))
      dev.copy(pdf, paste("P Value distrubiton correlation", filename, ".pdf"))
      dev.off()
      dev.copy(jpeg, paste("P Value distrubiton correlation", filename, ".jpeg"))
      dev.off() 
    }
    
  }}



##Sort FWTvFKO genes as sexually dimorphic or not
Sex.30, Sex.60, FKOvFWT.30, FKOvFWT.60, MKOvMWT.30, MKOvMWT.60

TM.FWTvFKOsex <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[4,]<TMcutoff
]))
TM_statsfun(TM.FWTvFKOsex, Pa = Post.MWTvFWT.60, Pb = Post.FKOvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Sexually dimorphic genes affected by CBG")

TM.FWTvFKOnosex <- data.frame(sort(pval[1, pval[2,]>TMcutoff & pval[4,]<TMcutoff
]))
TM_statsfun(TM.FWTvFKOnosex, Pa = Post.MWTvFWT.60, Pb = Post.FKOvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Genes affectd by FKO that are not sexually dimorphic")





TM.MWTvFWT60not30 <- data.frame(sort(pval[1, pval[1,]>TMcutoff & pval[2,]<TMcutoff
]))
TM_statsfun(TM.MWTvFWT60not30, Pa = Post.MWTvFWT.60, Pb = Post.MWTvFWT.30, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Sexually dymorphic genes at 60")

TM.MWTvFWT30not60 <- data.frame(sort(pval[1, pval[1,]<TMcutoff & pval[2,]>TMcutoff
]))
TM_statsfun(TM.MWTvFWT30not60, Pa = Post.MWTvFWT.30, Pb = Post.MWTvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "30", filename = "Sexually dymorphic genes at 30")

TM.MWTvFWT30and60 <- data.frame(sort(pval[1, pval[1,]<TMcutoff & pval[2,]<TMcutoff
]))
TM_statsfun(TM.MWTvFWT30and60, Pa = Post.MWTvFWT.30, Pb = Post.MWTvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Sexually dymorphic genes at both 30 and 60")

TM.FWTvFKOsex60 <- data.frame(sort(pval[1, pval[1,]>TMcutoff & pval[2,]<TMcutoff & pval[4,]<TMcutoff
]))
TM_statsfun(TM.FWTvFKOsex60, Pa = Post.MWTvFWT.60, Pb = Post.FKOvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Sexually dimorphic genes affected by CBG at 60 not 30")

TM.FKOsex60 <- data.frame(sort(pval[1, pval[1,]>TMcutoff & pval[2,]<TMcutoff & pval[4,]>TMcutoff
]))
TM_statsfun(TM.FKOsex60, Pa = Post.MWTvFWT.60, Pb = Post.FKOvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Sexually dimorphic genes unaffected by CBG at 60 not 30")

TM.FWTvFKOnosex60 <- data.frame(sort(pval[1, pval[1,]>TMcutoff & pval[2,]>TMcutoff & pval[4,]<TMcutoff
]))
TM_statsfun(TM.FWTvFKOnosex60, Pa = Post.MWTvFWT.60, Pb = Post.FKOvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Genes affectd by FKO that are not sexually dimorphic at 60")








Puberty*Sex Hyperdimensional analysis
{r setup}
knitr::opts_knit$set(root.dir = 'C:/Documents/University/PhDNSc Project/Results/Hammond CBG seq colab/Adrenal complete analysis/Hyperdimensional')


#Extracting individual coefficients for puberty and sexual dimorphisms

limma_fun2(coef = "F30v60", Sx = "F", Gn = "WT", Ag = "60", filename = "Puberty Female")
#Can see some sex differences in sexual maturation. Some genes in the FKO are similar to MWT

limma_fun2(coef = "M30v60", Sx = "F", Gn = "WT", Ag = "60", filename = "Puberty Male")
#As with the females, some sex differences are observed. Some genes in the FKO are similar to MWT
#Male and females unsurprisingly show changes in the expression of 1742 genes at 30 days and 4122 at 60 days, underlining an effect of puberty on sexual differentiation.

limma_fun2(coef = "MWTvFWT.30", Sx = "F", Gn = "WT", Ag = "30", filename = "Sexual dimorphism at 30")
#Clear sexual dimorphism but the FKO follows the FWT at 30 days and it appears many but not all of these DEGs are similar at 60 days.
#Only 1/4 remain sexually dimorphic by 60 days, others have different expression levels.

limma_fun2(coef = "MWTvFWT.60", Sx = "F", Gn = "WT", Ag = "60", filename = "Sexual dimorphism at 60")
#The FKO does not clearly match the FWT at 60 days.
#Half the genes that do not appear sexually dimorphic at 30 days changed in the FWT, but not the MWT nor FKO. The other half match the FWT and change in both the MWT and FKO. Hence half the genes are feminized towards the post-pubertal FWT and half are masculinized towards the post-pubertal MWT.


#3D PCA plot to see how Sex, Age, and genotype relate to each other in 3D space

#PCA function 2 dimensions
PCAfun <- function(X, Y){
  par(xpd = TRUE, mar = c(4,4,4,8))
  plot(pcs[, c(X, Y)], bg = des.col$col, pch = 21, cex = 1.5)
  legend("topright", inset = c(-.21,0), legend = as.character(unique(des.col$group)), pt.bg = unique(des.col$col), pch = 21)
  text(pcs[, c(X, Y)], row.names(des.col), cex = 0.6, pos = 2)
  dev.copy(pdf, file = paste("PCA ", X, " by ", Y, ".pdf", sep = ""))
  dev.off()
  dev.copy(jpeg, file = paste("PCA ", X, " by ", Y, ".jpeg", sep = ""))
  dev.off()
}

PCAfun("X3", "X2")
PCAfun("X4", "X3")
PCAfun("X3", "X4")
PCAfun("X2", "X4")

#PCA along 3 dimensions
library(plot3Drgl)
#"3" dimensional PCA but can choose angle
text3D(pcs$X2, pcs$X3, pcs$X4, labels = row.names(des.col), col = des.col$col, cex = 0.6, adj = 1, font = 1, phi = 30, theta = 50)
dev.copy(pdf, 'PCA X2 X3 X4 text.pdf')
dev.off()
dev.copy(jpeg, 'PCA X2 X3 X4 text.jpeg')
dev.off()
par(xpd = TRUE, mar = c(4,4,2,6))
scatter3D(pcs$X2, pcs$X3, pcs$X4, bg.col = des.col$col, pch = 21, cex = 1, colvar = FALSE, phi = 30, theta = 50)
legend("topright", inset = c(-0.2,0), legend = as.character(unique(des.col$group)), pt.bg = unique(des.col$col), pch = 21)
dev.copy(pdf, 'PCA X2 X3 X4.pdf')
dev.off()
dev.copy(jpeg, 'PCA X2 X3 X4.jpeg')
dev.off()

#3 dimensional interactive from scratch
rgl.fun1 <- function(x, y, z, colour){
  rgl.open()
  par3d(windowRect = c(20, 30, 800, 800))
  rgl.viewpoint(  zoom = .75 )
  rgl.bg(color = "white")
  rgl.spheres(x, y, z, color = colour, r = .02)
  rgl.lines(c(min(x), max(x)), c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(min(y),max(y)), c(0, 0), color = "red")
  rgl.lines(c(0, 0), c(0, 0), c(min(z),max(z)), color = "green")
}


rgl.fun1(x = pcs$X2, y = pcs$X3, z= pcs$X4, colour = des.col$col)
movie3d(spin3d(rpm = 6, axis = c(0, 1, 0)), duration = 10, type = "gif", dir = getwd(), movie = "PCA 3D")


#Venn diagrams to better understand what is happening in the FKO
##Female focused venn 30 days

C1 <- data.frame(merge(subset(Post.FKOvFWT.30, adj.P.Val<FDRcutoff), subset(Post.F30v60, adj.P.Val<FDRcutoff), by = "gene", all = T))
rownames(C1) <- C1$Row.names
C2 <- merge(C1, subset(Post.MWTvFWT.30, adj.P.Val<FDRcutoff), by = "gene", all = T)
C2[is.na(C2)] <- 0
nrow(C2) #Total number of genes in the set
C3 <- vennCounts(cbind(C2$adj.P.Val.x>0, C2$adj.P.Val.y, C2$`adj.P.Val`>0))
vennDiagram(C3, main = "At 30 days", names = c("FKOvFWT", "F30v60", "FWTvMWT"))
dev.copy(jpeg, 'Venn FKOvFWT-F30.60-FWTvMWT at 30 days.jpeg')
dev.off()
dev.copy(pdf, 'Venn FKOvFWT-F30.60-FWTvMWT at 30 days.pdf')
dev.off()

~1/3 of the DEGs in the FKO at 30 can be explained by innate sex differences, ~1/3 are different than both FWT and MWT, and ~1/4 are only different than FWT not MWT.
While there are 22X more sexually dymorphic DEGs; ~1/4 only in the FWT, 1/4 only in the FKO, and 1/2 in both.

##Female focused venn 60 days

C1 <- data.frame(merge(subset(Post.FKOvFWT.60, adj.P.Val<FDRcutoff), subset(Post.F30v60, adj.P.Val<FDRcutoff), by = "gene", all = T))
rownames(C1) <- C1$Row.names
C2 <- merge(C1, subset(Post.MWTvFWT.60, adj.P.Val<FDRcutoff), by = "gene", all = T)
C2[is.na(C2)] <- 0
nrow(C2) #Total number of genes in the set
C3 <- vennCounts(cbind(C2$adj.P.Val.x>0, C2$adj.P.Val.y, C2$`adj.P.Val`>0))
vennDiagram(C3, main = "At 60 days", names = c("FKOvFWT", "F30.60", "FWTvMWT"))
dev.copy(jpeg, 'Venn FKOvFWT-F30.60-FWTvMWT at 60 days.jpeg')
dev.off()
dev.copy(pdf, 'Venn FKOvFWT-F30.60-FWTvMWT at 60 days.pdf')
dev.off()

For the FKO there is a similar trend than at 30 days, but there hundreds of times more DEGs.
While the number of sexually dymorphic genes in the FWT has grown (2x) the number of genes that are sexually dymorphic in FKO has only increased by 10%.
Overall this suggests a lack of feminization at 60 days.

#Template matching (aka sorting) sexually dimorphic, FKO, and age effects into the various patterns observed in venn diagrams and heatmaps.
##Set up

library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(KEGGREST) 
library(gageData)
library(gage)
library(org.Rn.eg.db)
data(kegg.sets.rn) # dowload the rat KEGG set
data(sigmet.idx.rn)
kegg.sets.rn = kegg.sets.rn[sigmet.idx.rn] 

#cutoff value can be changes as required
TMcutoff <- 0.05

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.

pval <- data.frame(Post.MWTvFWT.30$adj.P.Val, Post.MWTvFWT.60$adj.P.Val, Post.F30v60$adj.P.Val, Post.M30v60$adj.P.Val, Post.FKOvFWT.30$adj.P.Val, Post.FKOvFWT.60$adj.P.Val, Post.MKOvMWT.30$adj.P.Val, Post.MKOvMWT.60$adj.P.Val)
rownames(pval) <- Post.MWTvFWT.30$gene
pval <- t(pval)
write.csv(pval, file = "pval to check.csv")


#Similar function as limma_fun2 but for the data resulting from template matching and with the addition of (e) gene ontology from the KEGG datatbase (f) gene ontology from the GO database
#Where: TM = result from sorting, Pa = Post... from individual comparison A, Pb = Post... from individual comaprison B, Gn = genotype, Pc = Post... from individual comaprison C OR FALSE is optional
TM_statsfun <- function(TM, Pa, Pb, Pc, Sx, Gn, Ag, filename){
  Namea <- deparse(substitute(Pa))
  Nameb <- deparse(substitute(Pb))
  Namec <- deparse(substitute(Pc))
  NameTM <- deparse(substitute(TM))
  df.2a <- subset(vdatweights, rownames(vdatweights) %in% rownames(TM)) #Subset data
  df.2b <- df.2a %>% rownames_to_column("gene")
  TM  <- data.frame(TM) %>% rownames_to_column("gene")
  
  
  #Adding KEGG terms if applicable to each row
  if(Pc == FALSE){
    df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }else{
    df.2c <-     df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., Pc, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          paste("logFC for ", Namec, sep = ""), paste("AveExpr for ", Namec, sep = ""), paste("t for ", Namec, sep = ""), paste("P.value for ", Namec, sep = ""), paste("adj.p.val for ", Namec, sep = ""), paste("B for ", Namec, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }
  assign(c(paste("dat.", NameTM, sep = "")), df.2c, envir = .GlobalEnv) #Save in Global environment
  write.csv(df.2c, file = paste("Spreadsheet ",  filename, " FDR corrected.csv", sep = ""))
  
  
  #Heatmap
  if(nrow(df.2a)<2){print("No DEGs")}else{
    
    Sc <- t(scale(t(df.2a))) #scale data
    htc <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = TRUE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htc)
    if(nrow(Sc)>50){
      pdf(file=paste("Heatmap 'relative' ", filename, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
    }
    else{
      pdf(file=paste("Heatmap relative ", filename, ".pdf", sep = "")
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, ".jpeg", sep = "")
      )
      draw(htc)
      dev.off() 
    }
    pdf(file=paste("Heatmap same-size ", filename, ".pdf", sep = "")
    )
    draw(htc)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, ".jpeg", sep = "")
    )
    draw(htc)
    dev.off()
    
    
    # You can feed KEGG  logFC or p values so it can rank them.
    df.2d <- df.2c[!duplicated(df.2c$gene),]
    df.2d$abs <- abs(df.2d[,8])
    df.2d <- df.2d %>%
      dplyr::select(gene, abs) %>%
      arrange(desc(abs)) %>%
      column_to_rownames("gene") 
    
    
    #Extract ontological information from KEGG database
    df.2e <- df.2d
    df.2e$entrezid = mapIds(org.Rn.eg.db,
                            keys=row.names(df.2d), 
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
    df.2e <- df.2e[complete.cases(df.2e),]
    df.2e <- df.2e %>%
      rownames_to_column("gene") %>%
      dplyr::select(entrezid, abs) %>%
      column_to_rownames("entrezid")
    
    KEGG <- gage(df.2e, gsets=kegg.sets.rn, same.dir = TRUE)
    lapply(KEGG, head, n = 30)
    KEGG.Res <- rbind(data.frame(KEGG$greater), data.frame(KEGG$less))
    KEGG.Res <- KEGG.Res[complete.cases(KEGG.Res),]
    write.csv(KEGG.Res, file=paste("Kegg log ", filename, ".csv", sep = ""))
    
  }}


#cor.test for df needed alter
cor.testdf <- function(x, Template)
{k<-cor.test(x,Template)$p.value}

#Similar to above, but produces extra heatmaps and tests how well the expected pattern matches the data using cor.test.
TM_statsfunB <- function(TM, Pa, Pb, Pc, Sx, Gn, Ag, filename, tTm, split){
  Namea <- deparse(substitute(Pa))
  Nameb <- deparse(substitute(Pb))
  Namec <- deparse(substitute(Pc))
  NameTM <- deparse(substitute(TM))
  df.2a <- subset(vdatweights, rownames(vdatweights) %in% rownames(TM)) #Subset data
  df.2b <- df.2a %>% rownames_to_column("gene")
  TM  <- data.frame(TM) %>% rownames_to_column("gene")
  
  
  #Adding KEGG terms if applicable to each row
  if(Pc == FALSE){
    df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }else{
    df.2c <-     df.2c <- left_join(TM, Pa, by = "gene", all.x = TRUE) %>% left_join(., Pb, by = "gene", all.x = TRUE) %>% left_join(., Pc, by = "gene", all.x = TRUE) %>% left_join(., df.2b, by = "gene", all.x = TRUE) %>% left_join(., Cat.search, by = "gene", all.x = TRUE) #merge individual contrasts, data, and catsearch
    df.2c <- df.2c[-c(2)]
    colnames (df.2c) <- c("gene",
                          paste("logFC for ", Namea, sep = ""), paste("AveExpr for ", Namea, sep = ""), paste("t for ", Namea, sep = ""), paste("P.value for ", Namea, sep = ""), paste("adj.p.val for ", Namea, sep = ""), paste("B for ", Namea, sep = ""),
                          paste("logFC for ", Nameb, sep = ""), paste("AveExpr for ", Nameb, sep = ""), paste("t for ", Nameb, sep = ""), paste("P.value for ", Nameb, sep = ""), paste("adj.p.val for ", Nameb, sep = ""), paste("B for ", Nameb, sep = ""),
                          paste("logFC for ", Namec, sep = ""), paste("AveExpr for ", Namec, sep = ""), paste("t for ", Namec, sep = ""), paste("P.value for ", Namec, sep = ""), paste("adj.p.val for ", Namec, sep = ""), paste("B for ", Namec, sep = ""),
                          Ratnames, "gene.1", "KeggPathway") #Rename cols
  }
  assign(c(paste("dat.", NameTM, sep = "")), df.2c, envir = .GlobalEnv) #Save in Global environment
  write.csv(df.2c, file = paste("Spreadsheet ",  filename, " FDR corrected.csv", sep = ""))
  
  
  
  if(nrow(df.2a)<2){print("No DEGs")}else{
    #Heatmap  
    Sc <- t(scale(t(df.2a))) #scale data
    #Arrange by row means
    rowMeans <- data.frame(rowMeans(Sc[,c(colnames(Sc) %in% rownames(des[des$Sex %in% Sx & des$Genotype %in% Gn & des$Age %in% Ag,]))])) #Row means match group of interest
    colnames(rowMeans) <- c("Means")
    rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
    Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means  
    
    htg <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htg)
    htc <- Heatmap(Sc, heatmap_legend_param = list(title = "Legend"), column_title = filename,
                   row_title = "Differentially Expressed Genes",
                   cluster_rows = TRUE, cluster_columns = FALSE,
                   show_column_names = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
                   col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                   top_annotation = HeatmapAnnotation(Sex = SGA.des$Sex, Genotype = SGA.des$Genotype, Age = SGA.des$Age, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"), Age = c("30" = "gray", "60" = "black")))
    ) #Heatmap specifics cluster rows
    draw(htc)
    if(nrow(Sc)>50){
      pdf(file=paste("Heatmap 'relative' ", filename, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htc)
      dev.off()
      
      pdf(file=paste("Heatmap 'relative' ", filename, " ", Sx, Gn, Ag, ".pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap 'relative' ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = ""), width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc))))
      )
      draw(htg)
      dev.off()
    }
    else{
      pdf(file=paste("Heatmap relative ", filename, ".pdf", sep = "")
      )
      draw(htc)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, ".jpeg", sep = "")
      )
      draw(htc)
      dev.off() 
      
      pdf(file=paste("Heatmap relative ", filename, " ", Sx, Gn, Ag, ".pdf", sep = "")
      )
      draw(htg)
      dev.off()
      
      jpeg(file=paste("Heatmap relative ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = "")
      )
      draw(htg)
      dev.off() 
      
    }
    pdf(file=paste("Heatmap same-size ", filename, ".pdf", sep = "")
    )
    draw(htc)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, ".jpeg", sep = "")
    )
    draw(htc)
    dev.off()
    
    pdf(file=paste("Heatmap same-size ", filename, " ", Sx, Gn, Ag, ".pdf", sep = "")
    )
    draw(htg)
    dev.off()
    
    jpeg(file=paste("Heatmap same-size ", filename, " ", Sx, Gn, Ag, ".jpeg", sep = "")
    )
    draw(htg)
    dev.off()
    
    # You can feed KEGG  logFC or p values so it can rank them.
    df.2d <- df.2c[!duplicated(df.2c$gene),]
    df.2d$abs <- abs(df.2d[,8])
    df.2d <- df.2d %>%
      dplyr::select(gene, abs) %>%
      arrange(desc(abs)) %>%
      column_to_rownames("gene") 
    
    
    #Extract ontological information from KEGG database
    df.2e <- df.2d
    df.2e$entrezid = mapIds(org.Rn.eg.db,
                            keys=row.names(df.2d), 
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
    df.2e <- df.2e[complete.cases(df.2e),]
    df.2e <- df.2e %>%
      rownames_to_column("gene") %>%
      dplyr::select(entrezid, abs) %>%
      column_to_rownames("entrezid")
    
    KEGG <- gage(df.2e, gsets=kegg.sets.rn, same.dir = TRUE)
    lapply(KEGG, head, n = 30)
    KEGG.Res <- rbind(data.frame(KEGG$greater), data.frame(KEGG$less))
    KEGG.Res <- KEGG.Res[complete.cases(KEGG.Res),]
    write.csv(KEGG.Res, file=paste("Kegg log ", filename, ".csv", sep = ""))
    
    #Correlation
    if(split == TRUE){
      df.2f <- data.frame(df.2a[,c(colnames(df.2a) %in% rownames(des[des$Age %in% Ag,]))])
      Correlation <- apply(df.2f, 1, cor.testdf, tTm)#Run correlation
      Correlation <- data.frame(Correlation)
      Cor.Hist <- hist(Correlation$Correlation, main = filename, xlim = c(0,1), breaks = seq(0, 1, l = 20))
      dev.copy(pdf, paste("P Value distrubiton correlation", filename, ".pdf"))
      dev.off()
      dev.copy(jpeg, paste("P Value distrubiton correlation", filename, ".jpeg"))
      dev.off()
    }else{
      df.2f <- df.2a
      Correlation <- apply(df.2f, 1, cor.testdf, tTm)#Run correlation
      Correlation <- data.frame(Correlation)
      Cor.Hist <- hist(Correlation$Correlation, main = filename, xlim = c(0,1), breaks = seq(0, 1, l = 20))
      dev.copy(pdf, paste("P Value distrubiton correlation", filename, ".pdf"))
      dev.off()
      dev.copy(jpeg, paste("P Value distrubiton correlation", filename, ".jpeg"))
      dev.off() 
    }
    
  }}




##Male and Male/Female CBG at 30 and 60 days
Sex.30, Sex.60, F30v60, M30v60, FKOvFWT.30, FKOvFWT.60, MKOvMWT.30, MKOvMWT.60

TM.MFKO.30 <- data.frame(sort(pval[1, pval[7,]<TMcutoff
]))
TM_statsfun(TM.MFKO.30, Pa = Post.MKOvMWT.30, Pb = Post.FKOvFWT.30, Pc = FALSE, filename = "CBG Male + Male and Female at 30")

TM.MFKO.60 <- data.frame(sort(pval[1, pval[8,]<TMcutoff
]))
TM_statsfun(TM.MFKO.60, Pa = Post.MKOvMWT.60, Pb = Post.FKOvFWT.60, Pc = FALSE, filename = "CBG Male + Male and Female at 60")


##Female CBG at 30 days
Sex.30, Sex.60, F30v60, M30v60, FKOvFWT.30, FKOvFWT.60, MKOvMWT.30, MKOvMWT.60

TM.SexnoCBG.30 <- data.frame(sort(pval[1, pval[1,]<TMcutoff & pval[5,]>TMcutoff
]))
tTm.SexnoCBG.30 <- rep(c(1,1,0,0), each = 4)
TM_statsfunB(TM = TM.SexnoCBG.30, Pa = Post.FKOvFWT.30, Pb = Post.MWTvFWT.30, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "30", filename = "Sexual dymporphism no CBG effect at 30 ", tTm = tTm.SexnoCBG.30, split = TRUE)
#Many genes in the 30 day old are still sexually dimorphic.

TM.FKOonly.30 <- data.frame(sort(pval[1, pval[5,]<TMcutoff & pval[7,]>TMcutoff
]))
TM_statsfun(TM = TM.FKOonly.30, Pa = Post.FKOvFWT.30, Pb = Post.MKOvMWT.30, Pc = FALSE, filename = "Only Female CBG at 30")
#Some at 30 but not many

TM.FKONotfeminized.30 <- data.frame(sort(pval[1, pval[1,]<TMcutoff & pval[5,]<TMcutoff & pval[7,]>TMcutoff
]))
tTM.FKONotfeminized.30 <- rep(c(0,1,1,1), each = 4)
TM_statsfunB(TM = TM.FKONotfeminized.30, Pa = Post.FKOvFWT.30, Pb = Post.MWTvFWT.30, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "30", filename = "FKO not Feminized at 30", tTm = tTM.FKONotfeminized.30, split = TRUE)
#Three at 30 match male

TM.FKOUnique.30 <- data.frame(sort(pval[1, pval[1,]>TMcutoff & pval[5,]<TMcutoff & pval[7,]>TMcutoff
]))
tTm.FKOUnique.30 <- rep(c(0,1,0,0), each = 4)
TM_statsfunB(TM = TM.FKOUnique.30, Pa = Post.FKOvFWT.30, Pb = Post.MWTvFWT.30, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "30", filename = "FKO Unique at 30", tTm = tTm.FKOUnique.30, split = TRUE)
#Some at 30, three look like a further extreme than the males, and the others diverge from both male and female.


##Puberty
Sex.30, Sex.60, F30v60, M30v60, FKOvFWT.30, FKOvFWT.60, MKOvMWT.30, MKOvMWT.60

TM.Puberty.MF <- data.frame(sort(pval[1, pval[2,]>TMcutoff & pval[3,]<TMcutoff & pval[4,]<TMcutoff
]))
TM_statsfun(TM.Puberty.MF, Pa = Post.F30v60, Pb = Post.M30v60, Pc = FALSE, filename = "Puberty Male and Female")

TM.Puberty.Fonly <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[3,]<TMcutoff
]))
TM_statsfun(TM.Puberty.Fonly, Pa = Post.F30v60, Pb = Post.M30v60, Pc = FALSE, filename = "Puberty Female only")

TM.Puberty.Monly <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[4,]<TMcutoff
]))
TM_statsfun(TM.Puberty.Monly, Pa = Post.M30v60, Pb = Post.F30v60, Pc = FALSE, filename = "Puberty Male only")

TM.PubertynoCBG <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[3,]<TMcutoff & pval[6,]>TMcutoff
]))
tTm.PubertynoCBG <- rep(c(0,0,0,0,1,1,0,0), each = 4)
TM_statsfunB(TM.PubertynoCBG, Pa = Post.F30v60, Pb = Post.FKOvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Puberty Female only NO CBG effects", tTm = tTm.PubertynoCBG, split = FALSE)
#Only effects of puberty in females no effects of CBG (FKO mostly = FWT)


##Female CBG at 60 days
Sex.30, Sex.60, F30v60, M30v60, FKOvFWT.30, FKOvFWT.60, MKOvMWT.30, MKOvMWT.60

TM.SexnoCBG.60 <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[6,]>TMcutoff
]))
tTm.SexnoCBG.60 <- rep(c(1,1,0,0), each = 4)
TM_statsfunB(TM = TM.SexnoCBG.60, Pa = Post.FKOvFWT.60, Pb = Post.MWTvFWT.60, Pc = FALSE, Sx = "F", Gn = "WT", Ag = "60", filename = "Sexual dymorphism no CBG effect at 60", tTm = tTm.SexnoCBG.60, split = TRUE)

TM.FKOonly.60 <- data.frame(sort(pval[1, pval[6,]<TMcutoff & pval[8,]>TMcutoff
]))
TM_statsfun(TM = TM.FKOonly.60, Pa = Post.FKOvFWT.60, Pb = Post.MKOvMWT.60, Pc = FALSE, filename = "Only Female CBG at 60")

TM.FKOonly.60not30 <- data.frame(sort(pval[1, pval[6,]<TMcutoff & pval[8,]>TMcutoff & pval[5,]>TMcutoff & pval[7,]>TMcutoff
]))
TM_statsfun(TM = TM.FKOonly.60, Pa = Post.FKOvFWT.60, Pb = Post.MKOvMWT.60, Pc = FALSE, filename = "Only Female CBG at 60 only")

TM.FKONotfeminized.60 <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[3,]<TMcutoff & pval[6,]<TMcutoff & pval[8,]>TMcutoff & pval[5,]>TMcutoff & pval[7,]>TMcutoff                       ]))
tTm.FKONotfeminized.60 <- rep(c(0,0,0,0,1,0,0,0), each = 4)
TM_statsfunB(TM.FKONotfeminized.60, Pa = Post.FKOvFWT.60, Pb = Post.F30v60, Pc = Post.MWTvFWT.60, Sx = "F", Gn = "WT", Ag = "60", filename = "FKO not Feminized at 60 only", tTm = tTm.FKONotfeminized.60, split = FALSE)
#These DEGs did not mature and remain similar to the male-like pattern.

TM.FKOMasculinized.60 <- data.frame(sort(pval[1, pval[2,]<TMcutoff & pval[3,]>TMcutoff & pval[4,]<TMcutoff & pval[6,]<TMcutoff & pval[8,]>TMcutoff & pval[5,]>TMcutoff & pval[7,]>TMcutoff     ]))
tTm.FKOMasculinized.60 <- rep(c(0,0,0,0,0,1,1,1), each = 4)
TM_statsfunB(TM.FKOMasculinized.60, Pa = Post.FKOvFWT.60, Pb = Post.M30v60, Pc = Post.MWTvFWT.60, Sx = "F", Gn = "WT", Ag = "60", filename = "FKO Masculinized at 60 only", tTm = tTm.FKOMasculinized.60, split = FALSE)
#It seems that these are masculinized

TM.FKNotmature.60 <- data.frame(sort(pval[1, pval[2,]>TMcutoff & pval[3,]<TMcutoff & pval[6,]<TMcutoff & pval[8,]>TMcutoff & pval[5,]>TMcutoff & pval[7,]>TMcutoff     ]))
tTm.FKONotmature.60 <- rep(c(0,0,0,0,1,0,1,1), each = 4)
TM_statsfunB(TM.FKNotmature.60, Pa = Post.FKOvFWT.60, Pb = Post.F30v60, Pc = Post.M30v60, Sx = "F", Gn = "WT", Ag = "60", filename = "FKO not Mature at 60 only", tTm = tTm.FKONotmature.60, split = FALSE)
#These DEGs failed to mature in the FKO post-puberty

TM.FKOUnique.60 <- data.frame(sort(pval[1, pval[3,]>TMcutoff & (pval[4,]>TMcutoff | pval[2,]>TMcutoff & pval[4,]<TMcutoff) & pval[6,]<TMcutoff & pval[8,]>TMcutoff & pval[5,]>TMcutoff & pval[7,]>TMcutoff                       ]))
tTm.FKOUnique.60 <- rep(c(0,0,0,0,0,1,0,0), each = 4)
TM_statsfunB(TM.FKOUnique.60, Pa = Post.FKOvFWT.60, Pb = Post.F30v60, Pc = Post.M30v60, Sx = "F", Gn = "WT", Ag = "60", filename = "FKO Unique at 60 only", tTm = tTm.FKOUnique.60, split = FALSE)
#Some are unique, some may classify above with more power.


#Citations

citation("limma")
citation("ggplot2")
citation("tidyverse")
citation("ComplexHeatmap")
citation("reshape2")
citation("KEGGREST")
citation("gageData")






#Not in use

#Box plots with points for individual genes of interest
##Loading packages and making the function

library(reshape2)
library(ggplot2)
#Setting up the design table to work with ggplots
tibdes <- des %>%
  rownames_to_column("sample_id")
tibdes$SGA <- SGA
#Setting up the data to work with ggplots
tibvdatw <- vdatweights %>%
  rownames_to_column("gene")

#Setting up the two factors

#function based on Stat540 seminar4 to reeorganize the values of a single gene to work with ggplots.
tibfun <- function(expressionMatrix) {
  expressionMatrix <- expressionMatrix %>%
    as.data.frame() %>%
    column_to_rownames("gene") %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>% 
    melt(id = "sample_id") %>% 
    as_tibble() %>% 
    select(sample_id,
           gene = variable, 
           expression = value) %>%
    left_join(tibdes, by = "sample_id") %>%
    ggplot(aes(x = SGA, y = expression)) +
    geom_boxplot() +
    geom_point(position = position_dodge2(width = .5)) +
    labs(x = "Genotype", y = "Gene expression (logcpm)", title = expressionMatrix) +
    theme(plot.title = element_text(hjust = 0.5))
  return(expressionMatrix)
}

#Make a box plot for a particular gene for males. Change gene =="gene of interest" to the gene of interest and what you save it as.
tibfun(filter(tibvdatw, gene == "Nr1d2")) 
dev.copy(pdf, 'Box+point plot Nr1d2.pdf')
dev.off()



#Kegg search for heatmaps
##Loading everything setting up

library(KEGGREST)
library(gageData)
#library(gage)
library(org.Rn.eg.db)
library(circlize)
library(tidyverse)
library(dplyr)
library(grid)
library(ComplexHeatmap)

#Specify the pathways to analyse
Pathendocrine = list("rno00140", "rno04925", "rno04927", "rno04080", "rno00100", "rno04915")
Pathmetabolism = list("rno00061", "rno00062", "rno00071")
Pathsignaling = list("rno04014", "rno04010", "rno04012", "rno04310", "rno04330", "rno04340", "rno04350", "rno04390", "rno04370", "rno04371", "rno04064", "rno04668", "rno04066", "rno04068", "rno04020", "rno04070", "rno04072", "rno04071", "rno04024", "rno04022", "rno04151", "rno04152", "rno04150")


#Function to make heatmaps from the genes that relate to a particular KEGG term. Path = A KEGG term or pathway of intrest, data = data only dataframe to search through
keggfun <- function(Path){
  Kegg.path <- keggGet(Path)
  Kegg.path.df <- data.frame(Kegg.path[[1]][["GENE"]]) %>% column_to_rownames("Kegg.path..1.....GENE...")
  Kegg.gene <- mapIds(org.Rn.eg.db,
                      keys=row.names(Kegg.path.df), 
                      column="SYMBOL",
                      keytype="ENTREZID",
                      multiVals="first")
  Kegg.gene <- Kegg.gene[!is.na(Kegg.gene)]
  Sub.dat <- subset(data, rownames(data) %in% Kegg.gene)
  Sc <- t(scale(t(Sub.dat))) #scale data
  if(nrow(Sc)==0) { print("No Gene")
  } else {
    heat <- Heatmap(Sc,
                    heatmap_legend_param = list(title = "Legend"),
                    row_title = "DEGs",
                    cluster_rows = TRUE, cluster_columns = FALSE, show_row_dend = FALSE,
                    show_column_names = TRUE, show_row_names = TRUE,
                    col = colorRamp2(c(-4, 0, 4), c("blue","black", "yellow")),
                    top_annotation = HeatmapAnnotation(SGA.des, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red")))
    )
    pdf(file=paste(Path, "FWT vs MWT.pdf", sep = ""), width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*39)*as.numeric(nrow(sub.df))))
    )
    draw(heat)
    dev.off()
    jpeg(file=paste(Path, "FWT vs MWT.jpeg", sep = ""), width=480, height=as.numeric(36+26+(19*as.numeric(nrow(sub.df))))
    )
    draw(heat)
    dev.off()
  }
}
lapply(Pathendocrine, DATA , keggfun)
lapply(Pathmetabolism, DATA, keggfun)
lapply(Pathsignaling, DATA, keggfun)










#limma contrasts collapsing groups.
Actually leads to weird results....

Colcm <- makeContrasts(levels=Maindm, FKOunique = SGF.KO-(SGM.WT+SGF.WT)/2, male-like = SGF.WT-(SGF.KO+SGM.WT)/2, Notmale-like = SGM.WT-(SGF.KO+SGF.WT)/2)

Colf<-lmFit(vdat, Maindm, weights=vdat$weights)
Colf<-contrasts.fit(Colf, Colcm)
Colf<-eBayes(Colf)

cutoff=1

#FKOunique
HitsFKOunique <- topTable(Colf, coef = "FKOunique", p.value=cutoff, n=nrow(vdat))
HitsFKOunique <- HitsFKOunique[order(row.names(HitsFKOunique)), ]
write.csv(HitsFKOunique, "HitsFKOunique.csv")
HistFKOunique <- hist(topTable(Colf, coef = "FKOunique", number = 100000)[,"P.Value"], main="Unique to FKO")
dev.copy(pdf, 'PvalFKOunique.pdf')
dev.off()
dev.copy(jpeg, 'PvalFKOunique.jpeg')
dev.off()

#male-like
Hitsmale-like <- topTable(Colf, coef = "male-like", p.value=cutoff, n=nrow(vdat))
Hitsmale-like <- Hitsmale-like[order(row.names(Hitsmale-like)), ]
write.csv(Hitsmale-like, "Hitsmale-like.csv")
Histmale-like <- hist(topTable(Colf, coef = "male-like", number = 100000)[,"P.Value"], main="male-like")
dev.copy(pdf, 'Pvalmale-like.pdf')
dev.off()
dev.copy(jpeg, 'Pvalmale-like.jpeg')
dev.off()

#Notmale-like
HitsNotmale-like <- topTable(Colf, coef = "Notmale-like", p.value=cutoff, n=nrow(vdat))
HitsNotmale-like <- HitsNotmale-like[order(row.names(HitsNotmale-like)), ]
write.csv(HitsNotmale-like, "HitsNotmale-like.csv")
HistNotmale-like <- hist(topTable(Colf, coef = "Notmale-like", number = 100000)[,"P.Value"], main="Not male-like")
dev.copy(pdf, 'PvalNotmale-like.pdf')
dev.off()
dev.copy(jpeg, 'PvalNotmale-like.jpeg')
dev.off()

library(circlize)
library(ComplexHeatmap)
#FDR cutoff
FDRcutoff <- 0.05

#Unique to FKO
FDR.FKOunique <- subset(HitsFKOunique, adj.P.Val<FDRcutoff)
#Select data based on FDR cutoff
dat.Con <- subset(vdatweights, row.names(vdatweights) %in% row.names(FDR.FKOunique))
#Scale the data (essentially take the average of each row and subtract every sample from it). This is essential when making a heatmap with more than one column.
Sc.FKOunique <- t(scale(t(dat.Con)))

#Getting the means of all repeat samples
rowMeans <- data.frame(rowMeans(Sc.FKOunique[,c(colnames(Sc.FKOunique) %in% rownames(SG.des[c(SG.des$Genotype %in% "KO", SG.des$Sex %in% "F"),]))]))
colnames(rowMeans) <- c("Means")
#Ordering rowMeans based on value
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE]
#Ordering data based on rowMeans
Sc.FKOunique <- Sc.FKOunique[order(match(rownames(Sc.FKOunique), rownames(rowMeans))),]

heat <- Heatmap(Sc.FKOunique, heatmap_legend_param = list(title = "Legend"),
                row_title = "Differentially Expressed Genes",
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = TRUE, show_row_names = FALSE,
                col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                top_annotation = HeatmapAnnotation(SG.des, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"))))

pdf("Heat FKOunique.pdf", width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*as.numeric(nrow(Sc.FKOunique))))
)
draw(heat)
dev.off()
jpeg("Heat FKOunique.jpeg", width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc.FKOunique))))
)
draw(heat)
dev.off()

#male-like
FDR.male-like <- subset(Hitsmale-like, adj.P.Val<FDRcutoff)
#Select data based on FDR cutoff
dat.Con <- subset(vdatweights, row.names(vdatweights) %in% row.names(FDR.male-like))
#Scale the data (essentially take the average of each row and subtract every sample from it). This is essential when making a heatmap with more than one column.
Sc.male-like <- t(scale(t(dat.Con)))

#Getting the means of all repeat samples
rowMeans <- data.frame(rowMeans(Sc.male-like[,c(colnames(Sc.male-like) %in% rownames(SG.des[c(SG.des$Genotype %in% "KO", SG.des$Sex %in% "F"),]))]))
colnames(rowMeans) <- c("Means")
#Ordering rowMeans based on value
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE]
#Ordering data based on rowMeans
Sc.male-like <- Sc.male-like[order(match(rownames(Sc.male-like), rownames(rowMeans))),]

heat <- Heatmap(Sc.male-like, heatmap_legend_param = list(title = "Legend"),
                row_title = "Differentially Expressed Genes",
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = TRUE, show_row_names = FALSE,
                col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                top_annotation = HeatmapAnnotation(SG.des, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"))))

pdf("Heat male-like.pdf", width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*as.numeric(nrow(Sc.male-like))))
)
draw(heat)
dev.off()
jpeg("Heat male-like.jpeg", width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc.male-like))))
)
draw(heat)
dev.off()

#Not male-like
FDR.Notmale-like <- subset(HitsNotmale-like, adj.P.Val<FDRcutoff)
#Select data based on FDR cutoff
dat.Con <- subset(vdatweights, row.names(vdatweights) %in% row.names(FDR.Notmale-like))
#Scale the data (essentially take the average of each row and subtract every sample from it). This is essential when making a heatmap with more than one column.
Sc.Notmale-like <- t(scale(t(dat.Con)))

#Getting the means of all repeat samples
rowMeans <- data.frame(rowMeans(Sc.Notmale-like[,c(colnames(Sc.Notmale-like) %in% rownames(SG.des[SG.des$Sex %in% "F",]))]))
colnames(rowMeans) <- c("Means")
#Ordering rowMeans based on value
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE]
#Ordering data based on rowMeans
Sc.Notmale-like <- Sc.Notmale-like[order(match(rownames(Sc.Notmale-like), rownames(rowMeans))),]

heat <- Heatmap(Sc.Notmale-like, heatmap_legend_param = list(title = "Legend"),
                row_title = "Differentially Expressed Genes",
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = TRUE, show_row_names = FALSE,
                col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
                top_annotation = HeatmapAnnotation(SG.des, col = list(Sex = c("F" = "deeppink", "M" = "Blue"), Genotype = c("WT" = "green", "KO" = "red"))))

pdf("Heat Notmale-like.pdf", width=7, height=as.numeric((7/958*70) + (7/958*45) + ((7/958*2)*as.numeric(nrow(Sc.Notmale-like))))
)
draw(heat)
dev.off()
jpeg("Heat Notmale-like.jpeg", width=480, height=as.numeric(36+26+(1*as.numeric(nrow(Sc.Notmale-like))))
)
draw(heat)
dev.off()



