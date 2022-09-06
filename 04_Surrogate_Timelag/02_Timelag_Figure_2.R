############################################################################
####
#### R script for Fujita (2019)
####
#### SpiecEasi
#### 2022.02.16 Kenta Suzuki
#### 2022.04.14 Hirokazu Toju 
#### R 4.1.2
#### 
############################################################################

setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

set.seed(123)

source('functions/functions.R')
source('functions/rEDM075_surrogate_codes.R')

# -- Create directory to save
dir <- make.dir('04_Surrogate_Timelag/02_output')

n.ex.rep <- 5 # number of tanks

library(doParallel)
library(foreach)
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)


############################################################################

d <- list()

for (j in 1:n.ex.rep) {
	CommDelay <- foreach(i=1:5) %dopar% read.table(sprintf('04_Surrogate_Timelag/01_output/Table/Surrogate_ASVs_CommDealy%s_All_th30_Tank%s_2.txt', i, j), header=T)
	
	NoDelay <- list(read.table(sprintf('04_Surrogate_Timelag/01_output/Table/Surrogate_ASVs_All_th30_Tank%s_2.txt', j), header=T))
	
	EnvDelay <- foreach(i=1:5) %dopar% read.table(sprintf('04_Surrogate_Timelag/01_output/Table/Surrogate_ASVs_EnvDealy%s_All_th30_Tank%s_2.txt', i, j), header=T)
	
	d[[j]] <- c(CommDelay, NoDelay, EnvDelay)
}

Delay <- -5:5
d2 <- list()

for (j in 1:n.ex.rep) {
	pr_X_0002 <- foreach(i=1:11, .combine='rbind') %dopar% d[[j]][[i]]$AS.result_Contr.pH[d[[j]][[i]]$ASV.ID=='X_0002'] # Best Cetobacterium
	pr_X_0027 <- foreach(i=1:11, .combine='rbind') %dopar% d[[j]][[i]]$AS.result_Contr.pH[d[[j]][[i]]$ASV.ID=='X_0027'] # Edwardsiella
	pr_X_0107 <- foreach(i=1:11, .combine='rbind') %dopar% d[[j]][[i]]$AS.result_Contr.pH[d[[j]][[i]]$ASV.ID=='X_0107'] # Aeromonas
	pr_X_0109 <- foreach(i=1:11, .combine='rbind') %dopar% d[[j]][[i]]$AS.result_Contr.pH[d[[j]][[i]]$ASV.ID=='X_0109'] # Worst Acinetobacter
	d2[[j]] <- cbind(Delay, pr_X_0002, pr_X_0027, pr_X_0107, pr_X_0109)
	colnames(d2[[j]]) <- c('Delay', 'pr_X_0002', 'pr_X_0027', 'pr_X_0107', 'pr_X_0109')
}

d3 <- rbind(d2[[1]], d2[[2]], d2[[3]], d2[[4]], d2[[5]])
Tank <- c(rep(1, times=11), rep(2, times=11), rep(3, times=11), rep(4, times=11), rep(5, times=11))

d4 <- data.frame(cbind(Tank, d3))
d4.min <- min(d4[,3:ncol(d4)])
d4.max <- max(d4[,3:ncol(d4)])


g1 <- ggplot(d4, aes(x=Delay, y= pr_X_0002, color=as.factor(Tank))) + geom_line() + labs(x= "Delay introduced to eel appetite score (days)", y= "Partial correlation with eel appetite score") + ggtitle("X_0002 (Cetobacterium)") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + geom_vline(xintercept=0, color='red', lty=2) + geom_hline(yintercept=0, color='red', lty=2) + #scale_y_continuous(limits = c(d4.min-0.05, d4.max+0.05)) + 
scale_x_continuous(breaks = Delay)

g2 <- ggplot(d4, aes(x=Delay, y= pr_X_0027, color=as.factor(Tank))) + geom_line() + labs(x= "Delay introduced to eel appetite score (days)", y= "Partial correlation with eel appetite score") + ggtitle("X_0027 (Edwardsiella)") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + geom_vline(xintercept=0, color='red', lty=2) + geom_hline(yintercept=0, color='red', lty=2) + 
#scale_y_continuous(limits = c(d4.min-0.05, d4.max+0.05))+ 
scale_x_continuous(breaks = Delay)


g3 <- ggplot(d4, aes(x=Delay, y= pr_X_0107, color=as.factor(Tank))) + geom_line() + labs(x= "Delay introduced to eel appetite score (days)", y= "Partial correlation with eel appetite score") + ggtitle("X_0107 (Aeromonas)") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + geom_vline(xintercept=0, color='red', lty=2) + geom_hline(yintercept=0, color='red', lty=2) + 
#scale_y_continuous(limits = c(d4.min-0.05, d4.max+0.05)) + 
scale_x_continuous(breaks = Delay)

g4 <- ggplot(d4, aes(x=Delay, y= pr_X_0109, color=as.factor(Tank))) + geom_line() + labs(x= "Delay introduced to eel appetite score (days)", y= "Partial correlation with eel appetite score") + ggtitle("X_0109 (Acinetobacter)") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + geom_vline(xintercept=0, color='red', lty=2) + geom_hline(yintercept=0, color='red', lty=2) + 
#scale_y_continuous(limits = c(d4.min-0.05, d4.max+0.05)) + 
scale_x_continuous(breaks = Delay)


g10 <- gridExtra::grid.arrange(g1, g2, g3, g4, nrow = 1)

ggsave(g10, filename=sprintf("%s/Timelag_Partial.Corr_ASVs.pdf", dir$figdir), h=3.5, w=16)


