#!/usr/bin/env Rscript

### When running --check-sex in plink, the result is a series of F coefficients for the homozygosity of the X chromosome.
### This script uses clustering (Mclust) to assign sex based on the F coefficient.

library(ggplot2)
library(mclust)
library(plyr)

### This reads in the command line argument passed by immchip.ced.qc.sh, which is the file of F coefficients for the X chromosome.
### args[1] is the location of the .sexcheck file.
### args[2] is the location of log_direc/ced
### args[3] is the directory to write the lists for removal and updating to.
### args[4] is the consortium
args <- commandArgs(trailingOnly=T)

cons <- args[4]

### This is the .sexcheck file. It has the columns FID, IID, PEDSEX, SNPSEX, STATUS, and F.
sexprob <- read.table(args[1], header = T, as.is = T)

### We remove duplicates.
sexprob <- sexprob[!duplicated(sexprob$FID), ]

### Write stdout to this file.
sink(file = paste0(args[2], '/', cons, '.sex.check.analysis.txt'), type = 'output')

cat('In this analysis we compare the recorded sex of', cons, 'subjects to the plink inferred sex based on homozygosity of the X chromosome, and also to the sex we infer using the Mclust package to cluster subjects by X-chromosome homozygosity.\nSee also the plots in ', cons, '.sex.check.homozygosity.plots.pdf\n')

### How many are annotated as a problem?
cat('How many subjects are flagged as having problems for discordance between recorded and inferred sex?')
table(sexprob$STATUS)

cat('\nCompare PEDSEX (the recorded sex; y-axis) with SNPSEX (the inferred sex; x-axis).\n2 indicates female\n1 indicates male\n0 indicates unknown')
table(sexprob$PEDSEX, sexprob$SNPSEX)

### Make a histogram of the homozygosity rate distribution for all samples.
pdf(paste0(args[2], '/', cons, '.sex.check.plots.pdf'))

rep <- ggplot(sexprob, aes(x=F)) + geom_histogram(binwidth = 0.05) + facet_grid(~ PEDSEX)
rep

### Now do the clustering.
sexprob <- sexprob[!is.na(sexprob$F), ]
myclust <- Mclust(sexprob$F, modelNames = 'V', G = 2)

sexprob$cluster <- factor(myclust$classification)
### I make sure the mean homozygosity of men is higher than the mean homozygosity of women. If not, I flip the cluster labels.
switch.sex <- F
if (mean(sexprob$F[sexprob$cluster == 1]) < mean(sexprob$F[sexprob$cluster == 2])) {
    sexprob$cluster <- ifelse(sexprob$cluster == 1, 2, 1)
    switch.sex <- T
}


### Compare pedsex with model-assigned sex
p <- ggplot(sexprob, aes(as.factor(cluster), fill = as.factor(PEDSEX))) + geom_bar()
p <- p + xlab('Model-assigned sex')
p <- p + theme(legend.position = 'top')
p <- p + scale_fill_discrete(name = 'Reported sex')
p


### Compare snpsex with model-assigned sex
p <- ggplot(sexprob, aes(as.factor(cluster), fill = as.factor(SNPSEX))) + geom_bar()
p <- p + xlab('Model-assigned sex')
p <- p + theme(legend.position = 'top')
p <- p + scale_fill_discrete(name = 'Plink-assigned sex')
p

### Get the Mclust-determine probabilities of each sex.
probabilities <- as.data.frame(myclust$z)

if (switch.sex) names(probabilities) <- c('PFemale', 'PMale') else c('PMale', 'PFemale')

sexprob <- cbind(sexprob, probabilities)

cat('\nCompare PEDSEX (the recorded sex; y-axis) with clustering-model-inferred sex (x-axis).\n2 indicates female\n1 indicates male\n0 indicates unknown')
table(sexprob$PEDSEX, sexprob$cluster)


p <- ggplot(sexprob, aes(x=PFemale)) + geom_histogram(binwidth = 0.05) + xlab('P(female)') + ggtitle('Fitted model')
p
p <- ggplot(sexprob, aes(x=as.factor(cluster), y=PFemale)) + geom_boxplot()
p <- p + ylab('Model assigned probability of being female') + xlab('Sex: 1 = M, 2 = F') + ggtitle('Fitted model')
p

p <- ggplot(sexprob, aes(x=PMale)) + geom_histogram(binwidth = 0.05) + xlab('P(male)') + ggtitle('Fitted model')
p
p <- ggplot(sexprob, aes(x=as.factor(cluster), y=PMale)) + geom_boxplot()
p <- p + ylab('Model assigned probability of being male') + xlab('Sex: 1 = M, 2 = F') + ggtitle('Fitted model')
p


cat('\nSubjects to remove written to ', args[3], '/', cons, '.sex.discordance.to.remove.txt', sep = '')
cat('\nSubjects to update written to ', args[3], '/', cons, '.sex.updates.txt', sep = '')

sink()
dev.off()

### The truly worrisome subjects are the ones that have a sex recorded, but have the other sex according to our model.
### We don't trust these, so we're just going to throw them out.
### We're also going to throw out the (one) subject for which we could not calculate the measure of homozygosity (F).
write.table(sexprob[sexprob$PEDSEX != sexprob$cluster & sexprob$PEDSEX != 0, c('FID', 'IID')], paste0(args[3], '/', cons, '.sex.discordance.to.remove.txt'), row.names = F, col.names = F, quote = F)

write.table(sexprob[is.na(sexprob$F), c('FID', 'IID')], paste0(args[3], '/', cons, '.sex.discordance.to.remove.txt'), row.names = F, col.names = F, quote = F, append = T)

### The subjects that do not have a recorded sex will simply be assigned to their sex according to the model.
write.table(sexprob[sexprob$PEDSEX != sexprob$cluster & sexprob$PEDSEX == 0, c('FID', 'IID', 'cluster')], paste0(args[3], '/', cons, '.sex.updates.txt'), row.names = F, col.names = F, quote = F)
