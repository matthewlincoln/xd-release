#!/usr/bin/env Rscript
### njc (adapted by mrl)
### 2016-10-25
### Renames subjects using our new system
### FIDs will be Fxxx-xxx where each xxx is a three-letter word
### IIDs will be Ixxx-xxx where each xxx is three digits

args <- commandArgs(trailingOnly = T)

# args[1] is the working directory where the .fam file is stored:
work_direc <- args[1]

# args[2] is the reference directory where the word list is stored:
ref_direc <- args[2]

set.seed(02241993)

# Read family data for each consortium:
# fam.file <- read.table(paste0(work_direc, "/all.pooled.qc.fam"), as.is = T)
fam.ced <- read.table(paste0(work_direc, "/ced/ced.pooled.qc.fam"), as.is = T)
fam.ibd <- read.table(paste0(work_direc, "/ibd/ibd.pooled.qc.fam"), as.is = T)
fam.ms <- read.table(paste0(work_direc, "/ms/ms.pooled.qc.fam"), as.is = T)
fam.sle <- read.table(paste0(work_direc, "/sle/sle.pooled.qc.fam"), as.is = T)
fam.t1d <- read.table(paste0(work_direc, "/t1d/t1d.pooled.qc.fam"), as.is = T)
fam.ra <- read.table(paste0(work_direc, "/ra/ra.pooled.qc.fam"), as.is = T)

# len.fam <- nrow(fam.file)
len.ced <- nrow(fam.ced)
len.ibd <- nrow(fam.ibd)
len.ms <- nrow(fam.ms)
len.sle <- nrow(fam.sle)
len.t1d <- nrow(fam.t1d)
len.ra <- nrow(fam.ra)

# Read three letter words:
words <- toupper(read.table(paste0(ref_direc, "/three.letter.words.new.lines.txt"), header = F, as.is = T)$V1)
# word.combs <- c(rep('', 1001*1001))
word.combs <- c(rep('', length(words)*length(words)))
# ns <- c(1:(len.ced + len.ibd + len.ms + len.sle.g + len.sle.o + len.t1d))
nums <- formatC(1:length(word.combs), width = 6, format = 'd', flag = '0')
num.combs <- paste0('I', substr(nums, 1, 3), '-', substr(nums, 4, 6))

position <- 1
for (i in words) {
  for (j in words) {
    w <- paste0('F', i, '-', j)
    word.combs[position] <- w
    position <- position + 1
  }
}

# Shuffle the IDs:
word.combs <- word.combs[sample(length(word.combs), size = length(word.combs), replace = F)]

# ind.fam <- c(1:len.fam)
ind.ced <- c(1:len.ced)
ind.ibd <- c(1:len.ibd + max(ind.ced))
ind.ms <- c(1:len.ms + max(ind.ibd))
ind.sle <- c(1:len.sle + max(ind.ms))
ind.t1d <- c(1:len.t1d + max(ind.sle))
ind.ra <- c(1:len.ra + max(ind.t1d))

# fam.recoding <- data.frame(fam.file$V1, fam.file$V2, word.combs[ind.fam], num.combs[ind.fam])
ced <- data.frame(fam.ced$V1, fam.ced$V2, word.combs[ind.ced], num.combs[ind.ced])
ibd <- data.frame(fam.ibd$V1, fam.ibd$V2, word.combs[ind.ibd], num.combs[ind.ibd])
ms <- data.frame(fam.ms$V1, fam.ms$V2, word.combs[ind.ms], num.combs[ind.ms])
sle <- data.frame(fam.sle$V1, fam.sle$V2, word.combs[ind.sle], num.combs[ind.sle])
t1d <- data.frame(fam.t1d$V1, fam.t1d$V2, word.combs[ind.t1d], num.combs[ind.t1d])
ra <- data.frame(fam.ra$V1, fam.ra$V2, word.combs[ind.ra], num.combs[ind.ra])

# Write renamed sample IDs:
# write.table(fam.recoding, paste0(work_direc, "/all.recoding.txt"),
# 	        row.names = F, col.names = F, quote = F)
write.table(ced, paste0(work_direc, "/ced/ced.recoding.txt"), row.names = F, col.names = F, quote = F)
write.table(ibd, paste0(work_direc, "/ibd/ibd.recoding.txt"), row.names = F, col.names = F, quote = F)
write.table(ms, paste0(work_direc, "/ms/ms.recoding.txt"), row.names = F, col.names = F, quote = F)
write.table(sle, paste0(work_direc, "/sle/sle.recoding.txt"), row.names = F, col.names = F, quote = F)
write.table(t1d, paste0(work_direc, "/t1d/t1d.recoding.txt"), row.names = F, col.names = F, quote = F)
write.table(ra, paste0(work_direc, "/ra/ra.recoding.txt"), row.names = F, col.names = F, quote = F)

fam.leftover <- word.combs[(max(ind.ra)+1):length(word.combs)]
ind.leftover <- num.combs[(max(ind.ra)+1):length(num.combs)]

# Write remaining family and individual names:
write.table(fam.leftover, paste0(work_direc, "/leftover.fam.names.txt"),
            row.names = F, col.names = F, quote = F)
write.table(ind.leftover, paste0(work_direc, "/leftover.ind.names.txt"),
            row.names = F, col.names = F, quote = F)