


setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

mip <- data.frame(readRDS("Table/mipMat.rds"))
mro <- data.frame(readRDS("Table/mroMat.rds"))

exc <- c("X_0319")

mip2 <- mip[rownames(mip)!=exc, colnames(mip)!=exc, ]
mro2 <- mro[rownames(mro)!=exc, colnames(mro)!=exc, ]

saveRDS(mip2, "./Table/mip.rds")
saveRDS(mro2, "./Table/mro.rds")