d <- read.table("mds.mds",header=T)
cols <- rep("gray",nrow(d))
cols[grepl("[1-817]",d$FID)] <- "grey"
cols[grepl("[1000-100000000]",d$FID)] <- 'red'
#cols[grepl("^NA",d$FID)] <- 'green'
#cols[grepl("^Y",d$FID)] <- 'blue'

pchs = rep(11,nrow(d))
pchs[grepl("[1-817]",d$FID)] <- 21
pchs[grepl("[1000-100000000]",d$FID)] <- 15
#pchs[grepl("^NA",d$FID)] <- 17
#pchs[grepl("^Y",d$FID)] <- 19
