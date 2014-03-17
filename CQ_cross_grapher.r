# grapher to generate plots of reference and alternate allele frequencies
# for all monkey samples combined
# Original version November 2013, final version February 5, 2014
# Christian Parobek
# Bootstrapping monkey replicates
# Wilcoxon test / t-test for difference in means



####################################################
#### CHECK WINDOW SIZE AND NUMBER OF REPLICATES ####
####################################################
win <- 20 # window size on just one side, so total window size is 2N+1
reps <- 1000 # number of reps for bootstrapping

# read in data
aotus_WR454 <- read.table("aotus_WR454.out", header=TRUE)
#aotus_86121 <- read.table("aotus_86121.out", header=TRUE)
#aotus_85823 <- read.table("aotus_85823.out", header=TRUE)
aotus_85986 <- read.table("aotus_85986.out", header=TRUE)
aotus_85922 <- read.table("aotus_85922.out", header=TRUE)
aotus_85850 <- read.table("aotus_85850.out", header=TRUE)
saimiri_4919 <- read.table("saimiri_4919.out", header=TRUE)
#saimiri_5076 <- read.table("saimiri_5076.out", header=TRUE)
saimiri_5081 <- read.table("saimiri_5081.out", header=TRUE)



####################################################
#### CALCULATE MEDIAN SLIDING WINDOW ####
####################################################
# declare some needed objects
preCQrefFreqs <- 0
preCQaltFreqs <- 0
postCQrefFreqs <- 0
postCQaltFreqs <- 0

# loop through all objects select the data objects # get returns an object from its name
for (i in ls(pattern="aotus*|saimiri*")) {
	
	# pre-CQ reference/sensitive frequencies
	s1_tot_read_ct <- get(i)$s1_ref + get(i)$s1_alt
	preCQrefFreqs <- preCQrefFreqs + get(i)$s1_ref / s1_tot_read_ct
	# post-CQ reference/sensitive frequencies
	s2_tot_read_ct <- get(i)$s2_ref + get(i)$s2_alt
	postCQrefFreqs <- postCQrefFreqs + get(i)$s2_ref / s2_tot_read_ct

}

# divide overall frequencies by length of monkeys
preCQrefFreqs <- preCQrefFreqs / length(ls(pattern="aotus*|saimiri*"))
postCQrefFreqs <- postCQrefFreqs / length(ls(pattern="aotus*|saimiri*"))

# Calculate sliding window and pick median

preCQmedian <- NULL
postCQmedian <- NULL
preCQmedian_sd <- NULL
postCQmedian_sd <- NULL

for (i in (1 + win):(length(preCQrefFreqs) - win)) {

preCQmedian <- append(preCQmedian, median(preCQrefFreqs[(i-win):(i+win)]))
postCQmedian <- append(postCQmedian, median(postCQrefFreqs[(i-win):(i+win)]))

preCQmedian_sd <- append(preCQmedian_sd, sd(preCQrefFreqs[(i-win):(i+win)]))
#/sqrt(length((i-win):(i+win)) # put in above line to make SD into SE
postCQmedian_sd <- append(postCQmedian_sd, sd(postCQrefFreqs[(i-win):(i+win)]))
#/sqrt(length((i-win):(i+win)) # put in above line to make SD into SE

}



####################################
#### BOOTSTRAPPING & STATISTICS ####
####################################
monkeys <- ls(pattern="aotus*|saimiri*")

# calculate frequency of reference allele for each monkey for preCQ and postCQ timepoints
preCQrefFreqs <- NULL
postCQrefFreqs <- NULL

for (i in monkeys) {

	preCQrefFreqs <- cbind(preCQrefFreqs, (get(i)$s1_ref) / (get(i)$s1_ref + get(i)$s1_alt))
	postCQrefFreqs <- cbind(postCQrefFreqs, (get(i)$s2_ref) / (get(i)$s2_ref + get(i)$s2_alt))
}

# declare variables that will store all replicates for preCQ and postCQ
preCQ <- NULL
postCQ <- NULL

# replicator
for (n in 1:reps) {

	# random number generator
	values <- sample(1:length(monkeys), length(monkeys), replace = TRUE)

	# declare variables that will hold mean values for each replicate
	preCQmeans <- NULL
	postCQmeans <- NULL

	# select the mean of "monkeys" x window
	for (i in (1 + win):(length(preCQrefFreqs[,1]) - win)) {

		preCQmeans <- append(preCQmeans, mean(preCQrefFreqs[(i-win):(i+win),values]))
		postCQmeans <- append(postCQmeans, mean(postCQrefFreqs[(i-win):(i+win),values]))

	}

	# layer mean calls into the aggregate variables, by row
	preCQ <- rbind(preCQ, preCQmeans)
	postCQ <- rbind(postCQ, postCQmeans)
}


# do the Wilcoxon Test
pValues <- NULL
for (i in 1:length(preCQmeans)) {
	pValues <- append(pValues, wilcox.test(preCQ[1:reps,i], postCQ[1:reps,i], alternative="greater", paired=FALSE, mu=0.20, conf.level=0.95)$p.value)
}
as.logical(pValues*length(pValues)<0.05)



##################
#### PLOTTING ####
##################
# setup plot parameters
pdf("chr1_final_figure.pdf")

options(scipen=20)
par(oma = c(8,2,8,0.5), las=1)
# oma = c(6,2,6,0.5)


# plot it
plot(preCQmedian ~ aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)],
xlab="",
ylab="Sensitive Allele Frequency",
pch=20,
cex=0.6,
col="black",
cex.lab=1,
type="n",
ylim=c(0,1),
axes=FALSE
)

axis(1, at=c(aotus_WR454$pos[(1+win)],aotus_WR454$pos[(length(aotus_WR454$pos)-win)]), labels=c("",""), lwd.ticks=0)
axis(side=1, at=c(100000,200000,300000,400000,500000,600000,700000,800000), labels=c("100kb","","300kb","","500kb","","700kb",""), lwd=0, lwd.ticks=1)
axis(side=2, at=c(0,0.25,0.5,0.75,1), labels=c("0","","0.5","","1.0"))

# plot pre-CQ points and Standard Deviation
segments(aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)], preCQmedian - preCQmedian_sd, aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)], preCQmedian + preCQmedian_sd, 
col="grey51",
lwd=0.4
)
points(preCQmedian ~ aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)],
pch=20,
cex=0.3,
col="black"
)

# plot post-CQ points and Standard Deviation
segments(aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)], postCQmedian - postCQmedian_sd, aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)], postCQmedian + postCQmedian_sd, 
col = rgb(1, 0, 0,0.3),
lwd=0.4
)
points(postCQmedian ~ aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)],
pch=20,
cex=0.3,
col="red"
)

# plot points for statistical significance
cond.astrx<-ifelse(as.logical(pValues*length(aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)])<0.05),"*",NA) # Logical if passes Bonferroni correction
text(aotus_WR454$pos[(1+win):(length(aotus_WR454$pos)-win)], 0.65, cond.astrx, cex=0.5) 
points(330000, 0, pch=8, cex=1.1) # add marker for location of pvcrt
mtext(expression(italic("pvcrt")), side=1, at=350000, cex=0.8)

# legend
legend(0,1, expression("Pre-CQ","Post-CQ"), col=c("black","red"), cex=1.1, pch=c(19,19),bg="white",horiz=FALSE, bty="n")

dev.off()

