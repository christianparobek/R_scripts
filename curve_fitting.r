# Graphing Derrick's K13 data
# February 22, 2014
# uses smooth.spline to smooth a curve plotted between data points
# took a long time to figure this one out

data <- read.table("kelch_curves.txt", header=TRUE)

par(mar=c(5,8,4,10))

plot(data$q34_rate ~ data$frequency, xlim=c(0,0.010), ylim=c(0,0.3), type="n", axes=FALSE, ylab="False Minor Alleles per 100 Bases", xlab="Min. Minor Allele Frequency", cex.lab=1.7)

axis(1, at=0:10/1000, cex.axis=1.5) # at=0:10/1000
axis(2, at=c(0,.1,.2,.3), labels=c(0,10,20,30), las=1, cex.axis=1.5)

ss <- smooth.spline(data$frequency, data$q33_rate, spar=0.01)
lines(predict(ss, (0:60)/10000), lwd=3, lty=3) # only sampling out to 0.6% here, because if we do it further (1.2%) it gets messy and confusing because the spline torques the curves away from zero, which is their actual value out this far. I believe it is in fact more honest (and looks better!) to draw it this way.

ss <- smooth.spline(data$frequency, data$q34_rate, spar=0.01)
lines(predict(ss, (0:60)/10000), lwd=3, lty=2) # ditto as above

ss <- smooth.spline(data$frequency, data$q35_rate, spar=0.01)
lines(predict(ss, (0:100)/10000), lwd=3, lty=1)

legend(0.005, 0.25, legend=c("q33","q34","q35"), lty=c(3,2,1), lwd=2, cex=1.7, box.col = "white")
