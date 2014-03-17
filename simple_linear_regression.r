# for jaymin & martha's project
# simple linear regression, with two datasets on a single frame
# 17 March 2014

# read in data
data_IDA <- read.table("jaymin_data_IDA.txt", header=TRUE)
data_IR <- read.table("jaymin_data_IR.txt", header=TRUE)

# fit linear model
fit_IDA <- lm(invasions ~ X._IDA, data=data_IDA)
fit_IR <- lm(invasions ~ X._IR, data=data_IR)

# set up plot
plot(data_IR$invasions ~ data_IR$X._IR, type="n", ylim=c(0,1000000), xlab="Something", ylab="Invasions")

# plot individual points
points(data_IDA$invasions ~ data_IDA$X._IDA, col="red")
points(data_IR$invasions ~ data_IR$X._IR, col="blue")

# plot regression lines
abline(fit_IDA)
abline(fit_IR)
