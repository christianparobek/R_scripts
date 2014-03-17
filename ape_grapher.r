# To read in a NEWICK file, plot it, and color it appropriately
# NEWICK file can be generated in MEGA
# Originally started for PvMSP/CSP project, then adapted for PCP, and Jaymin's Guate paper
# Now used for Steve T's K13 project
# 16 March, 2014

# load the ape package
library(ape)

# read in the file
tree <- read.tree("11Mar2014_boot_consensus_tree.nwk")
#tree <- read.tree("12Mar2014_powermarker_output.nwk")

# initialize a vector of the correct size for colors according to name
colors <- character(length(tree$tip.label))
# color samples according to name/geographic location
colors[grep("Basse", tree$tip.label)] <- "black" # black for West Africa
colors[grep("Kita", tree$tip.label)] <- "black"
colors[grep("San", tree$tip.label)] <- "black"
colors[grep("Navrongo", tree$tip.label)] <- "black"
colors[grep("Ziniare", tree$tip.label)] <- "black"
colors[grep("Bas-Congo", tree$tip.label)] <- "grey" # grey for Central Africa
colors[grep("Kinshasa-2006", tree$tip.label)] <- "grey"
colors[grep("Kinshasa-2007", tree$tip.label)] <- "grey"
colors[grep("Bandundu", tree$tip.label)] <- "grey"
colors[grep("Tororo", tree$tip.label)] <- "white" # white for East Africa
colors[grep("Kisumu", tree$tip.label)] <- "white"
colors[grep("Siaya", tree$tip.label)] <- "white"
colors[grep("Fukayosi", tree$tip.label)] <- "white"
colors[grep("Machinga", tree$tip.label)] <- "white"

# rename sites so don't have to edit the NEWICK file
newNames <- c("Machinga, Malawi", "Kisumu, Kenya", "Kinshasa, DR Congo (2007)", "Kita, Mali", "San, Mali", "Tororo, Uganda", "Ziniare, Burkina Faso", "Basse, The Gambia", "Kinshasa, DR Congo (2006)", "Fukayosi, Tanzania", "Siaya, Kenya", "Bandundu, DR Congo", "Bas-Congo, DR Congo", "Navrongo, Ghana")
tree$tip.label <- newNames

par(oma=c(0,0,0,0), mar=c(1,0,1,0))
plot(tree, type="unrooted", direction="rightwards", show.tip.label=TRUE, cex=1, font=1, tip.color="black", label.offset=0.3, edge.width=0.60, lab4ut="axial")

tiplabels(pch=21, col="black", bg=colors, lwd=0.8, cex=2)

#
