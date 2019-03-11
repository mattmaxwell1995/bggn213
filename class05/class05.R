#---
#title: "Class5 R Graphics"
#author: Matt Maxwell
#date: 1/25/19
#output: html_document 
#---



# Class 05 Graphics and plots with R
#This is some narrative text that I can style **bold** and **italc** and add links to [webpages]()




# Section 2A line plot

#save data file into easy to call variable
#read.table() and first call the zip file then the txt file 
#(i.e., tell R where to look)

weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

View(weight)

# ctrl + enter will get a readout of a function to visualize
#Barry's black bloxes are helpful on the lesson page.

plot(weight$Age, weight$Weight, typ ="b",
     pch =15, col=c("red"), cex =1.5, lwd=2, ylim=c(2,10),
     xlab="Age (months)", ylab="Weight(kg)",
     main = "Baby weight with age")



##Section 2B Bar plot
#same workflow as before, but need to use a sep = function I found in the read.table help page.

mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")

View(mouse)

## Need to utilize par() function to increase the barplot margin size
#par() has a ton of options for modifying your graph
par(mar=c(3.1, 11.1, 4.1, 2))

barplot(mouse$Count, names.arg = mouse$Feature,
        las = 1, horiz = TRUE, ylab = "",
        main = "Number of features in the mouse GRCme38 genome",
        xlim = c(0, 80000))




# Section 3A Providing color vectors


mf <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
View(mf)
# mfc <- read.delim("bimm143_05_rstats/male_female_counts.txt")

barplot(mf$Count, names.arg =mf$Sample, col=rainbow(nrow(mf)), las=2, ylab="Counts", xlab = "Gender", main = "Male & Female Rainblow Plot")

?hcl

## Section 3B
#read delim function saves us time when dealing with complex data
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
View(genes)

#how many genes
nrow(genes)
#how many genes are changing
#want to adjust my par, I think
par(mar=c(5, 5, 2, 1.5))
palette(c("blue","gray","red"))
plot(genes$Condition1, genes$Condition2, col=genes$State,
     xlab = "Expression Condition 1", ylab = "Expression Condition 2")
