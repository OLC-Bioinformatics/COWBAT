#!/usr/bin/env Rscript

# http://www.cbs.dtu.dk/courses/27626/Exercises/Alignment_exercise.php
# samtools view HG00418_A.bam | cut -f9 > OLC-795_insertsizes.txt

library("ggplot2")

# The arguments from the calling python script are passed in using the commandArgs function
# these arguments include the path to use as well as the name of the strain
cmdArgs <- commandArgs(trailingOnly=TRUE)

# The name will be used often enough that it needed its own variable
name <- cmdArgs[2]

# Create a string of the name and the necessary additions for the full name of the csv file
csvFile <- paste(cmdArgs[1], "/", name, "_insertsizes.csv", sep='')

# Read the tables of insert sizes
a = read.csv(csvFile)

# Remove negative numbers - as these data are the distance between paired reads, there will be the difference
# between read1 and read2, and then the difference between read2 and read1 - this will be negative, and is unncessary
a.v = a[a[,1]>0,1]

# Create a data frame with the same data as above - I found this was necessary for the plotting
b = data.frame(x = a[a[,1]>0,1])

# Only include insert sizes smaller than 2000. One issue I noticed was that sometimes reads are mapped in such a way that
# there are very large calculated differences between them. This removes those large inserts
c = data.frame(x = b[b[,1]<2000,1]) 

# Lets get rid of outliers and use the 5% and 95% intervals to calculate mean and standard deviation:
mn = quantile(a.v, seq(0,1,0.05))[2]
mx = quantile(a.v, seq(0,1,0.05))[20]

# Mean
calcMean <- mean(a.v[a.v >= mn & a.v <= mx])

# Standard deviation
calcSD <- sd(a.v[a.v >= mn & a.v <= mx])     

# Prepare strings to annotate the plot
annotateMean <- paste("Mean: " , round(calcMean, digits=1), " bp")
annotateSD <- paste("Standard Deviation" , round(calcSD, digits=1), " bp")
title <- paste("Distribution of library insert sizes for sequencing run of strain", name)

# Create a plot of data from c (0 < insertSizes < 2000)
p <-ggplot(c, aes(x=c[,1]))

# Add a histogram to the plot - binwidth makes bins of 20 units
p + geom_histogram(binwidth=20, fill="white", color="black") +
  # Annotate the plot all pretty-like
  annotate("text", x=750, y=Inf, vjust = 6, label = annotateMean) + 
  annotate("text", x=750, y=Inf, vjust = 7.5, label = annotateSD) +
  scale_x_continuous(name="Insert size (bp)") + 
  scale_y_continuous(name="Frequency") +
  ggtitle(title)

# String of the output pdf filename
filename <- paste(name, "_insert_sizes.pdf", sep='')

# saves the chart as name_insert_sizes.pdf, in landscape formate with letter paper dimensions
ggsave(filename, width=11, height=8.5)

# String of the output txt filename
textOutput <- paste(name, "_insert_sizes.txt", sep='')

#Creates a text file with name, mean and SD that will be parsed later to include a master file 
sink(textOutput) + cat(name) + cat("\t") + cat(calcMean) + cat("\t") + cat(calcSD)  + sink()
