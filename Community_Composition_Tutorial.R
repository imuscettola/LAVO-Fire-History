# Today we're going to be proceding with the samples from the same QIIME tutorial as last R class.
# We're going to bring them into R, to use a package called "phyloseq" to analyze and visualize them.
# https://joey711.github.io/phyloseq/

# Reminder: This is an R script. You can type directly into the R console, but it is best to execute commands from scripts.
# You can execute a line by typing [CMD+enter]. You can execute groups of lines by selecting them and doing the same.

# Comments are lines that start with a # sign. Even if they are executed, they are not run as code.
# Use comments throughout your scripts to remind yourself and future users of what you are doing at each step.

# The first thing you want to do when you open R is make sure you are in the right working directory.
# R will use whatever directory you set as its home base, from which it will look for any files you may have.
# Make sure this script and the final .biom file from our last week are in the same folder.
# Set the working directory under Misc in R, or under Session in RStudio, to be the folder where you have this script and the .biom folder.

# You can download the .biom file here: https://github.com/TheaWhitman/Soil_Micro_523/raw/master/feature-table-metaD-tax.biom


library("phyloseq")
library("ggplot2")
# This is the first line of code we'll run. This is used to "attach" the packages we will use in R.
# Phyloseq is a package of functions that helps us  analyze phylogenetic data.


?phyloseq
# In R, this is the command you run to get help or information on a function


ps = import_biom("feature-table-metaD-tax.biom", parseFunction=parse_taxonomy_greengenes)
# We are importing our .biom file and telling it the format the taxonomy is written in so it processes it correctly
# If you saved the file somewhere else, you'll need to direct it to the correct filepath
# Remember, it is best to have a folder for your code, and a separate folder for your data in your R project directory.
# (Don't worry about the warning message.)

ps
# This tells us a bit about our dataset
# We might want to check out the data a bit to see what it looks like.
# A key command for this is head()
# head() allows us to see only the top of something. That's great if it's, for example, a 380-row OTU table.

head(otu_table(ps))
# In the space below, try to look at the head of the sample data and our taxonomy table






# Maybe we're wondering how many sequences are in each sample. We can calculate and plot this:

d = sample_sums(ps)
# defining variable d as the column sums of our otu table.
d = data.frame(names(d),d)
# Creating a dataframe of our sample names
colnames(d)=c("Sample","Total")
# Naming the columns
head(d)
# displaying the top few values of d.


p = qplot(d$Total, geom="histogram", binwidth=60)
p
# This is our first use of ggplot, a figure-making package.
# You can make many different types of flexible, customized plots with ggplot.
# Plotting the Total values we calculated above, to see the total reads from each sample

#Some of our samples had very few sequences - some even had 0! (This is because we only used 1% of the total data, to #make analysis quicker.) For the purposes of this tutorial, we're going to want to get rid of the least abundant ones. #There's not a clear cutoff in the distribution we see above, so let's just take only samples with >200 sequences.
# phyloseq has a function to do this, called prune_samples

ps.cutoff = prune_samples(sample_sums(ps)>=200, ps)
ps.cutoff
# How many samples do we have now?
# How many samples are in the original ps object?
# In the space below, modify the code from above to see what happens if you change the cutoff from 200 to something else.







# Re-run the command with the cutoff set to something appropriate.


# Let's look at the taxonomic identity of the OTUs in our samples.
# We can use the plot_bar function from phyloseq, which actually draws on ggplot2.
plot_bar(ps.cutoff, fill="Phylum")



# Can you run the same bar graph command, but use a different taxonomic classification than Phylum? 
# What do you get?








# Okay, those plots are interesting, but they're the absolute abundances of each OTU
# We know that sequencing idiosyncracies, not real ecology, are likely driving differences.
# We might rather look at the relative abundance (fraction of total community)
# Phyloseq has a function to transform sample counts:
ps.norm = transform_sample_counts(ps.cutoff, function(x) x / sum(x) )
    
# And then we can make the same plot
plot_bar(ps.norm, fill="Phylum")



# We might also be interested in what information we have about the samples.
# You saw some of the categories when you ran head(sample_data(ps)) earlier.
# We can use the following command to look at just the column names of our sample data:
colnames(sample_data(ps.cutoff))

# Now we can group the bar charts we made above by different categories
# To do this, we add the facet_grid command to our plot_bar command:

plot_bar(ps.norm, fill="Phylum") + facet_grid(~Vegetation, scale="free") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
# Runs the plot_bar command, but grouped by Vegetation

# You can group by different variables.
# Try grouping by phylum by changing the facet_grid command




# Explore the data however you like!
# Use the ?plot_bar command to find out more options.




# Can you figure out how to give your graph a title?




# Look at the phyloseq stacked bar tutorial - can you figure out how to get rid of the black lines on the figure?
# (hint: look at the end of the tutorial)
# https://joey711.github.io/phyloseq/plot_bar-examples.html





# Now, try to import your group's .biom file. This may take some work.
# You may need to import the "mapping file" (which holds the sample data) separately, and add it to the OTU table and taxonomy table.






# Once you have it imported, you might want to save it as a phyloseq object so it's quick to import next time.
# You can do this as follows:
?saveRDS
saveRDS(YourPhyloseqObject, "YourPhyloseqObjectFilepath")
# This will save it in your working directory
# Remember you might want to save it in a separate data folder

# To import it next time, you can use
YourPhyloseqObject = readRDS("YourPhyloseqObjectFilepath")

# Note that the RDS format (R Data Structure) is not very "portable" - i.e.,
# You can't open the file using other non-R programs.
# Thus, it should just be used as a format for actively working with the data
# You can always go back to the original .biom file.



