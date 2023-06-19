# Import libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(magrittr)
source("main_functions.R")

# Read the command line argument as the file name
args <- commandArgs(trailingOnly = TRUE)

# FOR TESTING ONLY
args <- c('../data/GS-18S & 19S - PRISM BASE PEDIGREE UPDATE.csv',
         '../data/SorghumGSData_18_19.V3.csv',
         '../data/20201013_SB_AV_DArTag.csv')

logMessage(paste("Retrieved the following files: ", 
                  paste(args, collapse = ", ")))

# Check if a filename argument was provided
if (length(args) == 0) {
  stop("No file name provided.")
}

# Read the CSV files into data tables
data <- lapply(args, fread)

# Extract relevant columns from each data table
preY <- data[[1]][, c("Code", "Male", "Male pedigree", "Male Base Pedigree")]
preX <- data[[3]]
YADJ <- data[[2]]

# Convert column names to lowercase separated by underscores
setnames(preY, gsub(" ", ".", colnames(preY)))
setnames(preX, gsub(" ", ".", colnames(preX)))
setnames(YADJ, gsub(" ", ".", colnames(YADJ)))

# Preparation for X dataset
X <- t(preX)
X <- X[-c(1:19),]                       # TODO - verify this is standard for all
X <- ifelse(X=="1",2, ifelse(X=="2", 1, X))
rownames(X) <- gsub(rownames(X),pattern='.',replacement='-',fixed=TRUE)

# Preparation for Y dataset
YADJ$Year <- substr(YADJ$EntryList, 1, 2)
YADJ$Period <- paste(YADJ$Year, "S", sep='')

Y <- preY[!duplicated(preY$Code),]
Y <- merge(YADJ, Y, by='Code')
Y <- Y[substr(Y$Code, 1, 3)=="HAV",]    # TODO - verify this is standard for all

# Create an empty 2D matrix
index <- matrix(NA, nrow=dim(Y), ncol=2)
colnames(index) <- c('Group','Check')
index[,1] <- substr(Y$Male, 1,1)
index[,2] <- ifelse(substr(Y$Code,1,3)=="HAV", 0, 1)

# Bind index and Y, deduplicate dataset (rows withy same year, env, and genotype)
Y <- cbind(Y, index)[!duplicated(paste(Y$Year,'_',Y$Env,'_',Y$Male,sep='')),]

# Create two separate Y datasets
commonGidsMale <- intersect(rownames(X), Y$Male)
Y3 <- Y[Y$Male %in% commonGidsMale, ]
YN <- Y[!Y$Male %in% commonGidsMale,]

Y3.18 <- Y3[substr(Y3$EntryList,1,2)=="18",]
Y3.19 <- Y3[substr(Y3$EntryList,1,2)=="19",]

# Assigning pops to 2018
pops.18 <- unique(Y3.18$Male.Base.Pedigree)
tsp.18 <- grepl("VT", pops.18, fixed = TRUE)
Bs.18 <- pops.18[!tsp.18]
Rs.18 <- pops.18[tsp.18]
pops2.18 <- rbind(data.table(Male.Base.Pedigree = Bs.18, Pop = paste('B_', seq_along(Bs.18), sep = '')),
                  data.table(Male.Base.Pedigree = Rs.18, Pop = paste('R_', seq_along(Rs.18), sep = '')))

# Assigning pops to 2019
pops.19 <- names(table(Y3.19$Male.Base.Pedigree)[order(table(Y3.19$Male.Base.Pedigree), decreasing = TRUE)][1:8])
pops2.19 <- rbind(data.table(Male.Base.Pedigree = pops.19[grepl("BV", pops.19, fixed = TRUE)],
                             Pop = paste('B_', seq_along(pops.19[grepl("BV", pops.19, fixed = TRUE)]), sep = '')),
                  data.table(Male.Base.Pedigree = pops.19[!grepl("BV", pops.19, fixed = TRUE)],
                             Pop = paste('R_', seq_along(pops.19[!grepl("BV", pops.19, fixed = TRUE)]), sep = '')))

# Create the 'uids' column
Y3.18 <- cbind(Y3.18,1:dim(Y3.18)[1])
colnames(Y3.18)[dim(Y3.18)[2]] <- 'uids'

Y3.19 <- cbind(Y3.19,1:dim(Y3.19)[1])
colnames(Y3.19)[dim(Y3.19)[2]] <- 'uids'

# Merge pops with Y3 datasets
Y3.18 <- merge(Y3.18, pops2.18, by = 'Male.Base.Pedigree')[order(Y3.18$uids),]
Y3.19 <- merge(Y3.19[Y3.19$Male.Base.Pedigree %in% pops.19, ], pops2.19, by = 'Male.Base.Pedigree')[order(Y3.19$uids),]

# Means for 2018 and 2019
traits <- 6:8
TF.18 <- getMeansByGenotype(Y3.18, traits, 11)
TF.19 <- getMeansByGenotype(Y3.19, traits, "Male")

# Add a geno group column to TF.18 and TF.19 - maybe add to getMeansByGenotype()
GR.18 <- Y3.18$Male[Y3.18$EntryList %in% names(table(Y3.18$EntryList))[1]]
gps.18 <- ifelse(rownames(TF.18) %in% GR.18, "GR", "GB")
TF.18 <- data.frame(TF.18, Genotype.Group = gps.18)

GR.19 <- Y3.19$Male[Y3.19$EntryList %in% names(table(Y3.19$EntryList))[1]]
gps.19 <- ifelse(rownames(TF.19) %in% GR.19, "GR", "GB")
TF.19 <- data.frame(TF.19, Genotype.Group = gps.19)

# Split the datasets into genotype groups
Y3.18.B <- Y3.18[Y3.18$Group=="B",]
Y3.18.R <- Y3.18[Y3.18$Group=="R",]

Y3.19.B <- Y3.19[Y3.19$Group=="B",]
Y3.19.R <- Y3.19[Y3.19$Group=="R",]

# Visualizations of TF.18 and TF.19, saved in output/ folder

## Boxplot
filename <- "../output/boxplot.blues.across.envs.18.tiff"
# Create an empty list to store the boxplots made in loop
plots <- list()
for (t in 1:length(traits)) {
  plot_data <- data.frame(Group = TF.18[[length(traits) + 1]], Value = TF.18[[t]])
  
  plot <- ggplot(data = plot_data, aes(x = Group, y = Value)) +
    geom_boxplot() +
    labs(title = paste(Y3.18.B$Period[1], " ", colnames(TF.18)[1:length(traits)], "AP-AE"),
         ylab = "", xlab = "") +
    theme(plot.margin = margin(4, 4, 8, 1))
  
  plots[[t]] <- plot
}

# Combine the plots into a single image
combined_plot <- grid.arrange(grobs = plots, ncol = length(traits))

# Save the combined plot as a TIFF file
dir_path <- "../output/"
if (!dir.exists(dir_path)) { dir.create(dir_path) }
ggsave(filename, combined_plot, height = 6, width = 9, units = "in", dpi = 300)






