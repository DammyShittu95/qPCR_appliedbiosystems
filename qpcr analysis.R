

#https://liz-is.github.io/qpcr-analysis-with-r/aio.html - for help
BiocManager::install("pcr", "readxl","devtools")
BiocManager::install("readexcel", "qpcR")
library(pcr)
library(readxl)
library(readexcel)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(tidyverse)
library(devtools)
devtools::install_github("ewallace/tidyqpcr",build_vignettes = TRUE) ## Vignettes require cowplot package
library(tidyqpcr)

#install, import and view excel file while selecting correct tab within workbook
#either change path file or better yet, add data file to extdata folder (ie external data storage for this particular package)
#locate and read raw data

#qpcr2_test15 <- read_excel(
#  system.file("extdata", "12522 AD quant -Rtudio.xls",
#              package = "tidyqpcr"), sheet = "Results", #opens doc in terms of the package being used?
#  skip=6, range = "A1:O104", col_names = TRUE # could also use 'skip=6, col_names = TRUE' which skips first 6 rows of unnecessary info and uses header as column names
#)


#ensure extdata folder contains file: /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tidyqpcr/extdata
setwd("/Library/Frameworks/R.framework/Versions/4.2/Resources/library/tidyqpcr/extdata")
qpcr_test15 <- read_excel(
  system.file("extdata", "12522 AD quant -Rtudio.xls",
              package = "tidyqpcr"), sheet = "Results", #opens doc in terms of the package being used?
  skip=6, col_names = TRUE #skips first 6 rows of unnecessary info and uses header as column names; #argument: 'range= "A1:O104"' is supposed to exclude specific rows and columns in favour for just the data
)
View(qpcr_test15)

                #if we didn't include primer targets; we would create a new dataframe to explain our primers and targets and combine with the qPCR data
                #primer_key <- data.frame(row = c("A", "B", "C", "D", "E", "F", "G", "H"),
                       #primers = c(rep("7280;c84", 2), rep("cas9-267", 2), rep("dsx_T1", 2), rep("dsx_T1", 2)))


#make wells read properly and create plate layout
#separates well columns and row coordinates to make data look tidier
library(rlang)
library(qpcR)
library(tidyr)
tidy_qpcr <- separate(qpcr_test15, Well, into = c("well_row", "well_col"), 
                             sep = 1) #separates by 1 position
tidy_qpcr
View(tidy_qpcr)


#add a variable for the 'Task' column and name 'prep_type' then change 'Target Name' column to 'target_id' etc. to be compatible with package and view
prep_type=tidy_qpcr$Task

names(tidy_qpcr)[names(tidy_qpcr) == "Target Name"] <- "target_id" #alters column name from 'Target Name' to 'target_id'
  names(tidy_qpcr)[names(tidy_qpcr) == "Sample Name"] <- "sample_id" #alters column name from 'Sample Name' to 'sample_id'
  names(tidy_qpcr)[names(tidy_qpcr) == "Cт"] <- "Ct"
  names(tidy_qpcr)[names(tidy_qpcr) == "Cт SD"] <- "Ct_sd"
##note this works: rename(new_name = old_name) renames individual variables; rename_with() = renames columns
View(tidy_qpcr)
  

#------
#create_blank_plate_96well() #or create_blank_plate(well_row = LETTERS[1:8], well_col=1:12)
#works despite error message
#lower rows are removed when making plate
#ensure biological replicate group is added to cvs/ qpcr set up otherwise bar plot colour doesnt distinguish between replicates - can use separate function if tidy; see hidden code below
plate_qpcr <- 
  create_blank_plate_96well() %>%
  inner_join(tidy_qpcr) %>% #inner_join joins data frames together
    mutate(prep_type = Task) #adds a prep_type column (i.e. Task column from original data) for info on cDNA prep; can be +RT; -RT; NTC; 'tidy_qpcr$prep_type <- prep_type' also works
display_plate_qpcr(plate_qpcr) #displays plate layout
View(plate_qpcr)




#----
##we can filter unnecessary sample types and add column to show biological replicate if needed
#separate(plate_qpcr, plate_qpcr$sample_id, into = c("sample_id", "bio replicate"), sep =" ", remove = FALSE)
#filter(plate_norm, sample_id != "sybr + h20 only" & Sample.Name != "")  #--> != excludes data type, & adds more variable and conditions

#----
#now plot basic Ct info
#first make own data compatible to the function by replacing 'Undetermined' string to a Ct of 40; 
#make column numeric as factors are stored as character “labels” to underlying numeric levels; converting directly to numeric returns the levels, not actual values --> as.numeric(as.character(my_factor))
plate_qpcr$Ct = str_replace_all(plate_qpcr$Ct, "Undetermined", "40")
plate_qpcr$Ct <- as.numeric(plate_qpcr$Ct)
plate_qpcr$Ct_sd <- as.numeric(plate_qpcr$Ct_sd)
.data <- plate_qpcr


ggplot(data = subset(plate_qpcr, target_id != "h20 only" & Task != "STANDARD"), aes( #filters rows to show everything bar water controls etc
        x = sample_id, y = Ct, fill = target_id, colour = target_id, shape = sample_id, label = sample_id, 
        ymin=Ct, ymax=Ct + Ct_sd)) +  #removing ymin=Ct-Ct_sd stops error bar covering text
#add all aes() features in first line
#changed filter to subset(); #uses normalised mean delta delta values and filter water control
 #geom_point(position = position_dodge(width = 1)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 1), colour = "black") + ##colour is for the bar's outline
    #'stat = "identity"' in geom_bar skips aggregating rows and we provide y-values; geom_col doesn't aggregate by default;
    #'position = position_dodge' controls width between bars; shape separates overlapping bars/ values and gives each sample an identity
      geom_errorbar(position = position_dodge(width = 1), width=.2, colour = "black") + #adds error bars; remove '-sd' to remove lower error bar
      geom_text(position = position_dodge(width = 1), hjust = 1.5, colour = "black", size = 1.5, angle=90) +
      labs(y="Ct", title = "All reps, Unnormalised") + #labs() modify axis, legend, and plot labels
    facet_wrap(~target_id, nrow = 1) + #creates panels - can be used for features such as time points
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



#create normalised Ct data

#first, rebuild tidyqpcr's 'calculate_deltacq_bysampleid' function to calculate_deltact_bysampleid
# Fn2 OR function name w/out '()' => opens trace function to open a package/ function's source code"
#df (ie plate_qpcr) contains: sample_id;  value_name - default Ct ; target_id - target genes --> in this case the ref_target is tRNA - 'rg1262/63 hkg control'
#crucial --> sample_id should be the same for technical replicates, but differ for different biological and experimental replicates. See tidyqpcr vignettes for examples --> in my case biological replicates differ by WT1 & WT2


#define hkg ie ref. target gene before rebuilding function
ref_target_ids <- plate_qpcr[plate_qpcr$target_id %in% "rg1262/63 hkg control",] #this subsets the column 'target_id' according to the hkg/ target gene

calculate_deltact_bysampleid <- function (plate_qpcr, ref_target_ids, norm_function = mean) {
  plate_qpcr %>% dplyr::group_by(.data$sample_id) %>% dplyr::do(calculate_normvalue(.data,
              ref_ids = ref_target_ids, value_name = "Ct", id_name = "target_id",   #groups by sample an primer/target type, then normalises ref. target_ids (ie target genes) Ct values; this will be used to calculate delta Ct for each sample_id
              norm_function = norm_function)) %>% dplyr::rename(ref_Ct = .data$value_to_norm_by) %>%  #new_name = old_name renames individual variables or columns; in this case ref_cq for the targets
    dplyr::ungroup() %>% dplyr::mutate(delta_Ct = .data$Ct - .data$ref_Ct,   #calculates normalised ΔCT in each sample: ΔCt = Ct(sample1) - ref target/HKG's Ct(sample1)
               rel_abund = 2^-.data$delta_Ct) #these lines ungroup columns back to original state, then creates or 'mutate()' an additional column named 'rel_abund' containing the 2^-.data$delta_Ct value ie expression change-fold
}

#then create normalised (delta Ct) dataframe
plate_norm <- calculate_deltact_bysampleid(plate_qpcr,ref_target_ids = "rg1262/63 hkg control", norm_function = mean) #calculate the mean (own code changed it from median) Ct of reference gene separately for each sample; then subtract it from all other sample's Ct
View(plate_norm) #note: samples like -RT are seen as NaN ie impossible values


#group technical replicates and calculate mean delta Ct and mean rel. abundance
plate_norm_by_mean <- plate_norm %>%
  group_by(sample_id, target_id) %>%
  summarize(
    delta_Ct  = mean(delta_Ct, na.rm = TRUE),
    rel_abund = mean(rel_abund, na.rm = TRUE))
View(plate_norm_by_mean)

#now plot the ΔCt & rel. abundance (2^-ΔCt)
plot_norm_by_mean <- ggplot(data = subset(plate_norm_by_mean, target_id != "h20 only"), aes(x = sample_id, y = delta_Ct, fill = sample_id, colour = target_id, label = sample_id)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 1), colour = "black") +
  geom_text(position = position_dodge(width = 1), hjust = -0.2, colour = "black", angle=90) +
  labs(y = "delta Ct", title = "mean of replicates, Normalised") +
  facet_wrap(~target_id,nrow=1) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
plot_norm_by_mean

plot_rel_abun <- ggplot(data = subset(plate_norm_by_mean, target_id != "h20 only"), aes(x = sample_id, y = rel_abund, fill = sample_id, colour = target_id, label = sample_id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), colour = "black") +
  geom_text(position = position_dodge(width = 1), hjust = -0.2, colour = "black", angle=90) +
  labs(y = "2^-ΔCt, abundance relative \n to housekeeping gene", title = "mean of replicates, Relative abundance") +
  facet_wrap(~target_id,nrow=1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
plot_rel_abun



#----
#calculate ΔΔCt against a chosen reference sample ie WT or time point 0

# first, rebuild tidyqpcr's 'calculate_deltadeltacq_bytargetid' function to 'calculate_deltadeltact_bytargetid'
calculate_deltadeltact_bytargetid <- function (delta_Ct_df, ref_sample_ids, norm_function = mean, ddct_positive = TRUE)
  {
ddct_factor <- (-1)^ddct_positive #output ΔΔCt as positive if target is detected more in the sample, compared to reference sample
delta_Ct_df %>% dplyr::group_by(.data$target_id) %>%
  dplyr::do(calculate_normvalue(.data, 
       ref_ids = ref_sample_ids, value_name = "delta_Ct", id_name = "sample_id", 
       norm_function = norm_function)) %>% dplyr::rename(ref_delta_Ct = .data$value_to_norm_by) %>% 
    dplyr::mutate(delta_delta_Ct = ddct_factor * (.data$delta_Ct - .data$ref_delta_Ct), #-ΔΔCt= -1* (mean ΔCt(sample) - mean ΔCt(reference sample)) for each target
                fold_change = 2^.data$delta_delta_Ct) ##adds column for fold-change; (2^-ΔΔCt)
    }

#then use it to calculate ΔΔCt for each sample: ΔΔCt = mean delta_Ct(sample) - mean delta_Ct(reference sample) for each target
#relative change in gene target ΔCt values compared to reference sample; globally normalises quantification cycle (log2-fold) data across sample_id
plate_delta_delta <- plate_norm %>%
     calculate_deltadeltact_bytargetid(ref_sample_ids = "WT 1")

plate_delta_delta_mean <- plate_delta_delta %>%
  group_by(sample_id, target_id) %>%
  summarize(
    delta_delta_Ct = mean(delta_delta_Ct, na.rm = TRUE),
    fold_change = mean(fold_change, na.rm = TRUE))

View(plate_delta_delta) #displays calculations appended to columns of the original file (normalised plate_norm)
View(calculate_deltadeltact_bytargetid(plate_norm_by_mean, ref_sample_ids = "WT 1", norm_function = mean)) #shows calculations only
#View(plate_delta_delta_mean) #also shows summarised table


#----

#now plot data!
ΔΔCt_plot <- ggplot(data = subset(plate_delta_delta_mean, target_id != "h20 only" & target_id != "rg1262/63 hkg control" & delta_delta_Ct != "NaN" #remove unnecesary samples and controls
            ), aes(x = sample_id, y = delta_delta_Ct, fill = target_id, shape = sample_id, label = sample_id)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1), colour = "black") +
  geom_point(position = position_dodge(width = 1)) +
  geom_text(position = position_dodge(width = 1), hjust = 1.1, size = 3, angle =90) +   # hjust=ifelse(plate_delta_delta_mean$delta_delta_Ct < max(plate_delta_delta_mean$delta_delta_Ct)/1.5, -0.1, 1.1)  #hjust, if ΔΔCt >0 then adjust text to 1.3, if <0 adjust to -1
  labs(y = "ΔΔCt (log2 fold change)\n relative to WT biological replicate 1") + 
  facet_wrap(~target_id, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
ΔΔCt_plot

fold_change_plot <- ggplot(data = subset(plate_delta_delta_mean, target_id != "h20 only" & target_id != "rg1262/63 hkg control" & delta_delta_Ct != "NaN" #remove unnecesary samples and controls
    ), aes(x = sample_id, y = fold_change, fill = target_id, shape = sample_id, label = sample_id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), colour = "black") + 
  geom_point(position = position_dodge(width = 1)) +
  geom_text(position = position_dodge(width = 1), vjust = -1, hjust = 1.1, size = 3) + 
  labs(y = "fold change, 2^(-ΔΔCt) relative to \n WT biological replicate 1") + 
  facet_wrap(~target_id, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
fold_change_plot









#-------------
#try something else.......!!








Quantity = qpcr_test15.1_tidy$Quantity
res <- pcr_assess(qpcr_test15.1_tidy,
            amount = Quantity,
            reference_gene = 'target_id',
            method = 'efficiency')
res


BiocManager::install("ddCt")







## add grouping variable by generating replicates of the values in x
group_var <- rep(c('c84 1', 'c84 2', 'QFS1 1', 'QFS1 2', 'WT 1'), each = 4) #there are 4 samples of c84 (two biological replicates for targets 7280; cas9-267), 4 samples of qfs1 two biological replicates for targets dsx/ dsxT1)

# calculate all values and errors in one step
#here we locate Ct values of the four target genes 7280;c84, cas9-267, dsx & dsxT1
#this obtains the expression of these targets normalised against the hkg-control relative to WT
#remember orginal data doesn;t contain cas9-267 as a target for c84
res <- pcr_analyze(qpcr_test15,
                   group_var = group_var,
                   reference_gene = 'rg1262/63 hkg control',
                   reference_group = 'WT1')

res

