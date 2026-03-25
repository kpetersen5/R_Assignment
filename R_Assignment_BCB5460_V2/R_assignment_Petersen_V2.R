
# Data Inspection 
###To start I imported both the fang_et_al_genotypes and snp_position into R

fang_et_al_genotype_URL <- "https://raw.githubusercontent.com/EEOB-BioData/BCB546_Spring2024/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt"
snp_position_URL <-"https://raw.githubusercontent.com/EEOB-BioData/BCB546_Spring2024/main/assignments/UNIX_Assignment/snp_position.txt"
fang_et_al_genotypes <- read.delim(fang_et_al_genotype_URL, header = TRUE, sep = "\t", quote = "")
snp_position<- read.delim(snp_position_URL, header = TRUE, sep = "\t", quote = "")

library(tidyverse) # making sure tidyverse is on

###Type of tells you the type of vector the data is. In this case it is list for both data files
typeof(fang_et_al_genotypes) #List
typeof(snp_position) #List

###To determine the type of list we can use the class function 
class(fang_et_al_genotypes) # data frame
class(snp_position) # data frame

###Length shows the length of the data file 
length(fang_et_al_genotypes) # 986
length(snp_position) # 15

###dim shows the dimensions of the data frame: number of rows and number of columns
dim(fang_et_al_genotypes) # 2782 986
dim(snp_position) # 983 15

###The names of the columns are shown below only for the snp_positions. Because there were 986 column names I did not list them all here but the first three are listed below
###- Sample_ID 
###- JG_OTU
###- Group
colnames(snp_position) ### [1] "SNP_ID"               "cdv_marker_id"        "Chromosome"           "Position"             "alt_pos"             
### [6] "mult_positions"       "amplicon"             "cdv_map_feature.name" "gene"                 "candidate.random"    
### [11] "Genaissance_daa_id"   "Sequenom_daa_id"      "count_amplicons"      "count_cmf"            "count_gene"          
colnames(fang_et_al_genotypes)

### The size of the files 
object.size(fang_et_al_genotypes) # 22681376 bytes
object.size(snp_position) # 327392 bytes

### I used the command glimpse instead of str to look at the data because the glimpse output is cleaner and easier to read.This shows a neat view of the data in the files. To observe the data delete eval=FALSE 
glimpse(fang_et_al_genotypes)
glimpse(snp_position)

### Also because we are working in R, once the data has been uploaded it is easy to look at the full data frame in table using the view command
view(fang_et_al_genotypes)
view(snp_position)

### This also helps with understanding what the variables in the data sets. Below I have described the variables for each of the data sets that is important for the assignment 
### fang_et_al_genotypes:
###  - Group --> This tells you what organism that data came from. for our assignment we are intersted in: 
###  - Maize group --> ZMMIL, ZMMLR, and ZMMMR
###  - Teosinte group --> ZMPBA, ZMPIL, and ZMPJA
###  - The rest of the column titles are the SNP_IDs with the genotype for that specific snp
### snp_position:  
###  - SNP_ID --> This is the snp identifier
###  - Chromosome --> This is the chromosome on which the SNP is located
###  - Position --> This is the location (bp) of the SNP on that chromosome


# Data Processing 
### First I wanted to pull out the Maize and Teosinte groups from the fang_et_al_genotypes data frame. I created new vectors titled Maize genotypes that contained the groups ZMMIL, ZMMLR, and ZMMMR and Teosinte_genotypes that contained the groups ZMPBA, ZMPIL, and ZMPJA
Maize_genotype <- fang_et_al_genotypes %>% filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR"))
Teosinte_genotype <- fang_et_al_genotypes %>% filter(Group %in% c("ZMPBA", "ZMPIL", "ZMPJA"))

### Next I wanted to transpose the mazie and teosinte vectors so they could be joined with the SNP_positions data. This created two new files, Mazie_genotype_tran and Teosinte_genotype_tran
Maize_genotype_tran <- t(Maize_genotype)
Teosinte_genotype_tran <- t(Teosinte_genotype)

### I was running into issues with manipulating the Maize_genotype_tran and Teosinte_genotype_tran. I used the the typeof() function to determine what type of vector my new data was. Using these commands I determined that after using the transpose function the data is converted to a character vector
typeof(Maize_genotype_tran)
typeof(Teosinte_genotype_tran)

### I used the as.data.frame() function to change my vectors back to data frames for easier manipulation 
Maize_genotype_tran <- as.data.frame(Maize_genotype_tran)
Teosinte_genotype_tran <- as.data.frame(Teosinte_genotype_tran)


### These steps are to fix the column names so they can be merged with the snp_positions file
colnames(Maize_genotype_tran) <- Maize_genotype_tran[1,] # Takes the data in the first row and makes that the column names
colnames(Teosinte_genotype_tran) <- Teosinte_genotype_tran[1,] # Takes the data in the first row and makes that the column names
Maize_genotype_tran_2 <- Maize_genotype_tran[-(2:3),] # Deleting rows 2&3
Teosinte_genotype_tran_2 <- Teosinte_genotype_tran[-(2:3),] # Deleting rows 2&3

### Now I need to adjust the row names for the maize and teosinte data frames so Sample_ID is a row name with the correct data in the column. This will be needed to merge the Maize and Teosinte data with the snp_position data. I used the mutate function to add a new column to the data frame. This moved the sample_ID names from being row titles to being data in the data frame
Maize_genotype_tran_3 <- Maize_genotype_tran_2 %>% mutate(SAMPLE_ID = rownames(Maize_genotype_tran_2))
Maize_genotype_tran_3 <- select(Maize_genotype_tran_3, SAMPLE_ID, everything() ) # Moving the sample_ID to the front
Teosinte_genotype_tran_3 <- Teosinte_genotype_tran_2 %>% mutate(SAMPLE_ID = rownames(Teosinte_genotype_tran_2))
Teosinte_genotype_tran_3 <- select(Teosinte_genotype_tran_3, SAMPLE_ID, everything()) # Moving the sample_ID to the front

### Cutting down the snp_position files to have only the data I need
snp_position_2 <- snp_position[ ,-(5:15)]

### I used the Merge function to merge the snp_position data with the Maize and Teosinte data. I merged the data by the sample_ID column from the Maize data and the snp_ID column from the snp_positions data. I also made sure to keep the other columns in the snp_positions data because I need the chromosome number later on
Maize_merge <- merge(Maize_genotype_tran_3, snp_position_2, by.x= "SAMPLE_ID", by.y= "SNP_ID", all.y=TRUE)
Maize_merge_2 <- select(Maize_merge, SAMPLE_ID, Chromosome, Position, everything()) # Moving the Chromosome and Position number to the front
Teosinte_merge <- merge(Teosinte_genotype_tran_3, snp_position_2, by.x= "SAMPLE_ID", by.y= "SNP_ID", all.y=TRUE)
Teosinte_merge_2 <- select(Teosinte_merge, SAMPLE_ID, Chromosome, Position, everything()) # Moving the Chromosome and Position number to the front

### Valerie Yu suggested this code to help with sorting my data numerically. Previsouly my data was being sorted alpha numerically 
typeof(Maize_merge_2$Position) # this produced Character
typeof(Teosinte_merge_2$Position) # this produced Character
Maize_merge_2$Position <- as.numeric(Maize_merge_2$Position)
Teosinte_merge_2$Position <- as.numeric(Teosinte_merge_2$Position)
typeof(Maize_merge_2$Position) # this produced Double
typeof(Teosinte_merge_2$Position) # this produced Double

colnames(Maize_merge_2)[colnames(Maize_merge_2) == "SAMPLE_ID"] <- "SNP_ID"
colnames(Teosinte_merge_2)[colnames(Teosinte_merge_2) == "SAMPLE_ID"] <- "SNP_ID"

### To create a function to iterate over all the chromosomes I wanted to see if both maize and teosinte had the same number of unique chromosomes. If the vectors weren't the same length then the iteration in the next step wouldn't work. I created a vector with all the unique maize and teosinte chromosomes
chr_maize <- unique(Maize_merge_2$Chromosome)  
chr_teosinte <- unique(Teosinte_merge_2$Chromosome) 

### I am using an if else statement to compare the lengths of the vectors. The addition of the sort function orders the chromosomes so they can be compared evenly. !identical means not identical
### - "Maize and Teosinte vectors are equal" was printed so I can continue to the function and iteration steps
if (!identical(sort(chr_maize), sort(chr_teosinte))) {
  print("Maize and Teosinte vectors are not equal... STOP") 
} else {
  print("Maize and Teosinte vectors are equal")
}

### Writing a function to be able to iterate the formation of the 40 necessary files. Create_files is the main function that I will use for lapply further down. Iterate_chr is the function that is doing the majority of the work for the iteration process
Create_files <- function(chr) { # I am using chr as the argument so the overall function will iterate over each chromosome
  Iterate_chr <- function(sorted_data, desc = FALSE) { # Here I am using sorted_data as the argument that is defined int eh following lines
    sorted_data %>%
      filter(Chromosome == chr) %>%
      arrange(if(desc) desc(Position) else Position) %>%
      mutate(if(desc) across(everything(), ~str_replace_all(as.character(.), "\\?", "-")))
  }
  maize_inc <- Iterate_chr(Maize_merge_2, FALSE)
  maize_dec <- Iterate_chr(Maize_merge_2, TRUE)
  teosinte_inc <- Iterate_chr(Teosinte_merge_2, FALSE)
  teosinte_dec <- Iterate_chr(Teosinte_merge_2, TRUE)
  
  write_csv(maize_inc, paste0("maize_chr_", chr, "_increase.csv"))
  write_csv(maize_dec, paste0("maize_chr_", chr, "_decrease.csv"))
  write_csv(teosinte_inc, paste0("teosinte_chr_", chr, "_increase.csv"))
  write_csv(teosinte_dec, paste0("teosinte_chr_", chr, "_decrease.csv"))
}

## The data_sorted is going to filter out only the chromosome column that matches the chr that is currently being iterated over. Then it is going to arrange by position by ascending or descending values. I used the desc=FALSE to aid in the process. Essentially if the logic returns a TRUE then it will sorted in descending order. If the logic returns a FALSE then it will sort in ascending order. Then using the mutate command and across(everything()) command I was able to change the ? to - for the decreasing position files. The \\? allows the ? to act like a character and not a "special character" in R. Now I am using the Iterate_chr function to create 4 different categories of the sorted data. Using the write_csv command I am taking the four file versions created above and saving them as each of the 40 files needed with the chr number changing depending on which chromosome is being worked on by the function. The paste0 command just combines the naming convention in the parenthesis

### Finally I used lapply to run Create_files for all the chromosomes. Using the invisible command prevents lapply from printing the results on to the console. It keeps things looking cleaner
### **The output of the lapply is the 40 individual CSV files**
invisible(lapply(1:10, Create_files))


# Data Visualization

### Putting the data tables in the long format instead of the wide format helped To make it tidy and easier to plot. Using the pivot_longer command it flipped the data so all of the sample_IDs are in one column with their respective SNP and genotypes. I had to use the mutate(across(everything)), as.character() to tell R to treat the entire data table as a character because previosuly there was a mix of characters and integers and R doesn't like that. I also kept the the chromosome and position columns when flipping the data table 
Maize_long <- Maize_merge_2 %>% mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = -c(SNP_ID, Chromosome, Position), names_to = "Sample_ID", values_to = "Genotype")
Teosinte_long <- Teosinte_merge_2 %>% mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = -c(SNP_ID, Chromosome, Position), names_to = "Sample_ID", values_to = "Genotype")

### SNP Distribution by Chromosome 
### To compare the difference in the number of SNPs across chromosomes for Maize and Teosinte I created a double bar chart. I created a data table with the number of SNPs per chromosome and then sorted that data table by chromosome number 
Maize_snp_count <- Maize_long %>% count(Chromosome) %>% 
  mutate(Chromosome = factor(Chromosome, levels = c(1:10, "multiple", "unknown"))) %>% 
  arrange(Chromosome) %>% rename("SNPs" = n)
Teosinte_snp_count <- Teosinte_long %>% count(Chromosome) %>% 
  mutate(Chromosome = factor(Chromosome, levels = c(1:10, "multiple", "unknown"))) %>% 
  arrange(Chromosome) %>% rename("SNPs" = n)

### I combined the SNPs count data tables for Maize and Teosinte by adding a species column. Then combined them using the rbind command.Then I again sorted them by chromosome number 
Maize_snp_count$Species <- "Maize"
Teosinte_snp_count$Species <- "Teosinte"
Combined_snp_count <- rbind(Maize_snp_count, Teosinte_snp_count)
Combined_snp_count$Chromosome <- factor(Combined_snp_count$Chromosome, levels = c(1:10, "multiple", "unknown"))

library(ggplot2) # making sure ggplot2 is on

### I used ggplot to create the double bar chart. I used stat = identity so that it uses the count that is already done. I also used position = dodge to put the bar next to each other instead of overlapping. **With the bar chart you can see that Maize has more SNPs per chromosome than Teosinte and the larger the larger the chromosome number the less amount of SNPs present**
ggplot(Combined_snp_count, aes(x = Chromosome, y = SNPs, fill = Species)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "SNP Distribution by Chromosome", x = "Chromosome", y = "Number of SNPs")

## Data Visualization #1 - SNP Distribution per Chromosome
Maize_long$Species <- "Maize"
Teosinte_long$Species <- "Teosinte"
Combined_long <- rbind(Maize_long, Teosinte_long)
Combined_long$Chromosome <- factor(Combined_long$Chromosome, levels = c(1:10, "multiple", "unknown"))
Combined_long <- Combined_long[Combined_long$Chromosome %in% 1:10,] # I took out the multiple and unknown chromosomes because they had no position data

### I needed to change the poisiton values to numerical data in R so it could be plotted
Combined_long$Position <- as.numeric(as.character(Combined_long$Position))

### I used ploted both Maize and Teosinte on teh plots again this time overlapped. I used the facet wrap command to that all the chromosomes had a different plot but where all still on the same axis for easy comparison. **Again we see that maize has more SNPs than Teosinte. Chromosome 7 is also interesting becuase there are no SNPs in the middle of the chromosome** 
ggplot(Combined_long, aes(x = Position, fill = Species)) +
  geom_histogram(binwidth = 1000000, position = "identity", alpha = 0.5) +
  facet_wrap(~Chromosome, scales = "free_x") +
  labs(title = "SNP Distribution per Chromosome", x = "Position (Mbp)", y = "Number of SNPs") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

## Data Visualization #2 - Missing data and Amount of Heterozygosity
### I created a new column that was labeled either Missing, Homozygous, or Heterozygous using the ifelse logic statements
Combined_long$Type <- ifelse(Combined_long$Genotype %in% c("?/?"), "Missing",
                      ifelse(Combined_long$Genotype %in% c("A/A", "G/G", "T/T", "C/C"), "Homozygous",
                      "Heterozygous"))

### I plotted the resulting data in a stacked bar chart to compare Maize and Teosinte. **While the highest proportions of SNPs came from homozygous genotypes, Teosinte has a higher proportion of heterozygous genotypes than Maize**
ggplot(Combined_long, aes(x = Species, fill = Type)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", title = "Proportion of Homozygous, Heterozygous, and Missing Sites")

## Data Visualization #3 - Distribution of Genotype Type by Chromosome
### I plotted the proportion of each type of SNP (heterozygous, homozygous, or Missing) for each chromosome for both Maize and Teosinte. This will show is there are any chromosomes that are extremely different in terms of the type of SNP for either Teosinte or Maize. **Based on the plots there was not any discernible trends that would warrant further investigation** 
ggplot(Combined_long, aes(x = Chromosome, fill = Type)) +
  geom_bar(position = "fill") +
  facet_wrap(~Species, scales = "free_x") + 
  labs(title = "Distribution of Genotype Types by Chromosome",
       x = "Chromosome",
       y = "Proportion",
       fill = "Genotype Status")

### END