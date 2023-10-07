

####################################################################################
############  WELCOME TO CHRIS MULLIGAN'S FASTA FILE PROCESSING SCRIPT  ############
#########################  VERSION  2.2, 10/06/2023  ###############################
####################################################################################


# This R script does two things:
# 1. Cleans up sequence names from a fasta file and outputs an updated version
# 2. Appends taxonomic data from NCBI, with notes on how to get this data using a shell script/command line
# For best results, do this with the original complete fasta file (the script can handle multiple NCBI databases too)
      # (i.e. one big file with >seqname on one line then peptides on the next line)
# Some lines have been hard coded for a metagenomic dataset, but they can be removed as needed

# The script uses logical arguments to separate out species, accession numbers, etc.
# But it is by no means accurate, it simply cuts down the amount of manual edits needed
# DOUBLE CHECK YOUR FINAL OUTPUT FILES!



#First set up the work space:
### not run by TB ##dir.create("~/Desktop/R_cleanup") #This makes a new folder on your desktop to keep things clean
setwd("~/Sandbox/23_Rhizaria_SMTs/R_Annotate") 
# Move your fasta file in this folder

#Write your file name in the quotes (including the .fasta)
myFileName <- "7_Rhizaria_trimV3.copy.fasta"
# myFileName <- "filename.fasta"



# Before running, this script assumes that:
# 1. The sequence header contains the species name
# 2. The only thing ever in [brackets] is a species name, but species names outside of brackets are fine
# 3. If you've noted the outgroups, it's as |OUTGROUP at the end of the seqeunce header



# Did you modify the file names from the original NCBI output?
# 'Yes' changes the assumed sequence header to Accession|Species
fastaCondition <- ("No") #Write 'Yes' or 'No'

# Do you want the final sequence headers to be:
# 'Option_1' Species|Phylum_Class|Metadata|Accession
# 'Option_2' Accession|*ALL* taxonomy from NCBI
# 'Option_3' Species|LargeClade_Phylum_Family|Accession
TaxonomyCondition <- ("Option_3") #Write 'Option_1', 'Option_2', or 'Option_3'




# Most things are automated, i.e. you don't need to change any of the code beyond this point
# You will need to do some manual work in Section 2 in a web browser and command line





#Library packages ---
{
  require(Biostrings) #This package will take a while to install, go make some tea
  require(dplyr)
  require(reshape2)
  require(seqinr)
  require(stringr)
  require(tidyr)
  require(tidyverse)
}


 #Section 1: Format the pep fasta file ---------------------------------------------------------
{
fastaFile <- readDNAStringSet(myFileName) #Ignore warning

seq_name = names(fastaFile)
sequence = paste(fastaFile)
fasta.df <- data.frame(seq_name, sequence)

seqSplit <- as.data.frame(fasta.df$seq_name)
colnames(seqSplit) <- "seq_name"
}

{
#Isolate the species name ---
seqSplit$seq_name <- ifelse(grepl("\\|", seq_name), paste0(seqSplit$seq_name), 
                            paste0("|", seqSplit$seq_name))
seqSplit <- separate(seqSplit, col = seq_name, into = c('Bar', 'Etc'), sep='\\|')
seqSplit <- separate(seqSplit, col = Etc, into = c('PreBrac', 'Species_PostBrac'), sep='\\[')
seqSplit <- separate(seqSplit, col = Species_PostBrac, 
                     into = c('BracContents', 'PostBrac'), sep='\\]')

seqSplit[is.na(seqSplit)] <- ""

seqSplit[seqSplit == ''] <- NA
seqSplit$Species <- ifelse(is.na(seqSplit$BracContents), paste0(seqSplit$Bar), paste0(seqSplit$BracContents))

seqSplit$Bar <- NULL
seqSplit$BracContents <- NULL

seqSplit$Original <- fasta.df$seq_name


seqSplit$Species <- gsub("_metagenome", "", seqSplit$Species)
seqSplit$Species <- gsub("_metagenomics", "", seqSplit$Species)

seqSplit$Temp_Conditional <- ifelse(fastaCondition=="Yes", paste0("Yes"), paste0("No"))
seqSplit$Species <- ifelse(seqSplit$Temp_Conditional=="Yes", paste0(seqSplit$PreBrac), paste0(seqSplit$Species))

seqSplit$Temp_Conditional <- NULL
}

unique(seqSplit$Species) 
print("Are these all species? If not change it in your .fasta and start over or change the script")


#Isolate the accession number ---
{
TempAc <- as.data.frame(seqSplit$Original)
colnames(TempAc) <- "Original"
TempAc$Original <- gsub("WP_", "WP!", TempAc$Original) 
TempAc$Original <- gsub("XP_", "XP!", TempAc$Original)
rownames(TempAc) <- paste0(1:nrow(TempAc), "!")

TempAc$Original <- gsub("\\s", "_", TempAc$Original) 
TempAc <- TempAc %>% separate(Original, into = c('a', 'b', 'c'), sep= "\\|")
TempAc <- TempAc %>% separate(col = a, into = c('d', 'e', 'f', 'g'), sep= "_")
TempAc <- TempAc %>% separate(col = b, into = c('h', 'i', 'j', 'k'), sep= "_")
TempAc <- TempAc %>% separate(col = c, into = c('l', 'm', 'n', 'o'), sep= "_")


TempAc <- data.frame(lapply(TempAc, function(x) {gsub("!", "_", x)}))
TempAc$index <- gsub("!", "", rownames(TempAc))
TempAc$index <- sapply(TempAc$index, as.integer)

AccData <- TempAc %>% 
  pivot_longer(-index, names_to = NULL, values_to = "value") %>% 
  filter(str_count(value, "[0-9]") >= 5) %>% 
  full_join(TempAc %>% select(index), by = "index")
AccData <- as.data.frame(AccData)

AccData <- AccData %>% filter(grepl('[a-zA-Z]', value))

AccData <- AccData %>% 
  group_by(index) %>%
  mutate(value = toString(value)) %>%
  as.data.frame()

AccData <- AccData[!duplicated(AccData$index), ]
AccData$value <- gsub(", ", "_", AccData$value)

colnames(AccData) <- c("index", "Accession")
}

#Isolate/generate metadata ---
{
TempMeta <- as.data.frame(seqSplit$Original)
colnames(TempMeta) <- "Original"
TempMeta$index <- 1:nrow(TempMeta)
TempMeta$metag <- ifelse(grepl("icrobiom", TempMeta$Original), paste0("Microbiome"), 
                         ifelse(grepl("etagenom", TempMeta$Original), paste0("Metagenome"), NA))
TempMeta$env <- ifelse(grepl("hydrothermal|sediment|marine", 
                             TempMeta$Original), paste0("Environment"), NA)

TempMeta$Metadata <- ifelse(is.na(TempMeta$env), paste0(TempMeta$metag), paste0(TempMeta$env))
TempMeta$Metadata <- gsub("NA", NA, TempMeta$Metadata)

TempMeta$Human_value <- ifelse(grepl("Human", TempMeta$Original), paste0("Human_microbiome"), NA)
TempMeta$Human_value <- gsub("NA", NA, TempMeta$Human_value)

TempMeta$Metadata <- ifelse(is.na(TempMeta$Human_value), 
                            paste0(TempMeta$Metadata), paste0(TempMeta$Human_value))
TempMeta$Metadata <- gsub("NA", NA, TempMeta$Metadata)

}

#Clean up new dataframe ---
{
seqSplit$index <- 1:nrow(seqSplit)
seqSplit <- full_join(seqSplit, AccData, by = "index")
seqSplit <- full_join(seqSplit, TempMeta, by = "index")

colnames(seqSplit)

seqSplit$genus <- gsub("_.*", "", seqSplit$Species)
seqSplit <- seqSplit[,colSums(is.na(seqSplit)) < nrow(seqSplit)]
}

#Section 2: Getting taxonomy data from NCBI --------------------------------------------------
{
just_Genera <- unique(seqSplit$genus)
just_Genera <- gsub("_", " ", just_Genera)
just_Genera <- gsub(" sp.", "", just_Genera)
just_Genera <- paste0(just_Genera, ",")

write.table(just_Genera, file = "just_Genera.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
print("Follow the instructions outlined below")
}

#Paste contents of just_Genera.txt into this link (without changing any defaults):
#https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

#Save the results as tax_report.txt in the folder R_Annotate
#To save the file click 'Save in file' below the text box

{
NCBI_report <- read.delim("tax_report.txt", header = TRUE, sep = "\t")
NCBI_report <- NCBI_report[!is.na(NCBI_report$taxid),]
NCBI_report$taxid <- as.numeric(NCBI_report$taxid)

NCBI_IDs <- NCBI_report$taxid
NCBI_IDs <- paste0(NCBI_IDs, ",")
NCBI_IDs <- t(NCBI_IDs)

write.table(NCBI_IDs, file = "NCBI_taxIDs.txt", sep = "", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
print("Follow the instructions outlined below")
}

#Run the installation steps outlined here: https://www.ncbi.nlm.nih.gov/books/NBK179288/

#Create a file called taxid_to_taxonomy_v1.sh in a text editor with this code and save it to ~/Sandbox/R_Annotate/
#  efetch -db taxonomy -id MY_LIST -format xml | \
#  xtract -pattern Taxon -first TaxId -element Taxon -block '*/Taxon'  \
#  -unless Rank -equals 'no rank' -tab ',' -sep '_' -element Rank,ScientificName

{
my_Shell <- readr::read_delim("taxid_to_taxonomy_v1.sh", quote = "", delim = "\t", col_names = FALSE)
my_Shell

my_taxID_List <- paste(NCBI_IDs, collapse = " ")
my_taxID_List <- str_sub(my_taxID_List, end = -2)
my_taxID_List

new_Shell <- data.frame(lapply(my_Shell, function(x) {gsub("MY_LIST", paste0(my_taxID_List), x)}))
new_Shell

write_delim(new_Shell, file = "taxid_to_taxonomy_v2.sh", 
            col_names = FALSE, delim = "", quote_escape = "none")
print("Follow the instructions outlined below")
}

#Run the following in command line:
#  cd ~/Sandbox/R_Annotate/
#  export PATH=${PATH}:${HOME}/edirect
#  bash taxid_to_taxonomy_v2.sh > taxonomy_output.csv



#Section 3: Format the taxonomy .csv --------------------------------------------------------
{
taxID.df <- read.csv("taxonomy_output.csv", header = F, col.names = paste("V",1:50))
taxID.df <- taxID.df[, colSums(is.na(taxID.df)) < nrow(taxID.df)]
taxID.df$taxID <- taxID.df$V.1
taxID.df$taxID <- sub("[[:space:]].*", "", taxID.df$taxID)

taxID.df$V.2 <- gsub("clade_", "largeClade_", taxID.df$V.2)

taxID.all <- taxID.df
taxID.all$taxID <- paste0("!", taxID.all$taxID)
taxID.all <- taxID.all %>% mutate(across(everything(), ~ gsub(".*_", "", .)))
taxID.all$Mashed <- unite(taxID.all, Mashed)

MashedDF <- taxID.all$Mashed
MashedDF <- MashedDF %>% separate(col = Mashed, into = c('Mashed', 'taxID'), sep= "!")
MashedDF$Mashed <- gsub("[[:punct:]]+",'_', MashedDF$Mashed)
MashedDF$Mashed <- sub("_$", "", MashedDF$Mashed)

melt_taxa <- melt(taxID.df, id = c("taxID"))
melt_taxa$value <- gsub("[[:digit:]]+", "", melt_taxa$value)
melt_taxa$value <- gsub(" ", "", melt_taxa$value) 
melt_taxa$variable <- NULL
melt_taxa <- melt_taxa %>% separate(value, sep = "_", into = c("level", "data"))
melt_taxa <- melt_taxa[!duplicated(melt_taxa[c(1,2)]),]

cast_taxa <- dcast(melt_taxa, taxID ~ level)

taxa.DF <- cast_taxa[, c("largeClade", "phylum", "class", "order", "family", "genus", "taxID")]

taxa.DF$phylum <- ifelse(is.na(taxa.DF$phylum), paste0(taxa.DF$class), taxa.DF$phylum)
taxa.DF$taxID <- as.character(taxa.DF$taxID)

taxa.DF <- dplyr::full_join(taxa.DF, MashedDF, by = "taxID")
}

#Section 4: Combine and format -------------------------------------------------------------
{
shortID <- NCBI_report[, c("name", "taxid")]
shortID$name <- gsub(",", "", shortID$name)
colnames(shortID) <- c("genus", "taxID")

temp_split_ID <- full_join(seqSplit, shortID, by = "genus", relationship = "many-to-many")
temp_split_ID$taxID <- as.character(temp_split_ID$taxID)

with_taxa <- dplyr::full_join(temp_split_ID, taxa.DF, by = "taxID")



with_taxa$TempConditionFinal <- ifelse(TaxonomyCondition=="Option_1", paste0("Option_1"), 
                                       ifelse(TaxonomyCondition=="Option_2", paste0("Option_2"), 
                                              paste0("Option_3")))


with_taxa$final.name <- ifelse(with_taxa$TempConditionFinal=="Option_1", 
  paste(with_taxa$Species, paste(with_taxa$phylum, with_taxa$class, sep = "_"), 
        with_taxa$Metadata, with_taxa$Accession, sep = "|"), 
  ifelse(with_taxa$TempConditionFinal=="Option_2", 
         paste(with_taxa$Accession, paste(with_taxa$Mashed, with_taxa$Species, sep = "_"), sep = "|"), 
         paste(with_taxa$Species, paste(with_taxa$largeClade, with_taxa$phylum, with_taxa$family, sep = "_"), 
               with_taxa$Accession, sep = "|")))
  
with_taxa$TempConditionFinal <- NULL


with_taxa$final.name <- paste(">", with_taxa$final.name, sep = "")
with_taxa$final.name <- gsub("\\|NA", "", with_taxa$final.name)
with_taxa$final.name <- gsub("_NA", "", with_taxa$final.name)

with_taxa$final.name <- gsub("(\\|)+",'\\|', with_taxa$final.name)


with_taxa$final.name <- ifelse(grepl("outgroup", with_taxa$Original.y, ignore.case=TRUE), 
                               paste0(with_taxa$final.name, "|OUTGROUP"), 
                               paste0(with_taxa$final.name))

names(with_taxa)[names(with_taxa) == 'Original.x'] <- 'seq_name'
with_taxa <- dplyr::full_join(with_taxa, fasta.df, by = "seq_name", relationship = "many-to-many")


filename <- gsub("\\..*", "", myFileName)
}



#Some preliminary data visualizations ---
barplot(table(with_taxa$phylum), main = paste0(myFileName, " sample phylum frequency"), 
         col = "#800080", las = 1, horiz = TRUE, cex.names = 0.7)
hist(nchar(with_taxa$sequence), main = paste0(myFileName, " sequence read length"), 
     xlab = "length", col = "#800080")


write.csv(x = with_taxa, file = paste0(filename, "_allData.csv"))


#Section 5: Create final file and export --------------------------------------------------
{
output <- with_taxa[, c("final.name", "sequence")]

output <- output %>%
  arrange(final.name) %>%
  filter(!duplicated(sequence))

final <- data.frame(data = c(t(output)))

final <- data.frame(lapply(final, function(x) {gsub(" ", "_", x)}))


write.table(final, file = paste0(filename, "_cleaned.fasta"), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
}



#Package citations ---
  citation("Biostrings")
  citation("dplyr")
  citation("reshape2")
  citation("seqinr")
  citation("stringr")
  citation("tidyr")
  citation("tidyverse")
  print("Please credit Chris too, lol")
