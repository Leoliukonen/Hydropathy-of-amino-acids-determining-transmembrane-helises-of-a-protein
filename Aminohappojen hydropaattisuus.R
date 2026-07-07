# Write a program that translates a DNA sequence from FASTA format to an amino acid (AH) sequence
# And calculates the sliding window moving average of the resulting protein's hydrophobicity

# Load libraries
library(seqinr)
library(dplyr)
library(ggplot2)

# Save sequence to a variable
sequence <- read.fasta("EGFR var 1.fasta") # Here the mRNA sequence of a transmembrane protein
sekv <- as.data.frame(sequence)
s <- as.vector(sekv[[1]])

# Reverse complement the sequence, IF NECESSARY
comp <- c()
for (x in 1:length(s)){
  väli <- case_when(
    s[x] == "a" ~ "t",
    s[x] == "t" ~ "a",
    s[x] == "c" ~ "g",
    s[x] == "g" ~ "c",
    TRUE ~ "Error"
  )
  comp <- c(comp,väli)
}

# Change thymine to uracil
snp <- c()
for (x in 1:length(s)){
  väli3 <- case_when(
    s[x] == "a" ~ "a",
    s[x] == "t" ~ "u",
    s[x] == "c" ~ "c",
    s[x] == "g" ~ "g",
    TRUE ~ "Error"
  )
  snp <- c(snp,väli3)
}

# Create a codon vector, first finding the start codon
# Initialize vectors
aloitus_kodonit <- c()
terminaatio_kodonit <- c()

for(x in seq(1:length(snp)-2)){ # Iterate through the entire sequence
  if(paste0(snp[x],snp[x+1],snp[x+2]) == "aug"){ # Find all start candidates
    
    v1 <- c(as.integer(x))
    aloitus_kodonit <- c(aloitus_kodonit,v1)
    for(y in seq(x,length(snp),by=3)){ # Find the termination codon after the aug-codon
      
      if(paste0(snp[y],snp[y+1],snp[y+2]) %in% c("uaa","uag","uga")){ # Check through termination codons
        
        v2 <- c(as.integer(y)) # Combine to intermediate variable
        terminaatio_kodonit <- c(terminaatio_kodonit,v2)
        break
      }}
    
  }
  
}

p1 <- length(aloitus_kodonit) # Determine lengths
p2 <- length(terminaatio_kodonit)
aloitus_kodonit <- head(aloitus_kodonit,p2) # Remove extras
# Determine frame lengths

framet <- data.frame(aloitus_kodonit,terminaatio_kodonit)
# Add length column to the framet dataframe
framet <- mutate(framet, pituus = terminaatio_kodonit - aloitus_kodonit)

# Select the longest frame
longest_frame <- subset(framet,pituus == max(framet$pituus))

# Create codon vector
kodoni <- c()

# Add codons to the codon vector using a loop
for (i in seq(longest_frame$aloitus_kodonit,longest_frame$terminaatio_kodonit, by=3)){
  v3 <- c()
  v3 <- paste0(snp[i],snp[i+1],snp[i+2])
  kodoni <- c(kodoni,v3)
}

# Create AH (amino acid) vector
AH <- c()

for (x in 1:length(kodoni)) {
  
  väli1 <- case_when(
    kodoni[x] %in% c("uuu", "uuc") ~ "F",
    kodoni[x] %in% c("uua", "uug", "cuu", "cuc", "cua", "cug") ~ "L",
    kodoni[x] %in% c("ucu", "ucc", "uca", "ucg","agu","agc") ~ "S",
    kodoni[x] %in% c("uau", "uac") ~ "Y",
    kodoni[x] %in% c("ugu", "ugc") ~ "C",
    kodoni[x] == "ugg" ~ "W",
    kodoni[x] %in% c("ccu", "ccc", "cca", "ccg") ~ "P",
    kodoni[x] %in% c("cau", "cac") ~ "H",
    kodoni[x] %in% c("caa", "cag") ~ "Q",
    kodoni[x] %in% c("aga", "agg","cgu","cgc","cga","cgg") ~ "R",
    kodoni[x] %in% c("auu", "auc", "aua") ~ "I",
    kodoni[x] == "aug" ~ "M", # Start codon
    kodoni[x] %in% c("acu", "acc", "aca", "acg") ~ "T",
    kodoni[x] %in% c("aau", "aac") ~ "N",
    kodoni[x] %in% c("aaa", "aag") ~ "K",
    kodoni[x] %in% c("guu", "guc", "gua", "gug") ~ "V",
    kodoni[x] %in% c("gcu", "gcc", "gca", "gcg") ~ "A",
    kodoni[x] %in% c("gau", "gac") ~ "D",
    kodoni[x] %in% c("gaa", "gag") ~ "E",
    kodoni[x] %in% c("ggu", "ggc", "gga", "ggg") ~ "G",
    kodoni[x] %in% c("uaa", "uag", "uga") ~ "STOP",
    TRUE ~ "Unknown" # If codon is not found
  )
  AH <- c(AH,väli1)
}
print(AH)

# Determine the hydrophobicity of the formed protein
hydropathy <- c()

for (x in(1:length(AH))){ # Numbers are from the Kyte-Doolittle scale
  väli4 <- case_when(
    AH[x] == "L" ~ 3.8,
    AH[x] == "V" ~ 4.2,
    AH[x] == "I" ~ 4.5,
    AH[x] == "F" ~ 2.8,
    AH[x] == "C" ~ 2.5,
    AH[x] == "M" ~ 1.9,
    AH[x] == "A" ~ 1.8,
    AH[x] == "G" ~ -0.4,
    AH[x] == "T" ~ -0.7,
    AH[x] == "S" ~ -0.8,
    AH[x] == "W" ~ -0.9,
    AH[x] == "Y" ~ -1.3,
    AH[x] == "P" ~ -1.6,
    AH[x] == "H" ~ -3.2,
    AH[x] == "E" ~ -3.5,
    AH[x] == "Q" ~ -3.5,
    AH[x] == "D" ~ -3.5,
    AH[x] == "N" ~ -3.5,
    AH[x] == "K" ~ -3.9,
    AH[x] == "R" ~ -4.5,
    TRUE ~ NA
  )
  hydropathy <- c(hydropathy,väli4)
}

# Calculate hydropathy in sliding windows of 16
ikkuna <- 16 # Window size for sliding average

rolling_means <- rep(NA,ikkuna*0.5) # Fill empty values with NA

for(i in 1:(length(hydropathy)-(ikkuna+1))){
  väli5 <- mean(hydropathy[i:(i+(ikkuna-1))])
  rolling_means <- c(rolling_means,väli5)
}
# Fill empty values with NA
fill <- rep(NA,ikkuna*0.5)
rolling_means <- c(rolling_means,fill)

# Create table with amino acid sequence and mean hydropathy
AH <- head(AH,-1)
OAA <- seq(1:length(AH))
hbaa <- data.frame(OAA,rolling_means)

# Generate plot
plot_hydropathy <- ggplot(hbaa,aes(OAA,rolling_means))+
  theme_bw(base_size = 15)+
  geom_point()+
  geom_line()+
  labs(x = "Amino acid order", y = "Moving average of hydropathy")
plot_hydropathy