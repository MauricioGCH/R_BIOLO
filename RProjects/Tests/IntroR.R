library("seqinr")
library("shiny")
library("stringr")
library(Biostrings)
library(DT)
###########################
# Function to count
compte <- function(characters, letter){
  count <- 0
  for (i in characters){
    
    if (letter == i){
      count <- count + 1
    }
  }
  return(count)
}

complementary <- function(seq){
  comp <- c()
  for (i in seq){
    if ("a" == i){
      comp <- append(comp,"t");next
    }
    if ("t" == i){
      comp <- append(comp,"a");next  
    }
    if ("c" == i){
      comp <- append(comp,"g");next
    }
    if ("g" == i){
      comp <- append(comp,"c");next
    }
  }
  return(comp)
}

cDNAtomRNA <- function(seq){
  comp <- c()
  for (i in seq){
    
    if ("a" == i){
      comp <- append(comp,"u");next
      
    }
    if ("t" == i){
      comp <- append(comp,"a");next
    }
    if ("c" == i){
      comp <- append(comp,"g");next
    }
    if ("g" == i){
      comp <- append(comp,"c");next
    }
  }
  comp<-rev(comp)
  return(comp)
}

dna2peptide <- function(sequence, frame=0, sens='F'){
  if(sens=='R') {
    sequence <- rev(complementary(sequence))
  }
  codon.table <- read.table(file = "~/RProjects/Tests/TP-Student/Code_shortcut/Genetic_code.txt", 
                            col.names = c("codon", "aa", "letter"),
                            stringsAsFactors = F)
  
  codon.table$codon <- tolower(codon.table$codon) 
  
  rownames(codon.table) <- codon.table$codon 
  
  l <- 3 * ((length(sequence) - frame)%/%3) 
  
  c1 <- seq(1,l,3) + frame
  codons <- paste0(sequence[c1],sequence[c1+1],sequence[c1+2]) # adding through the first codon, then the second one, the third and it repeats.
  peptide <- codon.table[codons,"letter"]
  return(peptide)
}  
find.mysterious.proteins <- function(sequence,p.length){
  
  proteins <- c()
  Contador <-0
  for(j in c(0,1,2)){ #For all the reading frames
    
    for(i in c("R","F")){#For all the reading frames
      
      exp <- paste(dna2peptide(sequence,frame = j,sens=i),collapse = "")
      
      pattern <- paste0('M[ACDEFGHIKLMNPQRSTVWY]{',p.length,',}(X|$)')
      
      match<-str_extract_all(exp, pattern) # extract matches that correspond to the pattern
      
      if(str_length(match)>=p.length){
        
        proteins<-append(proteins,match)
      }
    }
  }
  proteins<-unlist(proteins)
  write.table(proteins, "~/RProjects/Tests/TP-Student/mysterious_prot_output.txt", append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE) # save result in a .txt
  
  return(proteins)
}


################################################################################
#CHAPTER 7 ALLIGNMENT

library(Biostrings)


summary(ref.prot.sequences)
summary(target.prot.sequences)

Allignment <- function(my.sequence,ref.prot.sequences,Fname){
  pair.alignment <- pairwiseAlignment(pattern = ref.prot.sequences,
                                      subject = my.sequence,
                                      substitutionMatrix = "BLOSUM62",
                                      type = "global")
  path <-paste0("~/RProjects/Tests/TP-Student/Students/",Fname,".txt")
  writePairwiseAlignments(pair.alignment, file=path)
  
}


#target.prot.sequences <- readAAStringSet('~/RProjects/Tests/TP-Student/Students/Protein_sequences.txt')
#my.sequence <-target.prot.sequences$seq31
#ref.prot.sequences <- readAAStringSet('~/RProjects/Tests/TP-Student/Students/FastaPKD1.txt')
#Allignment(my.sequence,ref.prot.sequences,"alignment.seq31.PKD1.txt")


Htmap <- function(){
  target.prot.sequences <- readAAStringSet('~/RProjects/Tests/TP-Student/Students/Protein_sequences.txt')
  sMatrix <- matrix(0,32,32)
  for(i in 1:32){
    for(j in 1:32){
      pair.alignment <- pairwiseAlignment(pattern = target.prot.sequences[i],
                                          subject = target.prot.sequences[j],
                                          substitutionMatrix = "BLOSUM62",
                                          type = "global")
      sMatrix[i,j] <- pair.alignment@score
    }
  }
  heatmap(sMatrix, Rowv = NA, Colv = NA, symm = T, main = "Alignment scores heatmap for 32 protein sequences", scale = "column" )
}
#Htmap()




