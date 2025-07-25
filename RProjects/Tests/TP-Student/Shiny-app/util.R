library(seqinr)
library(shiny)
library(stringr)
library(Biostrings)
library(DT)

compute.dna.metrics <- function(dna.sequences){
  # Bonus exercise 7.5: combine the 6 sapply calls
  dna.content <- data.frame(names=names(dna.sequences))
  
  dna.content$a <- sapply(dna.sequences, compte, 'a')
  dna.content$c <- sapply(dna.sequences, compte, 'c')
  dna.content$g <- sapply(dna.sequences, compte, 'g')
  dna.content$t <- sapply(dna.sequences, compte, 't')
  
  dna.content$length <- sapply(dna.sequences, length)
  dna.content$GC.content <- sapply(dna.sequences, GC)
  
  return(dna.content)
}

compte <- function(seq, letter){
  count <- 0
  for (i in seq){ # iterate and count
    
    if (letter == i){
      count <- count + 1
    }
  }
  return(count)
}

codon.table <- read.table(file = "~/RProjects/Tests/TP-Student/Code_shortcut/Genetic_code.txt", 
                          col.names = c("codon", "aa", "letter"),
                          stringsAsFactors = F)
codon.table$codon <- tolower(codon.table$codon) # use only lower case characters
rownames(codon.table) <- codon.table$codon # use codons as table keys

##################################################################################################

complementary <- function(seq){
  comp <- c()
  for (i in seq){ # iterate, compare and change accordingly
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
#####################################################################

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
  codons <- paste0(sequence[c1],sequence[c1+1],sequence[c1+2])
  peptide <- codon.table[codons,"letter"]
  return(peptide)
}
##############################################
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

###################################################################################

clean.sequence <- function(dirty.seq){
  return(gsub('[^A-Z]', replacement = "", toupper(dirty.seq)))
}

proteins.to.html <- function(list.of.proteins){
  # Every 60 amino acids, insert a new line (html code)
  proteins.cut <- gsub("(.{61,}?)", "\\1</br>", list.of.proteins)
  # Use console typo in html
  proteins.in.span <- paste0("<p style='font-family:monospace'>",proteins.cut,"</span>")
  return(paste0(proteins.in.span, collapse = "</br></br>"))
}

# more to come