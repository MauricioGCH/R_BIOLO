#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(seqinr)
library(Biostrings)
library(DT)

# 'util' file for genetic code logic and strings cleaning
source("util.R")

# server logic to read fasta files
shinyServer(function(input, output) {
  
  output$dna.content <- DT::renderDataTable({
    
    # save inFile for internal logic
    inFile <- input$dna.file
    
    if (is.null(inFile))
      return(NULL)
    
    dna.sequences <- read.fasta(inFile$datapath, seqtype = "DNA", forceDNAtolower = T)
    
    dna.content <- compute.dna.metrics(dna.sequences)
    
    return(dna.content)
  })
  
  output$aa.sequences <- renderUI({
    
    inFile.2 <- input$dna.file.2
    p.length <- input$protein.length
    if (is.null(inFile.2))
      return(NULL)
    
    dna.sequences <-  read.fasta(file = "~/RProjects/Tests/TP-Student/Students/Mysterious_seq.txt", 
                                 seqtype = "DNA", 
                                 forceDNAtolower = T)
      
      list.of.proteins <- find.mysterious.proteins(dna.sequences$seq0,p.length)
      
      HTML(proteins.to.html(list.of.proteins))

  
  
})
  output$align.out <- renderPrint({
    # Exercise 7.4 - Fill here: reference sequence given user's input, between PKD1 and PKD2
    
    sequence <- input$subject.seq
    if(input$select.ref.seq=="PKD1_protein"){
      ref.sequence <- readAAStringSet('~/RProjects/Tests/TP-Student/Students/FastaPKD1.txt')
    } 
    if(input$select.ref.seq=="PKD2_protein"){
      ref.sequence <- readAAStringSet('~/RProjects/Tests/TP-Student/Students/FastaPKD2.txt')
    }
      in.seq <- clean.sequence(input$subject.seq)
    if(is.null(in.seq))
      return(NULL)
      #Exercise 7.4 - Fill here: compute alignment with pairwiseAlignment given all the user parameters.
    pair.alignment <- pairwiseAlignment(pattern = ref.sequence,
                                        subject = sequence,
                                        substitutionMatrix = input$align.matrix,
                                        type = input$align.type,
                                        gapOpening = input$gOpening.length,
                                        gapExtension= input$gExtension.length)
      
      return(pair.alignment)
  })
  
})
