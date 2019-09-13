library(plotly)
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(stringr)

# Weird hack as R sometimes "forgets" its working directory
wd <- setwd(".")
setwd(wd)

# Make an empty Google analytics file (for dev version - not for production)
if (!file.exists('google-analytics.js'))
{
  file.create('google-analytics.js')
}



collatedFilename <- 'civicmine_collated.tsv'
sentencesFilename <- 'civicmine_sentences.tsv'
civicdbFilename <- 'nightly-ClinicalEvidenceSummaries.tsv'

collatedFilename <- normalizePath(collatedFilename)
fileInfo <- file.info(collatedFilename)
modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]

collated <- fread(collatedFilename,sep='\t',header=T,stringsAsFactors=T)
collated[collated$variant_group=='','variant_group'] <- '[unknown]'

evidencetypes <- sort(unique(as.character(collated$evidencetype)))
geneNames <- sort(unique(as.character(collated$gene_normalized)))
cancerNames <- sort(unique(as.character(collated$cancer_normalized)))
drugNames <- sort(unique(as.character(collated$drug_normalized)))
variantNames <- sort(unique(as.character(collated$variant_group)))



civicdb <- fread(civicdbFilename,sep='\t',header=T)
civicdb <- civicdb[,c('pubmed_id','entrez_id','doid','drugs','evidence_type','variant')]
civicdb$drugs <- tolower(civicdb$drugs)
civicdb$variant <- tolower(civicdb$variant)

civicdbFileInfo <- file.info(civicdbFilename)
civicdbModifiedDate <- strsplit(as.character(civicdbFileInfo$mtime), ' ')[[1]][1]

civicdb$combined <- paste(civicdb$evidence_type,civicdb$entrez_id,civicdb$doid,civicdb$drugs,sep='_')
collated$combined <- paste(collated$evidencetype,collated$gene_entrez_id,gsub("DOID:","",collated$cancer_id),collated$drug_normalized,sep='_')

collated$in_civic <- factor(collated$combined %in% civicdb$combined, labels=c("No","Yes"))

sentences <- fread(sentencesFilename,sep='\t',header=T)
sentences$pubmed_link <- paste("<a target=\"_blank\" href='https://www.ncbi.nlm.nih.gov/pubmed/", sentences$pmid, "'>", sentences$pmid, "</a>", sep='')
sentences$paper_in_civic <- factor(sentences$pmid %in% civicdb$pubmed_id, labels=c("No","Yes"))


citationTableExplanation <- "<br /><br /><br /><b>Citation Table:</b><br />Select a biomarker in the table above to see associated citations and sentences<br /><br />"
lastModifiedText = paste("CIViCmine updated on ",modifiedDate,". Comparing with CIViC updated on ",civicdbModifiedDate,".", sep="")

input <- data.frame(gene_input='EGFR',cancer_input='non-small cell lung carcinoma',drug_input="",data_table_rows_selected=1, stringsAsFactors=F)

# Some cleanup
collated[,combined:=NULL]
civicdb <- NULL

ui <- function(req) {
  fluidPage(
    tags$head(
      includeHTML("google-analytics.js"),
      includeHTML("metadata.html"),
      tags$style(".rightAlign{float:right; margin-left:5px; margin-bottom: 20px;}")
    ),
    titlePanel("",windowTitle="CIViCmine"),
    helpText(includeHTML("header.html")),
    verbatimTextOutput("click"),
    tabsetPanel(type = "tabs",
                tabPanel("Browse", 
                  sidebarPanel(
                               checkboxGroupInput("incivic_input", "In CIViC", c('Yes', 'No')),
                               selectizeInput(inputId = "gene_input", 
                                              label=p("Gene",actionLink("gene_clear", " (Clear)", style='font-size:70%')), 
                                              choices = c('',geneNames), 
                                              selected = '', 
                                              multiple = FALSE, 
                                              options = list(maxOptions = 2*length(geneNames))),
                               selectizeInput(inputId = "cancer_input", 
                                              label=p("Cancer",actionLink("cancer_clear", " (Clear)", style='font-size:70%')), 
                                              choices = c('', cancerNames), 
                                              selected = '', 
                                              multiple = FALSE, 
                                              options = list(maxOptions = 2*length(cancerNames))),
                               selectizeInput(inputId="drug_input", 
                                              label=p("Drug",actionLink("drug_clear", " (Clear)", style='font-size:70%')), 
                                              choices=c('', drugNames), 
                                              selected = '', 
                                              multiple = FALSE, 
                                              options = list(maxOptions = 2*length(drugNames))),
                               checkboxGroupInput('evidencetype_input', 'Evidence Type', choices = evidencetypes),
                               checkboxGroupInput('variant_input', 'Variant', choices = variantNames),
                               actionLink("selectall","Select All"),
                               
                               width=2
                               #verbatimTextOutput("gene_text")
                             ),
                   mainPanel(
                     HTML("<p></p>"),
                     splitLayout(cellWidths = c("31%", "31%", "31%"), 
                                 plotlyOutput("gene_piechart"),
                                 plotlyOutput("cancer_piechart"),
                                 plotlyOutput("drug_piechart")),
                     
                     downloadButton("download_collated_all", label = "Download All", class='rightAlign'),
                     downloadButton("download_collated_shown", label = "Download Shown", class='rightAlign'),
                     
                     DT::dataTableOutput("data_table"),
                     
                     HTML(citationTableExplanation),
                     
                     downloadButton("download_sentences_all", label = "Download All Sentences", class='rightAlign'),
                     downloadButton("download_sentences_above", label = "Download Sentences for Biomarkers Above", class='rightAlign'),
                     downloadButton("download_sentences_shown", label = "Download Sentences for Selected Biomarker", class='rightAlign'),
                     
                     DT::dataTableOutput("sentences_table"),
                     
                     helpText(lastModifiedText)
                   )
                  
                ),
                tabPanel("Help", helpText(includeHTML("help.html"))),
                tabPanel("About", helpText(includeHTML("about.html")))
                )
    
    
  )
}



server <- function(input, output, session) {
  
  tableData <- reactive({
    selected <- rep(TRUE,nrow(collated))
    
    if (!is.null(input$incivic_input) && length(input$incivic_input)==1) {
      if (input$incivic_input == 'Yes') {
        selected <- selected & collated$in_civic=='Yes'
      } else if (input$incivic_input == 'No') {
        selected <- selected & collated$in_civic=='No'
      }
    }
    
    if (input$gene_input!='') {
      selected <- selected & collated$gene_normalized==input$gene_input
    }
    
    if (input$cancer_input!='') {
      selected <- selected & collated$cancer_normalized==input$cancer_input
    }
    
    if (input$drug_input!='') {
      selected <- selected & collated$drug_normalized==input$drug_input
    }
    
    if (!is.null(input$evidencetype_input)) {
      selected <- selected & collated$evidencetype %in% input$evidencetype_input
    }
    
    if (!is.null(input$variant_input)) {
      selected <- selected & collated$variant_group %in% input$variant_input
    }
    
    table <- collated[selected,]
    
    if (nrow(table)>0) {
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$data_table <- DT::renderDataTable({
    table <- tableData()
    DT::datatable(table[,c('evidencetype','gene_normalized','cancer_normalized','drug_normalized','variant_withsub','in_civic','citation_count')],
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Evidence Type','Gene', 'Cancer','Drug','Variant','In CIViC','Citation #'),
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  dataTableProxy = dataTableProxy('data_table')
  
  output$sentences_table <- DT::renderDataTable({
    if(length(input$data_table_rows_selected)>0) {
      table <- tableData()
      row <- table[input$data_table_rows_selected,]
      entries <- sentences[sentences$matching_id==row$matching_id,]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    DT::datatable(entries[,c('pubmed_link','journal_short','year','paper_in_civic','section','subsection','formatted_sentence')],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'Paper in CIVIC', 'Section', 'Subsection', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  
  
  
  
  gene_event_val <- reactiveValues(count = 0)
  gene_piechart_data <- reactive({
    table <- tableData()
    piecounts <- NULL
    if (nrow(table) > 0) {
      
      piecounts <- aggregate(table$citation_count,by=list(table$gene_normalized),FUN=sum)
      colnames(piecounts) <- c('label','total_citation_count')
      piecounts <- piecounts[order(piecounts$total_citation_count,decreasing=T),]
      
      total <- sum(piecounts$total_citation_count)
      cutoff <- (1.5/100.0) * total
      piecounts$label <- as.character(piecounts$label)
      piecounts[piecounts$total_citation_count<cutoff,'label'] <- 'other'
      piecounts$label <- factor(piecounts$label, levels=c('other',as.character(piecounts[piecounts$total_citation_count>=cutoff,'label'])))
      
      piecounts <- aggregate(piecounts$total_citation_count,by=list(piecounts$label),FUN=sum)
      colnames(piecounts) <- c('label','total_citation_count')
    }
    piecounts
  })
  
  output$gene_piechart <- renderPlotly({
    piecounts <- gene_piechart_data()
    if (!is.null(piecounts)) {
      sourceName <- paste('gene_piechart_',gene_event_val$count,sep='')
      p <- plot_ly(piecounts, labels = ~label, values = ~total_citation_count, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
        layout(title = paste('Genes'),
               showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('gene_piechart_',gene_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    piecounts <- gene_piechart_data()
    if (length(d) > 0 && !is.null(piecounts))
    {
      selected <- piecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "gene_input","Gene", c("",geneNames), selected=selected)
        gene_event_val$count <- gene_event_val$count + 1
      }
    }
  })
  
  
  cancer_event_val <- reactiveValues(count = 0)
  cancer_piechart_data <- reactive({
    table <- tableData()
    piecounts <- NULL
    if (nrow(table) > 0) {
      
      piecounts <- aggregate(table$citation_count,by=list(table$cancer_normalized),FUN=sum)
      colnames(piecounts) <- c('label','total_citation_count')
      piecounts <- piecounts[order(piecounts$total_citation_count,decreasing=T),]
      
      total <- sum(piecounts$total_citation_count)
      cutoff <- (1.5/100.0) * total
      piecounts$label <- as.character(piecounts$label)
      piecounts[piecounts$total_citation_count<cutoff,'label'] <- 'other'
      piecounts$label <- factor(piecounts$label, levels=c('other',as.character(piecounts[piecounts$total_citation_count>=cutoff,'label'])))
      
      piecounts <- aggregate(piecounts$total_citation_count,by=list(piecounts$label),FUN=sum)
      colnames(piecounts) <- c('label','total_citation_count')
    }
    piecounts
  })
  
  output$cancer_piechart <- renderPlotly({
    piecounts <- cancer_piechart_data()
    if (!is.null(piecounts)) {
      sourceName <- paste('cancer_piechart_',cancer_event_val$count,sep='')
      p <- plot_ly(piecounts, labels = ~label, values = ~total_citation_count, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
        layout(title = paste('Cancers'),
               showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('cancer_piechart_',cancer_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    piecounts <- cancer_piechart_data()
    if (length(d) > 0 && !is.null(piecounts))
    {
      selected <- piecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "cancer_input","Cancer", c("",cancerNames), selected=selected)
        cancer_event_val$count <- cancer_event_val$count + 1
      }
    }
  })
  
  
  
  
  drug_event_val <- reactiveValues(count = 0)
  drug_piechart_data <- reactive({
    table <- tableData()
    piecounts <- NULL
    if (nrow(table) > 0) {
      table$drug_normalized <- as.character(table$drug_normalized)
      table[table$drug_normalized=='','drug_normalized'] <- 'N/A'
      
      piecounts <- aggregate(table$citation_count,by=list(table$drug_normalized),FUN=sum)
      colnames(piecounts) <- c('label','total_citation_count')
      piecounts <- piecounts[order(piecounts$total_citation_count,decreasing=T),]
      
      total <- sum(piecounts$total_citation_count)
      cutoff <- (1.5/100.0) * total
      piecounts$label <- as.character(piecounts$label)
      piecounts[piecounts$total_citation_count<cutoff,'label'] <- 'other'
      piecounts$label <- factor(piecounts$label, levels=c('other',as.character(piecounts[piecounts$total_citation_count>=cutoff,'label'])))
      
      piecounts <- aggregate(piecounts$total_citation_count,by=list(piecounts$label),FUN=sum)
      colnames(piecounts) <- c('label','total_citation_count')
    }
    piecounts
  })
  
  output$drug_piechart <- renderPlotly({
    piecounts <- drug_piechart_data()
    if (!is.null(piecounts)) {
      sourceName <- paste('drug_piechart_',drug_event_val$count,sep='')
      p <- plot_ly(piecounts, labels = ~label, values = ~total_citation_count, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
        layout(title = paste('Drugs'),
               showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('drug_piechart_',drug_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    piecounts <- drug_piechart_data()
    if (length(d) > 0 && !is.null(piecounts))
    {
      selected <- piecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "drug_input","Drug", c("",drugNames), selected=selected)
        drug_event_val$count <- drug_event_val$count + 1
      }
    }
  })
  
  output$download_collated_all <- downloadHandler(
    filename = function() {
      return("civicmine_collated.tsv")
    },
    content = function(file) {
      outdata <- collated
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_collated_shown <- downloadHandler(
    filename = function() {
      return("civicmine_collated_subset.tsv")
    },
    content = function(file) {
      outdata <- tableData()
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_sentences_all <- downloadHandler(
    filename = function() {
      return("civicmine_sentences.tsv")
    },
    content = function(file) {
      write.table(sentences, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_sentences_shown <- downloadHandler(
    filename = function() {
      return("civicmine_sentences_selectedbiomarker.tsv")
    },
    content = function(file) {
      if(length(input$data_table_rows_selected)>0) {
        table <- tableData()
        row <- table[input$data_table_rows_selected,]
        entries <- sentences[sentences$matching_id==row$matching_id,]
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_sentences_above <- downloadHandler(
    filename = function() {
      return("civicmine_sentences_multiplebiomarkers.tsv")
    },
    content = function(file) {
      table <- tableData()
      entries <- sentences[sentences$matching_id %in% table$matching_id,]
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )

  
  observeEvent(input$gene_clear, {
    updateSelectizeInput(session, "gene_input", selected = F)
  })  
  observeEvent(input$cancer_clear, {
    updateSelectizeInput(session, "cancer_input", selected = F)
  })  
  observeEvent(input$drug_clear, {
    updateSelectizeInput(session, "drug_input", selected = F)
  })  
}

shinyApp(ui, server)
