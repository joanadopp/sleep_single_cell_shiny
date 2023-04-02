library(shiny)
library(shinythemes)
library(shinyWidgets)
library(ggrepel)
library(tidyverse)
library(data.table)
library(R.utils)

LOCALPATH <- '/Users/Joana/app/'
print(list.dirs())
print(list.files())

if (Sys.info()['user']=='vibflysleep'){
  root_dir = LOCALPATH
} else {
  root_dir = "."
}

genes <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "zt_mms_long_sign.csv.gz")))$gene))
ct <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "zt_mms_long_sign.csv.gz")))$cluster))
genes1 <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "combined_nobatch_sign.csv.gz")))$gene))
ct1 <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "combined_nobatch_sign.csv.gz")))$cluster))
genes2 <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "ct_mms_sign.csv.gz")))$gene))
ct2 <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "ct_mms_sign.csv.gz")))$cluster))

ui <- navbarPage("Fly Sleep Single Cell", theme = shinytheme("flatly"),
                 
              tabPanel("about",
                       mainPanel(
                         h1("Explore expression of your gene of interest in your cluster of interest", style = "font-size:18px"),
                         h1("(1) across different circadian Zeitgeber times,", style = "font-size:18px"),
                         h1("(2) between sleep and wake states and", style = "font-size:18px"),
                         h1("(3) across different degrees of sleep pressure.", style = "font-size:18px"),
                         tags$br(),
                         tags$a(href="https://scope.aertslab.org/#/Fly_Brain_Sleep/Fly_Brain_Sleep%2FFly_Sleep.loom/welcome",
                                "Cell type annotations, gene expression and regulon activity are available on SCope (aertslab.org)", style = "font-size:18px"),
                         tags$br(),
                         h1("There are two versions of the shiny app. This one includes the expression of all significant genes in each of the 214 clusters across different circadian times, sleep/wake conditions, and levels of sleep pressure. 
                               The first version includes the expression of all 9995 genes in each of the 25 annotated clusters .", style = "font-size:18px"),
                         tags$a(href="https://joana-dopp.shinyapps.io/Fly_Sleep_Single_Cell_v1/",
                                "See version 1", style = "font-size:18px"),
                         tags$br(),
                         tags$br(),
                         tags$a(href="https://www.biorxiv.org/content/10.1101/2023.03.22.533150v1", "Prepint is available on bioRxiv", style = "font-size:18px")
                          )),
              
              tabPanel("circadian cyclers",
                       fluidRow(sidebarPanel(
                         #pickerInput("cell_type1",
                                     #"Choose cell type(s)",
                                     #choices = ct,
                                     #selected = 'EG_1',
                                     #multiple = TRUE,
                                     #options = list(`actions-box` = TRUE)),
                         
                         selectInput("transcript1",
                                     label = "Choose a transcript",
                                     choices = genes,
                                     selected = genes[1]))),
                       
                       plotOutput("lineplot", height = "500px", width = "1200px")
              ),
              
              tabPanel("sleep vs. wake",
                       fluidRow(sidebarPanel(
                         selectInput("cell_type2",
                                     label = "Choose a cell type",
                                     choices = ct1,
                                     selected = 'EG_1'),
                         
                         selectInput("transcript2",
                                     label = "Choose a transcript",
                                     choices = genes1,
                                     selected = genes1[1]))),
                       
                       plotOutput("volcano", width="1000px", height = "600px")
              ),
                       
              tabPanel("sleep drive",
                          tabsetPanel(
                            tabPanel("heatmap",
                                     selectInput("cell_type3",
                                                 label = "Choose a cell type",
                                                 choices = ct2,
                                                 selected = 'EG_1'),
                                     
                                     plotOutput("heatmap", height = "800px")),
                            
                            tabPanel("lineplot",
                                     #pickerInput("cell_type4",
                                                 #"Choose cell type(s)",
                                                 #choices = ct2,
                                                 #selected = 'EG_1',
                                                 #multiple = TRUE,
                                                 #options = list(`actions-box` = TRUE)),
                                     
                                     selectInput("transcript3",
                                                 label = "Choose a transcript",
                                                 choices = genes2,
                                                 selected = genes2[1]),
                                     
                                     plotOutput("lineplot2", height = "600px", width = "1200px"))
                          ))
)

server <- function(input, output) {
  #browser()
  
  zt_mms <- read.csv(gzfile(file.path(root_dir, "data/zt_mms_long_sign.csv.gz")))
  setDT(zt_mms)
  
  output$lineplot <- renderPlot({
    #ggplot(zt_mms[gene == input$transcript1 & cluster %in% input$cell_type1], aes(x=ZT, y=expression, color=cluster)) +
    ggplot(zt_mms[gene == input$transcript1], aes(x=ZT, y=expression, color=cluster)) +
      geom_line()+
      theme_classic(base_size=20)+
      ggtitle(input$transcript1)+
      ylab("mean normalized counts")+
      scale_x_continuous(breaks=c(2, 8, 14, 20))
    
  }, res=96)
  
  sleep_wake <- read.csv(gzfile(file.path(root_dir, "data/combined_nobatch_sign.csv.gz"))) 
  setDT(sleep_wake)
  
  output$volcano <- renderPlot({
    
    ggplot(sleep_wake[cluster %in% input$cell_type2], aes(x=logfoldchange, y=-log10(pval), col=pvalue)) +
      geom_point() + 
      geom_text_repel(data=sleep_wake[cluster %in% input$cell_type2] %>%
                        filter(gene %in% input$transcript2),
                      aes(label=gene), size=10, max.overlaps = 100, box.padding = 2) +
      ggtitle(input$transcript2)+
      scale_x_continuous(limits = c(-3, 3), name="log2 fold change") +
      scale_y_continuous(limits = c(0, 10), name="-log10 p-value") +
      theme_classic(base_size = 20)+
      guides(color = guide_legend(override.aes = aes(label = "")))
    
  }, res=96)
  
  ct_mms <- read.csv(gzfile(file.path(root_dir, "data/ct_mms_sign.csv.gz")))
  
  output$heatmap <- renderPlot({
    
    cluster_mms_sign <- subset(ct_mms, cluster == input$cell_type3, select=2:9)
    cluster_mms_sign <- data.frame(cluster_mms_sign[,-1], row.names = cluster_mms_sign[,1])
    cluster_mms_sign <- as.matrix(cluster_mms_sign)
    #z-score: for zscore calculation to work it needs to be matrix format with gene as the index
    cluster_mms_z <- t(apply(cluster_mms_sign, 1, scale))
    
    heatmap(cluster_mms_z,
            Colv = NA)
  }, res=96)
  
  colnames(ct_mms)[3] <- 1
  colnames(ct_mms)[4] <- 2
  colnames(ct_mms)[5] <- 3
  colnames(ct_mms)[6] <- 4
  colnames(ct_mms)[7] <- 5
  colnames(ct_mms)[8] <- 6
  colnames(ct_mms)[9] <- 7
  ct_mms_long <- gather(ct_mms, condition, expression, 3:9, factor_key=TRUE)
  ct_mms_long$condition <- as.numeric(as.character(ct_mms_long$condition))
  setDT(ct_mms_long)
  
  output$lineplot2 <- renderPlot({
    
    ggplot(ct_mms_long[gene == input$transcript3], aes(x=condition, y=expression, color=cluster)) +
      geom_line()+
      theme_classic(base_size=20)+
      ggtitle(input$transcript3)+
      ylab("mean normalized counts")+
      xlab("conditions ordered from low to high sleep drive")+
      scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7), labels=c("GBX-S ZT12-8", "S ZT12-20", "S ZT12-14", "SD ZT12-14", "SD ZT12-20", "SD ZT12-2", "SD ZT12-8"))+ 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  }, res=96)


}

shinyApp(ui = ui, server = server) 