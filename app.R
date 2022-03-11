
library(shiny)
library(tidyverse)
library(rmarkdown)
library(openxlsx)
library(readxl)
library(dqshiny)
library(ggpubr)
library(here)

## ---
library(showtext)
library(extrafont)

font_add_google("Montserrat", "Mont", regular.wt = 500)
showtext_auto()

my_theme <- theme_minimal() +
  theme(plot.title = element_text(family = "Mont", size = 22, colour = "#595959"),
        plot.subtitle = element_text(family = "Mont", size = 18, colour = "#595959"),
        legend.text = element_text(family = "Mont", size = 16, colour = "#595959"),
        axis.title = element_text(family = "Mont", size = 20, colour = "#595959"),
        strip.text = element_text(family = "Mont", size = 18, colour = "#595959"))


## ----

eDat <- readRDS(here::here("placenta2_miRM2_all_sncRNAs.RDS"))
pDat <- read_excel(here::here("snc_placenta_rob_pDat.xlsx"))
accessions <- read_excel(here::here("all_mir_accessions.xlsx"))

pDat2 <- pDat %>% 
  dplyr::select(sample, ID, caseID, tissue, condition, GA, trimester, sex, processingTime) 
pDat2$GA <- as.numeric(pDat2$GA)
pDat2$processingTime <- as.numeric(pDat2$processingTime)

eDat <- eDat %>% 
  filter(sncRNA == "miRNA") %>% 
  separate(Reference, c("precursor", "mature"), sep = ":") %>% 
  dplyr::select(-sncRNA)%>% 
  mutate(sum = rowSums(select_if(., is.numeric), na.rm = FALSE)) %>%
  group_by(mature) %>%
  slice_max(sum) %>%
  dplyr::select(-sum, -precursor) %>%
  distinct(mature, .keep_all = TRUE)

colnames(eDat)[2:31] <- pDat2$ID

precursor <- accessions %>% 
  dplyr::select(miRNA, precursor) %>% 
  distinct(miRNA, .keep_all = TRUE)
precursor <- eDat %>% 
  left_join(precursor, by = c("mature" = "miRNA")) %>% 
  dplyr::select(mature, precursor, everything())
  
eDat$mature <- str_sub(eDat$mature, 5)
precursor$mature <- str_sub(precursor$mature, 5)
precursor$precursor <- str_sub(precursor$precursor, 5)

e_plot <- eDat %>% 
  pivot_longer(cols = -c(mature), names_to = "ID", values_to = "RPM") %>% 
  inner_join(pDat2, by = "ID")

e_precursor <- precursor %>% 
  pivot_longer(cols = -c(mature, precursor), names_to = "ID", values_to = "RPM")  %>% 
  inner_join(pDat2, by = "ID") 

mirs <- eDat$mature

### -----

ui <- fluidPage(
  titlePanel("Expression of miRNAs in the Human Placenta"),
  sidebarLayout(
    sidebarPanel(
      autocomplete_input("mature", "Name of miRNA", mirs, value = "miR-210-3p", placeholder  = "miR-210-3p / let-7a-2-3p",)),
    mainPanel(h3(textOutput("mature_mir")),
      plotOutput("DotPlot"),
      br(),
      h3(textOutput("precursor_mir"),
      # plotOutput("Boxplots"),
      br(),
      plotOutput("bothMirsPlot"),
      br(),
      plotOutput("bothMirsSex"),
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$mature_mir <- renderText({
    paste("hsa-", input$mature, sep = "")
  })


### adding pre_mir as a reactive object
  
  pre_mir <- reactive({
    
    pre_mir <- e_precursor %>%
      filter(mature == input$mature)
    pre_mir <- pre_mir[1,2]
    pre_mir <- pre_mir$precursor
    
  })

  output$precursor_mir <- renderText({
    paste("Expression by trimester and sex for both hsa-", pre_mir(), "-5p", " and hsa-", pre_mir(), "-3p", sep = "")
  })
  
### dotplot expression for selected mature miRNA    
  output$DotPlot <- renderPlot({
    e_plot %>%
      filter(mature == input$mature) %>%
      dplyr::arrange(GA) %>%
      ggplot(aes(x = GA , y = log2(RPM + 1), colour = trimester)) +
      geom_point(size = 2.5) +
      geom_smooth(method = "lm", colour = "#666666", se = FALSE) + #doesn't work because x axis is categorical
      scale_colour_manual(values = c("#ecb3da", "#c36cac", "#981580")) +
      my_theme +
      # theme(axis.text.x = element_blank()) +
      theme(axis.text.x = element_text(size = 11)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0), breaks = seq(0, 40, 1)) +
      coord_cartesian(x = c(0,40)) +
      labs(x = "", y = "RPM\n", title = "\nSamples arranged in increasing order of gestational age", subtitle = "n = 30" , colour = "Trimester")
  })
  
### boxplot trimester expression for both precursor mirnas for selected mature mirna 
  
  output$bothMirsPlot <- renderPlot({
    
    e_precursor %>% 
      filter(precursor %in% pre_mir()) %>% 
      dplyr::arrange(GA) %>% 
      ggplot(aes(x = trimester, y = log2(RPM + 1))) +
      geom_boxplot(aes(fill = trimester), colour = "#595959", alpha = 0.9) +
      scale_fill_manual(values = c("#ecb3da", "#c36cac", "#981580")) +
      my_theme +
      theme(axis.text.x = element_blank()) +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "RPM\n", fill = "Trimester", x = "", title = "\nExpression by Trimester\nFor both mature -3p and -5p miRNAs\n") +
      facet_grid(~mature)
      
  })

  
### boxplot sex expression for both precursor mirnas for selected mature mirna  
  
  output$bothMirsSex <- renderPlot({
    
  e_precursor %>% 
      filter(precursor == pre_mir()) %>% 
      dplyr::arrange(GA) %>% 
      ggplot(aes(x = sex, y = log2(RPM + 1))) +
      geom_boxplot(aes(fill = sex), colour = "#595959", alpha = 0.9) +
      scale_fill_manual(values = c("#FEB700", "#E05D5D")) +
      my_theme +
      theme(axis.text.x = element_blank()) +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "RPM\n", fill = "Sex", x = "", title = "\nExpression by Sex\nFor both mature -3p and -5p miRNAs\n") +
      facet_grid(~mature)
  })
  
}

shinyApp(ui, server)
