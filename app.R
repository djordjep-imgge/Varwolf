library(shiny)
library(bslib)
library(htmltools)
library(tidyverse)
library(stringr)
library(tippy)
library(reactable)
library(shinycssloaders)
library(readxl)

with_tooltip <- function(value, tooltip, ...) {
  div(style = "cursor: help", class = "tag1 header",
      tippy(value, tooltip, ...))
}

#####################################-UI-#######################################

ui <- fluidPage(
  tabsetPanel(
    id = "display",
    type = "hidden",
    tabPanel("upload_page",
      div(
        fileInput(
          "upload",
          "",
          buttonLabel =  "Upload .avcf:",
          accept = ".avcf",
          width = "30%"
        ),
        style = "display:flex; 
                align-items:center;
                justify-content:center;
                height:100vh"
      )
    ),
    tabPanel("display_page",
      fluidRow(
        column(width = 1,
          h2(textOutput("name"),
            style = "display:flex; justify-content:center; font-weight: bold;"
          )
        ),
        column(width = 2,
          fileInput("HPO", "",
            buttonLabel = "Upload HPO gene list:",
            accept = ".xlsx"
          )
        ),
        column(width = 2, br(),
               selectInput("acmg",
                            "",
                            c("All" = "all",
                              "Benign" = "Benign", 
                              "Likely benign" = "Likely Benign",
                              "VUS" = "Uncertain significance",
                              "Likely pathogenic" = "Likely pathogenic",
                              "Pathogenic" = "Pathogenic"
                            ),
                            selected = "all",
               ),
               style = "display:flex;
                  justify-content:start;",
        ),
        column(width = 2, offset = 4,
               radioButtons("filter",
                            "Show filtered variants:",
                            c("Off " = "off", "On " = "on"),
                            selected = "off",
                            inline = TRUE
               ),
               style = "display:flex;
                       justify-content:end;
                       align-items:center;
                       height:10vh;"
        ),
        column(width = 1, br(),
          div(
            actionButton("reset", "Reset"),
            style = "display:flex; 
                      justify-content:center"
          )
        )
      ),
      fluidRow(
        withSpinner(reactableOutput("table"),
                      type = 4,
                      size = 1.5,
                      color = "#1f4d7a"
        )
      )
    )
  ),
  
#################-CSS-##################

  theme = bs_theme(bg = "hsl(210, 60%, 15%)",
                   fg = "white",
                   base_font = font_google("Roboto")
          ),
  tags$style('
    div[data-value="upload_page"]{
      background-image: url(varwolf-light-logo.png);
      background-size: 16%;
      background-repeat: no-repeat;
      background-attachment: fixed;
      background-position: 50% 20%;}'
  ),
  tags$head(
    tags$style(HTML("
    .tag2 {
      display: inline-block;
      padding: 0.2rem 0.75rem;
      border-radius: 20px;
      font-weight: 800;
    }

    .status-benign {
      background: hsl(120, 90%, 20%);
      color: hsl(120, 100%, 40%);
    }

    .status-likely-benign {
      background: hsl(90, 90%, 20%);
      color: hsl(90, 100%,40%);
    }

    .status-uncertain-significance {
      background: hsl(60, 90%, 20%);
      color: hsl(60, 100%, 40%);
    }

    .status-likely-pathogenic {
      background: hsl(30, 90%, 20%);
      color: hsl(0, 100%, 40%);
    }

    .status-pathogenic {
      background: hsl(0, 90%, 20%);
      color: hsl(0, 100%, 40%);
    }"))
  ),
  tags$head(
    tags$style(HTML("
    .tag1 {
      display: inline-block;
      padding: 0.2rem 0.8rem;
      border-radius: 10px;
      font-weight: 600;
    }

    .header {
      background: hsl(210, 60%, 40%);
      color: hsl(120, 100%, 100%);
    }
    }"))
  )
)

####################################-SERVER-####################################

server <- function(input, output) {

  options(shiny.maxRequestSize = 100 * 1024^2)

  table <- reactiveVal()
  tablereset <- reactiveVal()
  avcf_u <- reactiveVal()

  avcf <- reactive(read.csv(input$upload$datapath,
                            sep = "\t",
                            colClasses = "character")
  )
  
  output$name <- renderText(tools::file_path_sans_ext(input$upload$name))
  
  observeEvent(input$upload, {
    updateTabsetPanel(inputId = "display", selected = "display_page")
    avcf <- avcf() %>% unite("Position", c("Chromosome", "Position"),
                             remove = TRUE, sep = ":")
    avcf$ACMG.pred <- str_remove(avcf$ACMG.pred, "InterVar: ")
    avcf$ACMG.pred <- str_remove(avcf$ACMG.pred, " ")
    avcf$OMIM <- str_replace_all(avcf$OMIM, "_", " ")
    avcf$OMIM <- str_replace_all(avcf$OMIM, "&", ", ")
    avcf$Orphanet <- str_replace_all(avcf$Orphanet, "_", " ")
    avcf$Orphanet <- str_replace_all(avcf$Orphanet, "&", ", ")
    avcf$clinvar_clnsig <- str_replace_all(avcf$clinvar_clnsig, "_", " ")
    avcf$clinvar_clnsig <- str_replace_all(avcf$clinvar_clnsig, "&", ", ")
    avcf$clinvar_review <- str_replace_all(avcf$clinvar_review, "_", " ")
    avcf$clinvar_review <- str_replace_all(avcf$clinvar_review, "&", ", ")
    avcf$clinvar_trait <- str_replace_all(avcf$clinvar_trait, "_", " ")
    avcf$clinvar_trait <- str_replace_all(avcf$clinvar_trait, "&", ", ")
    avcf$GWAS <- str_replace_all(avcf$GWAS, "_", " ")
    avcf$GWAS <- str_replace_all(avcf$GWAS, "&", ", ")
    avcf$Consequence <- str_replace_all(avcf$Consequence, "_", " ")
    avcf$Consequence <- str_replace_all(avcf$Consequence, "&", ", ")
    avcf$BIOTYPE <- str_replace_all(avcf$BIOTYPE, "_", " ")
    avcf$BIOTYPE <- str_replace_all(avcf$BIOTYPE, "&", ", ")
    avcf$SIFT <- str_replace_all(avcf$SIFT, "_", " ")
    avcf$SIFT <- str_replace_all(avcf$SIFT, "\\(", " \\(")
    avcf$Polyphen <- str_replace_all(avcf$Polyphen, "_", " ")
    avcf$Polyphen <- str_replace_all(avcf$Polyphen, "\\(", " \\(")
    avcf <- avcf[, c(1:2, 8, 3:7, 17:19, 9:11, 15:16, 12:14, 20:39)]
    tablereset(avcf)
    table(avcf)
  })

  observeEvent(input$HPO, {
    table.hpo <- table()[table()$Gene %in% 
                 {read_excel(input$HPO$datapath)}$GENE_SYMBOL, ]
    table(table.hpo)
  })
  
  observeEvent(input$reset, {
    tablereset <- tablereset()
    table(tablereset)
  })
  
##############-Reactable-###############
  
  output$table <- renderReactable({
    reactable(
      if (input$acmg == "all") {
        if (input$filter == "off") {table()[table()$Filter == "PASS", ]}
        else {table()}
      } else {
        if (input$filter == "off") {
          table()[table()$Filter == "PASS" & table()$ACMG.pred == input$acmg, ]}
        else {table()[table()$ACMG.pred == input$acmg, ]}
      },
      height = "90vh",
      theme = reactableTheme(
        backgroundColor = "hsl(210,60%,15%)",
        borderColor = "hsl(210, 60%, 15%)",
        inputStyle = list(backgroundColor = "hsl(210, 60%, 20%)"),
        selectStyle = list(backgroundColor = "hsl(210, 60%, 20%)"),
        style = list(color = "white")
      ),
      language = reactableLang(
        pageInfo = "{rowStart}\u2013{rowEnd} of {rows} variants"),
      showPageSizeOptions = TRUE,
      pageSizeOptions = c(10, 100, 1000),
      defaultPageSize = 100,
      searchable = TRUE,
      rowStyle = list(background = "hsl(210, 60%, 10%)"),
      defaultColDef = colDef(
        align = "center",
        minWidth = 120,
        vAlign = "center",
        header = function(value) {div(value, class = "tag1 header")},
        headerStyle = list(background = "hsl(210, 60%, 15%)", color = "white")),
      columns = list(
        Position = colDef(
          name="Position",
          sticky="left",
          style = list(backgroundColor = "hsl(210, 60%, 15%)")
        ),
        dbSNP.ID = colDef(
          name = "dbSNP ID",
          sticky="left",
          style = list(backgroundColor = "hsl(210, 60%, 15%)")
        ),
        Genotype..ref.depth..alt.depth = colDef(
          header=with_tooltip("GT: RDP, ADP", "Genotype: Ref depth, Alt depth"),
          minWidth = 145,
          sticky="left",
          style = list(backgroundColor = "hsl(210, 60%, 15%)")
        ),
        Ref = colDef(
          name="Ref",
          minWidth = 80
        ),
        Alt = colDef(
          name="Alt",
          minWidth = 80
        ),
        Quality = colDef(
          name="Qual",
          minWidth = 80
        ),
        Filter = colDef(
          name="Filter",
          cell = function(value) {
            if (value == "PASS") "\u2714\ufe0f Pass" 
            else paste("\u274c", value)
          }),
        VCF.info = colDef(
          name="VCF Info",
          show = FALSE
        ),
        Gene = colDef(
          name = "Gene",
          style = list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        OMIM = colDef(
          name = "OMIM",
          minWidth = 200,
          style = list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        Orphanet = colDef(
          name = "Orphanet",
          minWidth = 200,
          style =  list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        EXON = colDef(
          name = "Exon"
        ),
        INTRON = colDef(
          name = "Intron"
        ),
        Consequence = colDef(
          name = "Class"
        ),
        IMPACT = colDef(
          name="Impact",
          cell = function(value) {
            if (value == "LOW") "Low"
            else if (value == "MODERATE") "Moderate"
            else if (value == "MODIFIER") "Modifier"
            else if (value == "HIGH") "High"
            else ""
          },
          style = function(value) {
            if (value == "LOW") {
              color <- "#09ba09"
            } else if (value == "MODERATE") {
              color <- "#e0cd20"
            } else if (value == "MODIFIER") {
              color <- "#04b8ff"
            } else {
              color <- "#ff0404"
            }
            list(color = color, fontWeight = "bold")
          }
        ),
        BIOTYPE = colDef(
          name = "Consequence",
          minWidth = 145
        ),
        AF = colDef(
          style = function(value) {
            if (value < 0.05) {
              font <- "bold"
            } else {
              font <- "normal"
            }
            list(fontWeight = font)
          },
          show = FALSE
        ),
        EUR_AF = colDef(
          style = function(value) {
            if (value < 0.05) {
              font <- "bold"
            } else {
              font <- "normal"
            }
            list(fontWeight = font)
          },
          show = FALSE
        ),
        gnomADe_AF = colDef(
          name = "Global AF",
          style = function(value) {
            if (value < 0.05) {
              font <- "bold"
            } else {
              font <- "normal"
            }
            list(fontWeight = font)
          }
        ),
        gnomADe_NFE_AF = colDef(
          name = "Europe AF",
          minWidth = 122,
          style = function(value) {
            if (value < 0.05) {
              font <- "bold"
            } else {
              font <- "normal"
            }
            list(fontWeight = font)
          }
        ),
        PUBMED = colDef(
          show = FALSE
        ),
        CADD_Phred = colDef(
          show = FALSE
        ),
        FATHMM = colDef(
          show = FALSE
        ),
        MetaSVM = colDef(
          cell = function(value) {
            if (value == "D") "Damaging"
            else if (value == "T") "Tolerated"
            else ""
          },
          style = function(value) {
            if (value == "T") {
              color <- "#09ba09"
            } else if (value == "D") {
              color <- "#ff0404"
            } else {
              color <- "#ff0404"
            }
            list(color = color, fontWeight = "bold")
          }
        ),
        MutationTaster = colDef(
          show = FALSE
        ),
        PROVEAN = colDef(
          show = FALSE
        ),
        clinvar_clnsig = colDef(
          name = "ClinVar",
          style = list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        clinvar_hgvs = colDef(
          name = "CV Name",
          style = list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        clinvar_review = colDef(
          name = "CV Submitters",
          minWidth = 160,
          style = list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        clinvar_trait = colDef(
          name = "CV Phenotype",
          minWidth = 160,
          style = list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        clinvar_var_source = colDef(
          show = FALSE
        ),
        GWAS = colDef(
          minWidth = 150,
          style =  list(backgroundColor = "hsl(210, 60%, 12.5%)")
        ),
        ACMG.pred = colDef(
          name = "ACMG",
          sticky="right",
          style = list(backgroundColor = "hsl(210, 60%, 15%)"),
          cell = function(value) {
            class <- paste0("tag2 status-",
                            tolower(str_replace_all(value, " ", "-")))
            div(value, class = class)
          }
        ),
        Scores = colDef(
          show = FALSE
        )
      ) 
    )
  })
}

shinyApp(ui = ui, server = server)