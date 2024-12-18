library(shiny)
library(DT)

shinyUI(navbarPage(title="Multicenter NGS studies between theory and practice",
#                   theme=shinytheme("slate"),
                   tabPanel("Minimum sequencing requirements",
                   fluidRow(
                       column(4,wellPanel(
                         h3("Target region:"),
                         radioButtons('targetRegion','',choices = c("Region of Interest","Intersecting Region"),selected = "Region of Interest",inline = T),
                         hr(),
                         
                         h3("Evaluation Method:"),
                         radioButtons('evaluationMethod','',choices = c("Absolute","Relative"),selected = "Absolute",inline = T),
                         hr(),
                         
                         radioButtons('cov','Coverage',choices = c("Customize","Default"),selected = "Default",inline = T),
                         uiOutput("covUI1"),
                         uiOutput("covUI1r"),
                         hr(),
                           
                         radioButtons('mq','Mapping Quality',choices = c("Customize","Default"),selected = "Default",inline = T),
                         uiOutput("mqUI1"),
                         uiOutput("mqUI1r"),
                         
                         radioButtons('bq','Base Quality',choices = c("Customize","Default"),selected = "Default",inline = T),
                         uiOutput("bqUI1"),
                         uiOutput("bqUI1r"),
                         
                         radioButtons('bn','Background Noise',choices = c("Customize","Default"),selected = "Default",inline = T),
                         uiOutput("bnUI1"),
                         uiOutput("bnUI1r"),
                         hr(),
                           
                           radioButtons('weight','Weighting',choices = c("Customize","Default"),selected = "Default",inline = T),
                           uiOutput("weightUI1"),
                           fluidRow(
                               column(4,uiOutput("weightUI2")),
                               column(4,uiOutput("weightUI3")),
                               column(4,uiOutput("weightUI4"))
                           ),
                           hr(),
                           
                           radioButtons('color','Color',choices = c("Customize","Default"),selected = "Default",inline = T),
                           uiOutput("colorUI1"),
                           uiOutput("colorUI2"),
                           hr(),
                           
                           radioButtons('datasets','Data sets',choices = c("Customize","Default"),selected = "Default",inline = T),
                           fluidRow(
                             column(4,uiOutput("datasetUI1")),
                             column(4,uiOutput("datasetUI2")),
                             column(4,uiOutput("datasetUI3"))
                           ),
                           fluidRow(
                             column(4,uiOutput("datasetUI4")),
                             column(4,uiOutput("datasetUI5")),
                             column(4,uiOutput("datasetUI6"))
                           ),
                           fluidRow(
                             column(4,uiOutput("datasetUI7")),
                             column(4,uiOutput("datasetUI8")),
                             column(4,uiOutput("datasetUI9"))
                           ),
                           fluidRow(
                             column(4,uiOutput("datasetUI10")),
                             column(4,uiOutput("datasetUI11"))
                           ),
                           hr(),
                           
                           actionButton("do1", "Update Analysis")
                       )),
                       column(8,tabsetPanel(
                           tabPanel("Plot",plotOutput("plot")),
                           tabPanel("Table",DT::dataTableOutput('table'))
                       )))),
                   
                   tabPanel("Mutation frequency",
                            fluidRow(
                                column(4,wellPanel(
                                    radioButtons('change_impact','Mutation frequency in MDS (select \'0\' to ignore a gene)',choices = c("Cosmic","Haferlach et al. 2014","Papaemmanuil et al. 2013","Customize"),selected = "Cosmic",inline = T),
                           fluidRow(
                               column(4,uiOutput("asxl1UI")),
                               column(4,uiOutput("bcorUI")),
                               column(4,uiOutput("cblUI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("cebpaUI")),
                               column(4,uiOutput("dnmt3aUI")),
                               column(4,uiOutput("etv6UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("ezh2UI")),
                               column(4,uiOutput("flt3UI")),
                               column(4,uiOutput("gata2UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("idh1UI")),
                               column(4,uiOutput("idh2UI")),
                               column(4,uiOutput("jak2UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("kdm6aUI")),
                               column(4,uiOutput("kitUI")),
                               column(4,uiOutput("kmt2aUI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("krasUI")),
                               column(4,uiOutput("mplUI")),
                               column(4,uiOutput("npm1UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("nrasUI")),
                               column(4,uiOutput("rad21UI")),
                               column(4,uiOutput("runx1UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("sf3b1UI")),
                               column(4,uiOutput("smc1aUI")),
                               column(4,uiOutput("smc3UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("srsf2UI")),
                               column(4,uiOutput("stag2UI")),
                               column(4,uiOutput("tet2UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("tp53UI")),
                               column(4,uiOutput("u2af1UI")),
                               column(4,uiOutput("wt1UI"))
                           ),
                           fluidRow(
                               column(4,uiOutput("zrsr2UI"))
                           ),
                           hr(),
                           hr(),
                           actionButton("do2", "Update Analysis")
                                )),
                           column(8,tabsetPanel(
                               tabPanel("Plot",plotOutput("plot1")),
                               tabPanel("Table",DT::dataTableOutput('table1'))
                           ))))
                   

)
)