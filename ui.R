# ui.R

library(auth0)
library(DT)
library(shinyjs)
library(shinyBS)

source("config/init_server_globals.R")
source("ui/dataSelectorTab/dataSelectorItem.R")
source("ui/dataSelectorTab/sessionLoaderItem.R")
source("ui/signalViewerTab/signalViewerItemGeneSignal.R")
source("ui/signalViewerTab/signalViewerItemAreaSignal.R")
source("ui/expressionViewerTab/expressionViewerItemKnownGenes.R")
source("ui/expressionViewerTab/expressionViewerItemCalculator.R")
source("ui/analysisTab/analysisItemDiffExpr.R")
source("ui/analysisTab/analysisItemMdsPca.R")
source("ui/analysisTab/analysisItemCorrelation.R")
source("ui/analysisTab/analysisItemClustering.R")
source("ui/analysisTab/analysisItemGoPathway.R")
source("ui/genomeBrowserTab/genomeBrowserItem.R")
source("ui/helpTab/helpItemDoc.R")
source("ui/helpTab/helpItemFaq.R")
source("ui/helpTab/helpItemAbout.R")

#auth0_ui(
function(request) {
    assign("request",request,envir=.GlobalEnv)
    fluidPage(
    shinyjs::useShinyjs(),
    includeCSS("www/seqcvibe.css"),
    tags$head(
        tags$link(
            rel="stylesheet",
            type="text/css",
            href="seqcvibe.css"
        ),
        tags$link(
            rel="stylesheet",
            type="text/css",
            href="pace.css"
        ),
        tags$script(
            src="pace.js"
        ),
        tags$script(HTML(
            "$(function() {
                setTimeout(function() {
                var vals = [0,1e-16,1e-8,1e-4,0.001,0.005,0.01,0.02,0.03,0.04,
                    0.05,0.1,0.2,0.3,0.4,0.5,1];
                $('#pvalue').data('ionRangeSlider').update({
                    'values': vals,
                    'from': 10
                });
                $('#fdr').data('ionRangeSlider').update({
                    'values': vals,
                    'from': 10
                });
                var fNatVals = [0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,
                    0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,
                    6,7,8,9,10]
                $('#fcNatural').data('ionRangeSlider').update({
                    'values': fNatVals,
                    'from': 10,
                    'to': 20
                });
                // Logout button functionality
                $('[data-value=Logout]').on('click',function() {
                    Shiny.setInputValue('preLogoutHelper','42');
                });
            },5)})"
        )),
        tags$script(
            'Shiny.addCustomMessageHandler("clearUrl",',
            'function(msg) {',
            '    setTimeout(function() {',
            '        history.pushState({}, "Page Title", msg.path);',
            '     },msg.delay);',
            '});'
        )
    ),
    #conditionalPanel(
    #   condition="input.universeLoading",
    #   div(
    #       class="fullscreen",div(
    #           class="splash","Starting SeqCVIBE..."
    #       )
    #   )
    #),
    navbarPage(
        id="seqcnavbar",
        title="SeqCVIBE",
        navbarMenu("Data management",icon=icon("database"),
            tabPanel("Data selector",icon=icon("pencil-square-o"),
                dataSelectorTabPanel()
            ),
            tabPanel("Sessions",icon=icon("television"),
                sessionLoaderTabPanel()
            )
        ),
        navbarMenu("Signal viewer",icon=icon("signal"),
            tabPanel("Gene signal",icon=icon("bar-chart"),
                geneSignalTabPanel()
            ),
            tabPanel("Area signal",icon=icon("area-chart"),
                areaSignalTabPanel()
            )
        ),
        navbarMenu("Expression viewer",icon=icon("eye"),
            tabPanel("Known genes",icon=icon("table"),
                expressionExplorerTabPanel()
            ),
            tabPanel("Calculator",icon=icon("calculator"),
                expressionCalculatorTabPanel()
            )
        ),
        navbarMenu("Analysis",icon=icon("flask"),
            tabPanel("Differential expression",icon=icon("star-half-o"),
                differentialExpressionTabPanel()
            ),
            tabPanel("Clustering analysis",icon=icon("sitemap"),
                clusteringTabPanel()
            ),
            tabPanel("Correlation analysis",icon=icon("line-chart"),
                correlationTabPanel()
            ),
            tabPanel("MDS/PCA analysis",icon=icon("cube"),
                mdsPcaTabPanel()
            )#,
            #tabPanel("GO/Pathway analysis",
            #    goPathwayTabPanel()
            #)
        ),
        tabPanel("Genome browser",icon=icon("map-o"),
            genomeBrowserTabPanel()
        ),
        navbarMenu("Help",icon=icon("question"),
            tabPanel("Documentation",icon=icon("book"),
                docTabPanel()
            ),
            tabPanel("FAQ",icon=icon("question-circle"),
                faqTabPanel()
            ),
            tabPanel("About",icon=icon("user"),
                #fluidRow(column(12,includeHTML("www/about.html")))
                aboutTabPanel()
            )
        ),
        tabPanel("Logout",icon=icon("sign-out"),
            hidden(textInput("preLogoutHelper",label=""))
        )
    ),
    tags$div(id="spinnerContainer",
        tags$div(class="spinner",
            tags$div(class="rect1"),
            tags$div(class="rect2"),
            tags$div(class="rect3"),
            tags$div(class="rect4"),
            tags$div(class="rect5"),
        )
    )
)}
#,info=auth0_info("config/_auth0.yml"))
#,info=a0_info)
