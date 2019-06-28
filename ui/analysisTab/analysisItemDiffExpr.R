differentialExpressionTabPanel <- function() {
    fluidRow(column(3,
        fluidRow(column(12,
            wellPanel(
                h4("Pipeline settings"),
                tabsetPanel(
                    id="pipSettingsTabset",
                    tabPanel(
                        title="General",
                        fluidRow(br()),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDePipeline",
                                label="Select pipeline",
                                choices=list(
                                    "DESeq"="deseq",
                                    "edgeR"="edger",
                                    "voom"="limma"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            htmlOutput("rnaDePipelineControl")
                        )),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeNormalizeWhen",
                                label="Apply filters when",
                                choices=list(
                                    "Before normalization"="prenorm",
                                    "After normalization"="postnorm"
                                ),
                                selected="postnorm"
                            )
                        )),
                        fluidRow(column(12,
                            selectInput(
                                inputId="rnaDeMTC",
                                label="Multiple testing correction",
                                choices=list(
                                    "False Discovery Rate (FDR)"=c(
                                        "Benjamini-Hochberg"="BH",
                                        "Benjamini-Yekutieli"="BY"
                                    ),
                                    "Family Wise Error Rate (FWER)"=c(
                                        "Holm"="holm",
                                        "Hommel"="hommel",
                                        "Bonferroni"="bonferroni"
                                    )
                                )
                            )
                        )),
                        fluidRow(br()),
                        div(style="font-weight:bold","Custom regions"),
                        fluidRow(column(12,
                            checkboxInput(
                                inputId="includeCustomRegions",
                                label=paste("Include regions from expression ",
                                    "calculator"),
                                value=FALSE
                            )
                        ))
                    ),
                    tabPanel(
                        title="Filtering",
                        fluidRow(br()),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeGeneFilter",
                                label="Basic gene filter",
                                choices=list(
                                    "Median expression"="median",
                                    "Mean expression"="mean",
                                    "Quantile"="quantile",
                                    "Known genes"="known"
                                )
                            ),
                            conditionalPanel(
                                condition="input.rnaDeGeneFilter=='quantile'",
                                textInput(
                                    inputId="rnaDeQuantileFilter", 
                                    label="Quantile", 
                                    value="0.25"
                                )
                            ),
                            conditionalPanel(
                                condition="input.rnaDeGeneFilter=='known'",
                                selectizeInput(
                                    inputId="rnaDeKnownFilter",
                                    label="List of genes for filtering",
                                    multiple=TRUE,
                                    choices=NULL,
                                    options=list(
                                        placeholder="Select genes",
                                        selectOnTab=TRUE
                                    )
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(style="font-weight:bold","Gene length filter"),
                            checkboxInput(
                                inputId="rnaDeGeneLengthFilter",
                                label="Gene length filtering",
                                value=FALSE
                            ),
                            conditionalPanel(
                                condition="input.rnaDeGeneLengthFilter",
                                textInput(
                                    inputId="rnaDeGeneLengthFilterValue",
                                    label="Gene length in bp",
                                    value="500"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(style="font-weight:bold","Gene biotype filters"),
                            checkboxInput(
                                inputId="rnaDeBiotypeFilter",
                                label="Biotype expression filtering",
                                value=FALSE
                            ),
                            conditionalPanel(
                                condition="input.rnaDeBiotypeFilter",
                                htmlOutput("checkboxBiotypeListRna")
                            )
                        ))
                    )
                ),
                fluidRow(br()),
                fluidRow(column(8,
                    htmlOutput("rnaDeSettingsError")
                ),column(4,
                     div(
                         class="pull-right",
                         style="display:inline-block",
                         actionButton(
                            inputId="performDeAnalysis",
                            label="Engage!",
                            icon=icon("rocket")
                        )
                     )
                ))
            ),
            wellPanel(
                h4("Results table settings"),
                tabsetPanel(
                    id="rnaDeTableSettings",
                    tabPanel(
                        title="Score",
                        fluidRow(br()),
                        fluidRow(column(12,
                            disabled(radioButtons(
                                inputId="statThresholdType",
                                label="Statistical score threshold",
                                inline=TRUE,
                                choices=c(
                                    "p-value"="pvalue",
                                    "FDR"="fdr"
                                ))
                            ),
                            conditionalPanel(
                                condition="input.statThresholdType=='pvalue'",
                                disabled(sliderInput(
                                    inputId="pvalue",
                                    label="p-value",
                                    min=0,
                                    max=1,
                                    value=0.05,
                                    step=0.01
                                ))
                            ),
                            conditionalPanel(
                                condition="input.statThresholdType=='fdr'",
                                disabled(sliderInput(
                                    inputId="fdr",
                                    label="FDR",
                                    min=0,
                                    max=1,
                                    value=0.05,
                                    step=0.01
                                ))
                            )
                        )),
                        fluidRow(column(12,
                            div(
                                style="font-weight:bold",
                                "Fold change threshold"
                            ),
                            div(
                                class="small",
                                helpText(paste("Select the fold change ",
                                    "display scale (natural or log2) from the ",
                                    "'Scale' tab above."))
                            ),
                            fluidRow(column(12,
                                disabled(selectizeInput(
                                    inputId="rnaDeCurrentContrast",
                                    label="Select contrast for fold change",
                                    choices=NULL
                                ))
                            )),
                            conditionalPanel(
                                condition=
                                    "input.rnaDeValueScaleRadio=='natural'",
                                disabled(sliderInput(
                                    inputId="fcNatural",
                                    label="Fold change (natural)",
                                    min=0,
                                    max=10,
                                    value=c(0.5,2),
                                    step=0.5
                                ))
                            ),
                            conditionalPanel(
                                condition="input.rnaDeValueScaleRadio=='log2'",
                                disabled(sliderInput(
                                    inputId="fcLog",
                                    label="Fold change (log2)",
                                    min=-5,
                                    max=5,
                                    value=c(-1,1),
                                    step=0.5
                                ))
                            )
                        ))
                    ),
                    tabPanel(
                        title="Filters",
                        fluidRow(br()),
                        fluidRow(column(12,
                            disabled(selectizeInput(
                                inputId="rnaDeShowSpecificGenes",
                                label="Show selected genes",
                                multiple=TRUE,
                                choices=NULL,
                                options=list(
                                    placeholder="Select genes",
                                    selectOnTab=TRUE
                                )
                            ))
                        )),
                        fluidRow(column(12,
                            htmlOutput("setDeChrs")
                        )),
                        fluidRow(column(12,
                            div(style=
                                "font-weight:bold","Show selected biotypes"),
                            disabled(checkboxInput(
                                inputId="rnaDeAnalyzedBiotypeFilter",
                                label="Biotype expression filtering",
                                value=FALSE
                            )),
                            conditionalPanel(
                                condition="input.rnaDeAnalyzedBiotypeFilter",
                                htmlOutput("checkboxBiotypeListAnalyzedRna")
                            )
                        ))
                    ),
                    tabPanel(
                        title="Scale",
                        fluidRow(br()),
                        fluidRow(column(12,
                            disabled(radioButtons(
                                inputId="rnaDeValueCompRadio",
                                label="Select summary value components",
                                choices=list(
                                    "Counts"="counts",
                                    "RPKM"="rpkm",
                                    "RPGM"="rpgm"
                                )
                            ))
                        )),
                        fluidRow(column(12,
                            disabled(radioButtons(
                                inputId="rnaDeValueScaleRadio",
                                label="Select summary value scale",
                                choices=list(
                                    "Natural"="natural",
                                    "log2"="log2"
                                )
                            ))
                        )),
                        fluidRow(column(12,
                            disabled(radioButtons(
                                inputId="rnaDeValueAverageRadio",
                                label="Select summary value averaging",
                                choices=list(
                                    "Mean"="mean",
                                    "Median"="median"
                                )
                            ))
                        )),
                        fluidRow(column(12,
                            disabled(radioButtons(
                                inputId="rnaDeValueDeviationRadio",
                                label="Select summary value deviation",
                                choices=list(
                                    "Standard deviation"="sd",
                                    "Median Absolute Deviation"="mad",
                                    "Interquartile Range"="IQR"
                                )
                            ))
                        ))
                    )
                )
            )
        ))
    ),column(9,
        fluidRow(column(10,
            #div(
            #    id="maplot-container",
            #    tags$img(
            #        id="loading-spinner",
            #        src="pacman.gif",
            #    )
            #)
            #plotlyOutput("rnaDeMAPlot",height="640px")
            plotOutput(
                outputId="rnaDeMAPlot",
                click="rnaDeMAPlotClick",
                dblclick="rnaDeMAPlotDblClick",
                brush=brushOpts(
                    id="rnaDeMAPlotBrush",
                    resetOnNew=TRUE
                ),
                height="600px"
            )
        ),column(2,
            wellPanel(
                h4("Plot behaviour"),
                fluidRow(column(12,
                    checkboxInput(
                        inputId="toggleRnaDeTableUpdate",
                        label="Real time table update",
                        value=FALSE
                    ),
                    checkboxInput(
                        inputId="toggleRnaDeZoom",
                        label="Toggle zoom",
                        value=FALSE
                    )               
                )),
                fluidRow(column(12,
                    disabled(actionButton(
                        inputId="resetRnaDeZoom",
                        label="Reset zoom",
                        icon=icon("home"),
                        class="btn-sm pull-right"
                    ))
                ))
            ),
            wellPanel(
                fluidRow(column(12,
                    h4("Plot colors"),
                    htmlOutput("maPlotColours")
                ))
            )
        )),
        fluidRow(column(4,
            div("")
        ),column(2,
            downloadButton(
				outputId="exportRnaDeMAPlotGG2",
				label="Export ggplot2",
				#icon=icon("file-image-o"),
				class="pull-right"
			)
        ),column(2,
            downloadButton(
				outputId="exportRnaDeMAPlotPNG",
				label="Export PNG",
				#icon=icon("file-image-o"),
				class="pull-right"
			)
        ),column(2,
			downloadButton(
				outputId="exportRnaDeMAPlotPDF",
				label="Export PDF",
				#icon=icon("file-pdf-o"),
				class="pull-right"
			)
        ),column(2,
            div("")
        )),
        fluidRow(br()),
        fluidRow(column(12,
            wellPanel(
                tabsetPanel(
                    id="rnaDeAnalysisResults",
                    tabPanel(
                        title="Summary",
                        fluidRow(br()),
                        fluidRow(column(12,
                            htmlOutput("rnaDeAnalysisSummary")
                        ))
                    ),
                    tabPanel(
                        title="Annotation",
                        fluidRow(br()),
                        fluidRow(column(12,
                            htmlOutput("rnaDeAnalysisAnnotation")
                        ))
                    ),
                    tabPanel(
                        title="Flags",
                        fluidRow(br()),
                        fluidRow(column(12,
                            htmlOutput("rnaDeAnalysisFlags")
                        ))
                    ),
                    tabPanel(
                        title="All",
                        fluidRow(br()),
                        fluidRow(column(12,
                            htmlOutput("rnaDeAnalysisAll")
                        ))
                    )
                )
            )
        ))
    ))
}
