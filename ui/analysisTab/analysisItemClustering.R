clusteringTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            h4("Clustering settings"),
            tabsetPanel(
                id="clusteringSettings",
                tabPanel(
                    title="Genes",
                    fluidRow(br()),
                    fluidRow(column(12,
                        h4("Gene settings"),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaClusterGeneList",
                                label="Select genes for clustering",
                                choices=list(
                                    "Selected genes/custom regions"="select",
                                    "Custom list"="custom",
                                    "Differentially expressed genes"="degenes"
                                )
                            ),
                            conditionalPanel(
                                condition=
                                    "input.rnaClusterGeneList=='select'",
                                selectizeInput(
                                    inputId="selectClusteringGeneName",
                                    label="Select genes to cluster",
                                    multiple=TRUE,
                                    choices=NULL,
                                    options=list(
                                        placeholder="Type some gene names",
                                        selectOnTab=TRUE
                                    )
                                )
                            ),
                            conditionalPanel(
                                condition=
                                    "input.rnaClusterGeneList=='custom'",
                                tags$textarea(
                                    id="rnaClusteringCustomList",
                                    class="bsv-textarea",
                                    rows=5,cols=30,
                                    placeholder=paste("Paste gene names ",
                                        "separated by newlines",sep="")
                                )
                            ),
                            conditionalPanel(
                                condition=
                                    "input.rnaClusterGeneList=='degenes'",
                                div(class="small",
                                    helpText(paste("Genes currently present ",
                                        "in the differentially expressed ",
                                        "genes table in the 'Differential ",
                                        "expression' panel will be clustered. ",
                                        "Please be careful with the number of ", 
                                        "genes to be clustered (e.g. ",
                                        "p-value=1). Very long computations ",
                                        "will be interrupted for application ",
                                        "safety.",sep=""))
                                )
                            )
                        ))
                    )),
                    fluidRow(column(12,
                        h4("Clustering variables"),
                        radioButtons(
                            inputId="rnaClusteringVariablesRadio",
                            label="Select clustering variable",
                            choices=list(
                                "Individual replicates"="replicates",
                                "Condition averages"="averages",
                                "Fold changes"="fc"
                            )
                        ),
                        conditionalPanel(
                            condition="input.rnaClusteringVariablesRadio=='fc'",
                            htmlOutput("rnaClusteringFcControl")
                        )
                    )),
                    fluidRow(column(12,
                        h4("Expression measurement"),
                        radioButtons(
                            inputId="rnaClusteringMeasureRadio",
                            label="Select expression measure",
                            choices=list(
                                "Normalized counts"="counts",
                                "RPKM"="rpkm",
                                "RPGM"="rpgm"
                            )
                        ),
                        radioButtons(
                            inputId="rnaClusteringScaleRadio",
                            label="Select expression scale",
                            choices=list(
                                "Natural"="natural",
                                "log2"="log2"
                            )
                        )
                    ))
                ),
                tabPanel(
                    title="Heatmap",
                    fluidRow(br()),
                    fluidRow(column(12,
                        h4("Heatmap settings"),
                        selectizeInput(
                            inputId="selectClusteringDistance",
                            label="Select distance metric",
                            choices=c(
                                "Euclidean"="euclidean",
                                "Maximum"="maximum",
                                "Manhattan"="manhattan",
                                #"Canberra"="canberra",
                                "Minkowski"="minkowski",
                                "Pearson correlation"="pearson",
                                "Spearman correlation"="spearman",
                                "Cosine"="cosine"
                            )
                        ),
                        selectizeInput(
                            inputId="selectClusteringLinkage",
                            label="Select linkage function (dendrogram)",
                            choices=c(
                                "Average"="average",
                                "Complete"="complete",
                                "Single"="single",
                                "McQuitty"="mcquitty",
                                "Median"="median",
                                "Centroid"="centroid",
                                "Ward v1"="ward1",
                                "Ward v2"="ward2"
                            )
                        ),
                        radioButtons(
                            inputId="rnaClusterWhat",
                            label="Select clustering dimensions",
                            choices=list(
                                "Cluster both rows and columns"="both",
                                "Cluster rows"="row",
                                "Cluster columns"="column",
                                "No clustering (only show image)"="none"
                            )
                        )
                    )),
                    fluidRow(column(6,
                        textInput(
                            inputId="kRowInput",
                            label="# of row clusters",
                            value="1"
                        )
                    ),column(6,
                        textInput(
                            inputId="kColInput",
                            label="# of column clusters",
                            value="1"
                        )
                    )),
                    fluidRow(column(12,
                        colourInput(
                            inputId="heatmapColourExtremeDown",
                            label="Colour for extreme low values",
                            value="#1A9641"
                        ),
                        colourInput(
                            inputId="heatmapColourMildDown",
                            label="Colour for mild low values",
                            value="#A6D96A"
                        ),
                        colourInput(
                            inputId="heatmapColourMiddle",
                            label="Colour for neutral values",
                            value="#FFFFBF"
                        ),
                        colourInput(
                            inputId="heatmapColourMildUp",
                            label="Colour for mild high values",
                            value="#FDAE61"
                        ),
                        colourInput(
                            inputId="heatmapColourExtremeUp",
                            label="Colour for extreme high values",
                            value="#D7191C"
                        )
                    )),
                    fluidRow(column(12,
                        checkboxInput(
                            inputId="checkColorSaturation",
                            label="Saturate colors (outlier protection)",
                            value=FALSE
                        ),
                        conditionalPanel(condition="input.checkColorSaturation",
                            div(class="small",
                                helpText(paste("Use this option when there is ",
                                    "a sufficient number of genes to cluster ",
                                    "(e.g. >100) so as to compute meaningful ",
                                    "saturation values.",sep=""))
                            ),
                            fluidRow(column(6,
                                textInput(
                                    inputId="colorSaturationLowQuantile",
                                    label="Low quantile",
                                    value="0.05"
                                )
                            ),column(6,
                                textInput(
                                    inputId="colorSaturationHighQuantile",
                                    label="High quantile",
                                    value="0.95"
                                )
                            ))
                        )
                    ))
                )
            ),
            fluidRow(br()),
            fluidRow(column(8,
                htmlOutput("rnaClusteringSettingsError")
            ),column(4,
                 div(
                     class="pull-right",
                     style="display:inline-block",
                     actionButton(
                        inputId="performRnaClustering",
                        label="Engage!",
                        icon=icon("rocket")
                    )
                 )
            ))
        )
    ),column(9,
        fluidRow(column(12,
            htmlOutput("heatmapOutput")
        )),
        fluidRow(br()),
        fluidRow(br()),
        fluidRow(column(8,
            div("")
        ),column(2,
            downloadButton(
                outputId="exportRnaHeatmapPNG",
                label="Export PNG",
                #icon=icon("file-image-o"),
                class="pull-right"
            )
        ),column(2,
            downloadButton(
                outputId="exportRnaHeatmapPDF",
                label="Export PDF",
                #icon=icon("file-pdf-o"),
                class="pull-right"
            )
        ))
    ))
}

