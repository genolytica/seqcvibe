expressionExplorerTabPanel <- function() {
    fluidRow(column(3,
        fluidRow(column(12,
            wellPanel(
                h4("Gene settings"),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaExpressionGeneList",
                        label="Select genes",
                        choices=list(
                            "Select from list"="select",
                            "Custom list"="custom",
                            "All genes"="all"
                        )
                    ),
                    conditionalPanel(
                        condition=
                            "input.rnaExpressionGeneList=='select'",
                        selectizeInput(
                            inputId="selectExpressionGeneName",
                            label="Select genes to display",
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
                            "input.rnaExpressionGeneList=='custom'",
                        tags$textarea(
                            id="rnaExpressionCustomList",
                            class="bsv-textarea",
                            rows=5,cols=30,
                            placeholder=paste("Paste gene names ",
                                "separated by newlines",sep="")
                        )
                    )
                ))
            )
        )),
        fluidRow(column(12,
            wellPanel(
                h4("Expression settings"),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaExpressionMeasureRadio",
                        label="Select RNA-Seq expression measure",
                        choices=list(
                            "Raw counts"="raw",
                            "DESeq Normalized counts"="norm",
                            "RPKM"="rpkm",
                            "RPGM"="rpgm"
                        )
                    )
                )),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaExpressionScaleRadio",
                        label="Select expression measure scale",
                        choices=list(
                            "Natural"="natural",
                            "log2"="log2"
                        )
                    )
                )),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaExpressionAverageRadio",
                        label="Select expression measure averaging",
                        choices=list(
                            "Mean"="mean",
                            "Median"="median"
                        )
                    )
                )),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaExpressionDeviationRadio",
                        label="Select expression deviation measure",
                        choices=list(
                            "Standard deviation"="sd",
                            "Median Absolute Deviation"="mad",
                            "Interquartile Range"="IQR"
                        )
                    )
                ))
            )
        ))
    ),column(9,
        fluidRow(column(12,
            wellPanel(
                htmlOutput("rnaExpressionTables")
            )
        ))
    ))
}
