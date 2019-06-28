expressionCalculatorTabPanel <- function() {
    fluidRow(column(3,
        fluidRow(column(12,
            wellPanel(
                h4("Region settings"),
                radioButtons(
                    inputId="customExpressionFromGeneExplorer",
                    label="Custom region settings",
                    choices=list(
                        "From gene explorer"="fromge",
                        "Custom regions"="custom"
                    )
                ),
                conditionalPanel(condition=
                    "input.customExpressionFromGeneExplorer=='custom'",
                    fluidRow(column(8,
                        textInput(
                            inputId="customNameExpr", 
                            label="", 
                            value="",
                            placeholder="Name"
                        )
                    ),column(4,
                        selectizeInput(
                            inputId="customStrandExpr",
                            label="",
                            choices=c("","+","-","*"),
                            options=list(placeholder="Strand")
                        )
                    )),
                    fluidRow(column(4,
                        htmlOutput("setChrsExpr")
                    ),column(4,
                        textInput(
                            inputId="customStartExpr", 
                            label="", 
                            value="",
                            placeholder="Start"
                        )
                    ),column(4,
                        textInput(
                            inputId="customEndExpr",
                            label="",
                            value="",
                            placeholder="End"
                        )
                    )),
                    fluidRow(column(12,
                        div(class="small table-container",
                            DT::dataTableOutput(
                                "customRegionListExpr")
                        )
                    )),
                    fluidRow(column(12,
                        htmlOutput("customRegionExprError")
                    )),
                    fluidRow(br()),
                    fluidRow(column(4,""
                    ),column(4,
                        div(
                            class="pull-right",
                            style="display:inline-block",
                            actionButton("removeRnaCustomRegion",
                                "Remove",icon=icon("minus"))
                        )
                    ),column(4,
                        div(
                            class="pull-right",
                            style="display:inline-block",
                            actionButton("addRnaCustomRegion",
                                "Add",icon=icon("plus"))
                        )
                    ))
                ),
                fluidRow(br()),
                fluidRow(column(8,
                    htmlOutput("customRnaCalcError")
                ),column(4,
                    div(
                        class="pull-right",
                        style="display:inline-block",
                        actionButton(
                            inputId="calculateCustomRegionRna",
                            label="Engage!",
                            icon=icon("rocket")
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
                        inputId="rnaCustomMeasureRadio",
                        label="Select RNA-Seq expression measure",
                        choices=list(
                            "Raw counts"="raw",
                            #"DESeq Normalized counts"="norm",
                            "RPKM"="rpkm",
                            "RPGM"="rpgm"
                        )
                    )
                )),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaCustomScaleRadio",
                        label="Select expression measure scale",
                        choices=list(
                            "Natural"="natural",
                            "log2"="log2"
                        )
                    )
                )),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaCustomAverageRadio",
                        label="Select expression measure averaging",
                        choices=list(
                            "Mean"="mean",
                            "Median"="median"
                        )
                    )
                )),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaCustomDeviationRadio",
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
                htmlOutput("rnaCustomTables")
            )
        ))
    ))
}
