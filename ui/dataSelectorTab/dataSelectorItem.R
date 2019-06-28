dataSelectorTabPanel <- function() {
    fluidRow(column(6,
        fluidRow(column(12,
            wellPanel(
                h4("Select input data"),
                fluidRow(column(5,
                    htmlOutput("dataSource")
                ),column(5,
                    htmlOutput("dataDataset")
                ),column(2,
                    htmlOutput("dataGenome")
                )),
                fluidRow(column(12,
                    htmlOutput("dataSelectHint")
                )),
                fluidRow(column(12,
                    htmlOutput("dataCustomSamples")
                )),
                fluidRow(column(6,
                    div(
                        class="pull-left",
                        actionButton(
                            inputId="loadDataset",
                            label="Load selected data",
                            icon=icon("truck")
                        )
                    )
                ),column(3,
                    div(
                        class="pull-right",
                        actionButton(
                            inputId="clearDataset",
                            label="Clear dataset",
                            icon=icon("trash")
                        )
                    )
                ),column(3,
                    div(
                        class="pull-right",
                        disabled(actionButton(
                            inputId="createDataset",
                            label="Create dataset",
                            icon=icon("table"),
                            class="btn-sample-select"
                        ))
                    )
                ))
            )
        ))
    ),column(6,
        fluidRow(column(12,
            wellPanel(
                h4("Current dataset"),
                fluidRow(column(12,
                    div(
                        class="small table-contaner",
                        DT::dataTableOutput("currentDatasetTable")
                    )
                ))
            )
        )),
        fluidRow(column(12,
            wellPanel(
                h4("Messages"),
                fluidRow(column(12,
                    div(
                        class="message-box",
                        htmlOutput("dataSelectorMessages")
                    )
                ))
            )
        ))
    ))
}
