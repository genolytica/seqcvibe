geneExpressionTabPanel<- function() {
    fluidRow(column(12,
        wellPanel(
            h4("Coming soon...")
         # h4("Boxplot settings"),
         #        id="expressionSettings",
         #            fluidRow(br()),
         #            fluidRow(column(12,
         #                h4("Gene settings"),
         #                fluidRow(column(12,
         #                    radioButtons(
         #                        inputId="ExprGeneList",
         #                        label="Select genes for clustering",
         #                        choices=list(
         #                            "Selected genes/custom regions"="select",
         #                            "Custom list"="custom",
         #                            "Differentially expressed genes"="degenes"
         #                        )
         #                    ),
         #                    conditionalPanel(
         #                        condition=
         #                            "input.ExprGeneList=='select'",
         #                        selectizeInput(
         #                            inputId="selectExprGeneName",
         #                            label="Select genes to cluster",
         #                            multiple=TRUE,
         #                            choices=NULL,
         #                            options=list(
         #                                placeholder="Type some gene names",
         #                                selectOnTab=TRUE
         #                            )
         #                        )
         #                    ),
         #                    conditionalPanel(
         #                        condition=
         #                            "input.ExprGeneList=='custom'",
         #                        tags$textarea(
         #                            id="rnaExprCustomList",
         #                            class="bsv-textarea",
         #                            rows=5,cols=30,
         #                            placeholder=paste("Paste gene names ",
         #                                "separated by newlines",sep="")
         #                        )
         #                    ),
         #                    conditionalPanel(
         #                        condition=
         #                            "input.ExprGeneList=='degenes'",
         #                        div(class="small",
         #                            helpText(paste("Genes currently present ",
         #                                "in the differentially expressed ",
         #                                "genes table in the 'Differential ",
         #                                "expression' panel will be clustered. ",
         #                                "Please be careful with the number of ", 
         #                                "genes to be clustered (e.g. ",
         #                                "p-value=1). Very long computations ",
         #                                "will be interrupted for application ",
         #                                "safety.",sep=""))
         #                        )
         #                    )
         #                ))
         #            )),
         #            fluidRow(column(12,
         #                h4("Expression measurement"),
         #                radioButtons(
         #                    inputId="geneExprMeasureRadio",
         #                    label="Select expression measure",
         #                    choices=list(
         #                        "Normalized counts"="counts",
         #                        "RPKM"="rpkm",
         #                        "RPGM"="rpgm"
         #                    )
         #                ),
         #                radioButtons(
         #                    inputId="geneExprScaleRadio",
         #                    label="Select expression scale",
         #                    choices=list(
         #                        "Natural"="natural",
         #                        "log2"="log2"
         #                    )
         #                )
         #            ))
        )
    ))
}