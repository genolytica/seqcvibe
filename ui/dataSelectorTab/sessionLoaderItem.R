sessionLoaderTabPanel <- function() {
    fluidRow(column(9,
        fluidRow(column(12,
            wellPanel(
                h4("Session manager"),
                fluidRow(column(12,
                    div(
                        class="center",
                        fluidRow(column(6, 
                            textInput(
                                inputId="description", 
                                label="Session description", 
                                placeholder="Session name"
                            )), column(6, 
                            disabled(bookmarkButton(
                                id="bookmarkBtn", 
                                label = "Add session", 
                                title = "Snapshot current status of application."
                            ))
                        )),
                        DT::dataTableOutput("urlTable"),
                        tags$style(type='text/css', "#bookmarkBtn { width:100%; margin-top: 25px;}")
                    )
                ))#,
                #fluidRow(column(12,
                #    div(
                #        class="pull-left",
                #        disabled(actionButton(
                #            inputId="deleteBM",
                #            label="Delete session",
                #            title="Delete selected sessions",
                #            icon=icon("minus")
                #        ))
                #    )
                #))
            )
        ))
    ))
}
