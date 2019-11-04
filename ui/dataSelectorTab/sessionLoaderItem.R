sessionLoaderTabPanel <- function() {
    fluidRow(column(9,
        fluidRow(column(12,
            wellPanel(
                h4("Application Bookmarks"),
                fluidRow(column(12,
                    div(
                        class="center",
                        fluidRow(column(6, 
                            textInput(
                                inputId = "description", 
                                label = "Bookmark description", 
                                placeholder = "Save bookmark as..."
                            )), column(6, 
                            bookmarkButton(
                                id="bookmarkBtn", 
                                label = "Add Bookmark", 
                                title = "Snapshot current status of application."
                            ))),
                        DT::dataTableOutput("urlTable"),
                        tags$style(type='text/css', "#bookmarkBtn { width:100%; margin-top: 25px;}")
                    )
                )),
                fluidRow(column(12,
                    div(
                        class="pull-left",
                        disabled(actionButton(
                            inputId="deleteBM",
                            label="Delete Bookmark",
                            title="Delete selected bookmarks",
                            icon=icon("minus")
                        ))
                    )
                ))
            )
        ))
    ))
}
