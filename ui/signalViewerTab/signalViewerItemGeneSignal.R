geneSignalTabPanel <- function() {
    fluidRow(column(3,
        fluidRow(column(12,
            wellPanel(
                h4("Plot data"),
                tabsetPanel(
                    id="geneType",
                    tabPanel(
                        title="Genes",
                        fluidRow(column(12,
                            selectizeInput(
                                inputId="geneGeneName",
                                label="",
                                multiple=TRUE,
                                choices=NULL,
                                options=list(
                                    placeholder="Select genes",
                                    selectOnTab=TRUE#,
                                    #create=FALSE#,
                                    #valueField="symbol",
                                    #labelField="symbol",
                                    #searchField="symbol",
                                    #load=I("function(query,callback) {
                                    #   if (!query.length) return callback();
                                    #   $.ajax({
                                    #       url: 'http://rest.genenames.org/search/symbol/'+encodeURIComponent(query+'*'),
                                    #       type: 'GET',
                                    #       dataType: 'json',
                                    #       error: function() {
                                    #           callback();
                                    #       },
                                    #       success: function(res) {
                                    #           callback(res.response.docs);
                                    #       }
                                    #   });
                                    #}")
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(class="small table-container",
                                DT::dataTableOutput("knownGeneList")
                            )
                        ))
                    ),
                    tabPanel(
                        title="Custom",
                        fluidRow(column(8,
                            textInput(
                                inputId="customName", 
                                label="", 
                                value="",
                                placeholder="Name"
                            )
                        ),column(4,
                            selectizeInput(
                                inputId="customStrand",
                                label="",
                                choices=c("","+","-"),
                                options=list(placeholder="Strand")
                            )
                        )),
                        fluidRow(column(4,
                            htmlOutput("setChrs")
                        ),column(4,
                            textInput(
                                inputId="customStart", 
                                label="", 
                                value="",
                                placeholder="Start"
                            )
                        ),column(4,
                            textInput(
                                inputId="customEnd",
                                label="",
                                value="",
                                placeholder="End"
                            )
                        )),
                        fluidRow(column(12,
                            div(class="small table-container",
                                DT::dataTableOutput("customRegionList")
                            )
                        )),
                        fluidRow(column(12,
                            htmlOutput("customRegionError")
                        )),
                        fluidRow(br()),
                        fluidRow(column(4,""
                        ),column(4,
                            div(
                                class="pull-right",
                                style="display:inline-block",
                                actionButton("removeCustomRegion",
                                    "Remove",icon=icon("minus"))
                            )
                        ),column(4,
                            div(
                                class="pull-right",
                                style="display:inline-block",
                                actionButton("addCustomRegion",
                                    "Add",icon=icon("plus"))
                            )
                        ))
                    )
                )
            )
        )),
        fluidRow(column(12,
            wellPanel(
                h4("Plot options"),
                h5("Flanking",style="font-size:1.2em;margin-top:20px;"),
                fluidRow(column(6,
                    textInput(
                        inputId="upstreamFlank", 
                        label="Upstream", 
                        value=2000
                    )
                ),
                column(6,
                    textInput(
                        inputId="downstreamFlank",
                        label="Downstream",
                        value=2000
                    )
                )),
                fluidRow(column(12,
                    htmlOutput("regionFlankError")
                )),
                fluidRow(br()),
                fluidRow(column(12,
                    htmlOutput("geneExplorerColours")
                )),
                fluidRow(br()),
                fluidRow(column(8,
                    radioButtons(
                        inputId="geneSumStatType",
                        label="Gene profile averaging",
                        choices=list(
                            "Mean"="mean",
                            "Median"="median",
                            "Trimmed mean"="trimmed"
                        )
                    )
                ),column(4,
                    textInput(
                        inputId="geneTrimPct", 
                        label="Trim fraction", 
                        value="0.1"
                    )
                )),
                fluidRow(br()),
                fluidRow(column(8,
                    htmlOutput("geneExplorerError")
                ),column(4,
                     div(
                         class="pull-right",
                         style="display:inline-block",
                         actionButton(
                            inputId="createGeneProfile",
                            label="Engage!",
                            icon=icon("rocket")
                        )
                     )
                ))
            )
        ))
    ),column(9,
        fluidRow(column(12,
            plotOutput("geneProfile",height="640px")
        )),
        fluidRow(column(6,
            div("")
        ),column(2,
            downloadButton(
                outputId="exportGeneGG2",
                label="Export ggplot2",
                #icon=icon("file-image-o"),
                class="pull-right"
            )
        ),column(2,
            downloadButton(
                outputId="exportGenePNG",
                label="Export PNG",
                #icon=icon("file-image-o"),
                class="pull-right"
            )
        ),column(2,
            downloadButton(
                outputId="exportGenePDF",
                label="Export PDF",
                #icon=icon("file-pdf-o"),
                class="pull-right"
            )
        )),
        fluidRow(br()),
        fluidRow(column(12,
            wellPanel(
                h4("Messages"),
                div(
                    class="message-box",
                    htmlOutput("geneExplorerMessages")
                )
            )
        ))
    ))
}
