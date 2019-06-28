expressionCalculatorTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    customRegions <- allReactiveVars$customRegions
    customRnaRegions <- allReactiveVars$customRnaRegions
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    
    addTheCustomRnaRegion <- eventReactive(input$addRnaCustomRegion,{
        regionValidate <- c(
            name=FALSE,
            strand=FALSE,
            chrom=FALSE,
            start=FALSE,
            end=FALSE,
            startBeforeEnd=FALSE,
            regionExistsName=FALSE,
            regionExistsLocus=FALSE
        )
        errMsg <- c(
            name="",
            strand="",
            chrom="",
            start="",
            end="",
            startBeforeEnd="",
            regionExistsName="",
            regionExistsLocus=""
        )
        if (input$customNameExpr=="") {
            regionValidate["name"] <- TRUE
            errMsg["name"] <- "A region name is required!"
        
        }
        if (input$customStrandExpr=="") {
            regionValidate["strand"] <- TRUE
            errMsg["strand"] <- "A region strand is required!"
        }
        if (input$customChrExpr=="") {
            regionValidate["chrom"] <- TRUE
            errMsg["chrom"] <- "A region chromosome is required!"
        
        }
        if (input$customStartExpr=="" 
            || as.integer(input$customStartExpr)<=0
            || is.na(suppressWarnings(
                as.numeric(input$customStartExpr)))) {
            regionValidate["start"] <- TRUE
            errMsg["start"] <- 
                "Region start base must be a valid positive integer!"
        }
        if (input$customEndExpr=="" 
            || as.numeric(input$customEndExpr)<=0
            || is.na(suppressWarnings(
                as.integer(input$customEndExpr)))) {
            regionValidate["end"] <- TRUE
            errMsg["end"] <- 
                "Region end base must be a valid positive integer!"
        }
        if (!(regionValidate["start"] || regionValidate["end"])) {
            # Check start before end
            s <- as.integer(input$customStartExpr)
            e <- as.integer(input$customEndExpr)
            if (e-s<=0) {
                regionValidate["startBeforeEnd"] <- TRUE
                errMsg["startBeforeEnd"] <- paste("Region start must be ",
                    "smaller than region end and the region length must ",
                    "be greater than 1!",sep="")
            }
            
            # Check the new region does not already exists
            cc <- customRnaRegions$chromosome
            ss <- customRnaRegions$start
            ee <- customRnaRegions$end
            tt <- customRnaRegions$strand
            if (input$customNameExpr %in% customRnaRegions$name) {
                regionValidate["regionExistsName"] <- TRUE
                errMsg["regionExistsName"] <- paste("Custom regions must ",
                    "have each a unique name! ",input$customNameExpr," ",
                    "already exists.",sep="")
            }
            if (length(ss)>0) {
                for (i in 1:length(ss)) {
                    if (input$customChrExpr==cc[i] 
                        && input$customStartExpr==ss[i] 
                        && input$customEndExpr==ee[i] 
                        && input$customStrandExpr==tt[i])
                    regionValidate["regionExistsLocus"] <- TRUE
                    errMsg["regionExistsLocus"] <- 
                        "Custom regions must have unique coordinates!"
                }
            }
        }
        
        if (any(regionValidate)) {
            output$customRegionExprError <- renderUI({
                div(class="error-message",paste(errMsg,collapse="\n"))
            })
            return("")
        }
        else {
            output$customRegionExprError <- renderUI({div()})
            customRnaRegions$chromosome <- c(customRnaRegions$chromosome,
                input$customChrExpr)
            customRnaRegions$start <- c(customRnaRegions$start,
                input$customStartExpr)
            customRnaRegions$end <- c(customRnaRegions$end,
                input$customEndExpr)
            customRnaRegions$strand <- c(customRnaRegions$strand,
                input$customStrandExpr)
            customRnaRegions$name <- c(customRnaRegions$name,
                input$customNameExpr)
        }
    })
    
    removeTheCustomRnaRegion <- eventReactive(
        input$removeRnaCustomRegion,{
        output$customRegionExprError <- renderUI({div()})
        j <- as.integer(input$customRegionListExpr_rows_selected)
        if (length(j)>0) {
            ii <- customRnaRegions$name[j]
            customRnaRegions$chromosome <- 
                customRnaRegions$chromosome[-j]
            customRnaRegions$start <- customRnaRegions$start[-j]
            customRnaRegions$end <- customRnaRegions$end[-j]
            customRnaRegions$strand <- customRnaRegions$strand[-j]
            customRnaRegions$name <- customRnaRegions$name[-j]
        }
    })
    
    calcCustomRnaCounts <- eventReactive(input$calculateCustomRegionRna,{
        if (is.null(currentMetadata$final)) {
            output$customRnaCalcError <- renderUI({
                div(class="error-message",paste("You must create a ",
                    "dataset first!",sep=""))
            })
        }
        else {
            output$customRnaCalcError <- renderUI({div()})
            if (input$customExpressionFromGeneExplorer=="fromge")
                customRnaRegions <- customRegions
            customInputRegions <- makeGRangesFromDataFrame(
                df=data.frame(
                    chromosome=as.character(customRnaRegions$chromosome),
                    start=as.numeric(customRnaRegions$start),
                    end=as.numeric(customRnaRegions$end),
                    strand=as.character(customRnaRegions$strand),
                    gene_id=as.character(customRnaRegions$name)
                ),
                keep.extra.columns=TRUE
            )
            currentCustomRnaTables$lengths <- width(customInputRegions)
            names(currentCustomRnaTables$lengths) <- 
                names(customInputRegions)
            #customInputRegions <- as.list(customInputRegions)
            names(customInputRegions) <- as.character(customRnaRegions$name)
            progress <- shiny::Progress$new()
            progress$initialize(
                session,
                min=0,
                max=length(customInputRegions)
            )
            progress$set(message="",value=0)
            on.exit(progress$close())
            
            updateProgress <- function(value=NULL,detail=NULL) {
                if (is.null(value)) {
                    value <- progress$getValue()
                    value <- value + 1
                }
                progress$set(value=value,detail=detail)
            }
            
            customCounts <- tryCatch({
                getCustomCounts(
                    coords=customInputRegions,
                    samples=as.character(currentMetadata$final$sample_id),
                    config=currentMetadata$final,
                    progressFun=updateProgress,
                    rc=RC
                )},
                error=function(e) {
                },finally=""
            )
            
            for (c in unique(as.character(currentMetadata$final$class))) {
                sc <- as.character(currentMetadata$final[
                    currentMetadata$final$class==c,"sample_id"])
                currentCustomRnaTables$tables[[c]] <- 
                    customCounts[,sc,drop=FALSE]
            }
        }
    })
    
    return(list(
        addTheCustomRnaRegion=addTheCustomRnaRegion,
        removeTheCustomRnaRegion=removeTheCustomRnaRegion,
        calcCustomRnaCounts=calcCustomRnaCounts
    ))
}

expressionCalculatorTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    
    createCurrentCustomRnaTables <- reactive({
        for (c in currentMetadata$class)
            currentCustomRnaTables$tables[[c]] <- NULL
    })
    
    handleCustomRnaSampleSelection <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            observeEvent(input[[paste("clearCustomRnaCountSelection_",x,
                sep="")]],{
                proxy <- dataTableProxy(paste("customCountRnaTable_",x,
                    sep=""))
                selectRows(proxy,NULL)
            })
        })
        lapply(c,function(x) {
            observeEvent(input[[paste("invertCustomRnaCountSelection_",x,
                sep="")]],{
                N <- input[[paste("customCountRnaTable_",x,"_rows_all",
                    sep="")]]
                sel <- input[[paste("customCountRnaTable_",x,
                    "_rows_selected",sep="")]]
                if (length(sel)>0) {
                    N <- N[-sel]
                    proxy <- dataTableProxy(paste("customCountRnaTable_",x,
                        sep=""))
                    selectRows(proxy,N)
                }
            })
        })
    })
    
    handleCustomRnaDownloadSelection <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            output[[paste("exportCustomRnaCountSelection_",x,sep="")]] <- 
                downloadHandler(
                    filename=function() {
                        tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                        paste(x,"_custom_expression_",tt,".txt", sep='')
                    },
                    content=function(con) {
                        sel <- input[[paste("customCountRnaTable_",x,
                            "_rows_selected",sep="")]]
                        if (length(sel)>0)
                            write.table(currentCustomRnaTables$displayTables[[
                                x]][sel,,drop=FALSE],file=con,sep="\t",
                                quote=FALSE,col.names=NA)
                    }
                )
            output[[paste("exportCustomRnaCountAll_",x,sep="")]] <- 
                downloadHandler(
                    filename=function() {
                        tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                        paste(x,"_custom_expression_",tt,".txt", sep='')
                    },
                    content=function(con) {
                        write.table(currentCustomRnaTables$displayTables[[x]],
                            file=con,sep="\t",quote=FALSE,col.names=NA)
                    }
                )
        })
    })
    
    return(list(
        createCurrentCustomRnaTables=createCurrentCustomRnaTables,
        handleCustomRnaSampleSelection=handleCustomRnaSampleSelection,
        handleCustomRnaDownloadSelection=handleCustomRnaDownloadSelection
    ))
}

expressionCalculatorTabPanelRenderUI <- function(output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    customRnaRegions <- allReactiveVars$customRnaRegions
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    
    output$customRegionListExpr <- DT::renderDataTable(
        data.frame(
            chromosome=customRnaRegions$chromosome,
            start=customRnaRegions$start,
            end=customRnaRegions$end,
            strand=customRnaRegions$strand,
            name=customRnaRegions$name
        ),
        class="display compact",
        rownames=FALSE,
        options=list(
            dom='t',
            paging=FALSE
        )
    )
    
    output$rnaCustomTables <- renderUI({
        M <- currentMetadata$final
        if (is.null(M))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please create a dataset first using the",
                    "'Data selector' tab."))
            )
        else {
            M <- currentMetadata$final
            s <- unique(as.character(M$source))
            d <- unique(as.character(M$dataset))
            c <- unique(as.character(M$class))
            le <- currentCustomRnaTables$lengths
            lapply(c,function(x,M,le,D) {
                tab <- currentCustomRnaTables$tables[[x]]
                if (is.null(tab)) {
                    tab <- data.frame(
                        name=character(0),
                        average=numeric(0),
                        deviation=numeric(0),
                        sample=character(0)
                    )
                    meta <- NULL
                }
                else {
                    #tab <- as.matrix(tab)
                    switch(input$rnaCustomMeasureRadio,
                        raw = {
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio))
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio))
                        },
                        rpkm = {
                            tab <- round(edgeR::rpkm(tab,gene.length=le,
                                lib.size=unlist(D$libsize[colnames(tab)]),),
                                digits=6)
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio),digits=6)
                        },
                        rpgm = {
                            for (ss in colnames(tab))
                                tab[,ss] <- round(tab[,ss]/le,digits=6)
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio),digits=6)
                        }
                    )
                    switch(input$rnaCustomScaleRadio,
                        natural = {
                            tab <- tab
                            average <- average
                            deviation <- deviation
                        },
                        log2 = {
                            tab <- round(log2(tab+1),digits=6)
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio),digits=6)
                        }
                    )
                    # Use alternative ids if present
                    if (!is.null(M$alt_id))
                        colnames(tab) <- as.character(M[which(
                            as.character(M$class)==x),"alt_id"])
                    # Attach some annotation
                    meta <- data.frame(
                        name=rownames(tab),
                        average=average,
                        deviation=deviation
                    )
                    currentCustomRnaTables$displayTables[[x]] <- cbind(meta,tab)
                }
                output[[paste("customCountRnaTable",x,sep="_")]] <- 
                    DT::renderDataTable(
                        cbind(meta,tab),
                        class="display compact",
                        rownames=FALSE,
                        options=list(
                            searchHighlight=TRUE,
                            pageLength=10,
                            lengthMenu=c(10,20,50,100)
                        )
                    )
            },M,le,loadedData[[s]][[d]])
            lapply(c,function(x,s,d) {
                list(
                    div(
                        class="small table-container",
                        h4(paste("Expression for ",x)),
                        DT::dataTableOutput(paste("customCountRnaTable",x,
                            sep="_")),
                        div(
                            style="display:inline-block; margin:5px;",
                            actionButton(
                                paste("clearCustomRnaCountSelection_",x,
                                    sep=""),
                                paste("Clear selection"),
                                class="btn-xs"
                            )
                        ),
                        div(
                            style="display:inline-block; margin:5px;",
                            actionButton(
                                paste("invertCustomRnaCountSelection_",x,
                                    sep=""),
                                paste("Invert selection"),
                                class="btn-xs"
                            )
                        ),
                        div(
                            style="display:inline-block; margin:5px;",
                            downloadButton(
                                paste("exportCustomRnaCountSelection_",x,
                                    sep=""),
                                paste("Export selection"),
                                class="btn-xs"
                            )
                        ),
                        div(
                            style="display:inline-block; margin:5px;",
                            downloadButton(
                                paste("exportCustomRnaCountAll_",x,sep=""),
                                paste("Export all"),
                                class="btn-xs"
                            )
                        )
                    ),
                    hr()
                )
            },s,d)
        }
    })
    
    output$setChrsExpr <- renderUI({
        selectizeInput(
            inputId="customChrExpr",
            label="", 
            choices=c("",
                getValidChromosomes("hg19")
            ),
            options=list(
                placeholder="Chrom"
            )
        )
    })
    
}

expressionCalculatorTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    customRegions <- allReactiveVars$customRegions
    customRnaRegions <- allReactiveVars$customRnaRegions
    
    expressionCalculatorTabPanelReactiveEvents <- 
        expressionCalculatorTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    addTheCustomRnaRegion <- 
        expressionCalculatorTabPanelReactiveEvents$addTheCustomRnaRegion
    removeTheCustomRnaRegion <- 
        expressionCalculatorTabPanelReactiveEvents$removeTheCustomRnaRegion
    calcCustomRnaCounts <- 
        expressionCalculatorTabPanelReactiveEvents$calcCustomRnaCounts
    
    expressionCalculatorTabPanelReactiveExprs <- 
        expressionCalculatorTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    createCurrentCustomRnaTables <- 
        expressionCalculatorTabPanelReactiveExprs$createCurrentCustomRnaTables
    handleCustomRnaSampleSelection <-
        expressionCalculatorTabPanelReactiveExprs$handleCustomRnaSampleSelection
    handleCustomRnaDownloadSelection <-
        expressionCalculatorTabPanelReactiveExprs$handleCustomRnaDownloadSelection
        
    expressionCalculatorTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        addTheCustomRnaRegion()
        removeTheCustomRnaRegion()
    })
    
    observe({
        handleCustomRnaSampleSelection()
        handleCustomRnaDownloadSelection()
    })
    
    observe({
        if (length(customRnaRegions$name)!=0 
            || (input$customExpressionFromGeneExplorer=="fromge"
                && length(customRegions$name)!=0))
            shinyjs::enable("calculateCustomRegionRna")
        else
            shinyjs::disable("calculateCustomRegionRna")
    })
    
    observe({
        if (length(customRnaRegions$name)==0)
            shinyjs::disable("removeRnaCustomRegion")
        else
            shinyjs::enable("removeRnaCustomRegion")
    })
    
    observe({
        tryCatch({
            shinyjs::disable("calculateCustomRegionRna")
            calcCustomRnaCounts()
        },error=function(e) {
            return()
        },
        finally={
            if (length(customRnaRegions$name)==0)
                shinyjs::disable("calculateCustomRegionRna")
            else
                shinyjs::enable("calculateCustomRegionRna")
        })
    })
}
