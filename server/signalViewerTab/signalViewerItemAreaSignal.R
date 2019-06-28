areaSignalTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    customRegions <- allReactiveVars$customRegions
    currentOpts <- allReactiveVars$currentOpts
    customArea <- allReactiveVars$customArea
    customRegionsInArea <- allReactiveVars$customRegionsInArea
    areaPlots <- allReactiveVars$areaPlots
    areaExplorerMessages <- allReactiveMsgs$areaExplorerMessages
    
    updateCurrentArea <- eventReactive({
        input$customChrA
        input$customStartA
        input$customEndA
    },{
        customArea$chromosome <- input$customChrA
        customArea$start <- input$customStartA
        customArea$end <- input$customEndA
        customArea$name <- paste(customArea$chromosome,":",
            customArea$start,"-",customArea$end,sep="")
        customArea$strand <- "*"
        if (!any(
            isEmpty(input$customChrA),
            isEmpty(input$customStartA),
            isEmpty(input$customEndA),
            as.integer(input$customStartA)>=as.integer(input$customEndA)
        ))
            areaExplorerMessages <- updateMessages(
                areaExplorerMessages,
                type="INFO",
                msg=paste(getTime("INFO"),"Area to plot: ",customArea$name,
                    sep="")
            )
    })
    
    updateFlanksA <- eventReactive({
        input$upstreamFlankA
        input$downstreamFlankA
    },{
        flankValidate <- c(
            upstream=FALSE,
            downstream=FALSE
        )
        errMsg <- c(
            upstream="",
            downstream=""
        )
        u <- input$upstreamFlankA
        d <- input$downstreamFlankA
        if (u!=currentOpts$flank[1]) {
            u <- as.integer(u)
            if (is.na(suppressWarnings(u)) || isEmpty(u) || u<=0) {
                flankValidate["upstream"] <- TRUE
                errMsg["upstream"] <- paste("Upstream flanking region",
                    "must be a positive integer!")
            }
        }
        if (d!=currentOpts$flank[2]) {
            d <- as.integer(d)
            if (is.na(suppressWarnings(d)) || isEmpty(d) || d<=0) {
                flankValidate["downstream"] <- TRUE
                errMsg["downstream"] <- paste("Downstream flanking region",
                    "must be a positive integer!")
            }
        }
        if (any(flankValidate)) {
            output$areaFlankError <- renderUI({
                div(class="error-message",paste(errMsg,collapse="\n"))
            })
        }
        else {
            output$areaFlankError <- renderUI({div()})
            if (u!=currentOpts$flank[1]) {
                currentOpts$flank[1] <- u
                areaExplorerMessages <- updateMessages(
                    areaExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Upstream flanking region",
                        " changed! New upstream flanking: ",u,sep="")
                )
            }
            if (d!=currentOpts$flank[2]) {
                currentOpts$flank[2] <- d
                areaExplorerMessages <- updateMessages(
                    areaExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Downstream flanking ",
                        "region changed! New downstream flanking: ",d,
                        sep="")
                )
            }
        }
    })
    
    addTheCustomRegionInArea <- eventReactive(input$addCustomRegionInArea,{
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
        if (input$customNameInArea=="") {
            regionValidate["name"] <- TRUE
            errMsg["name"] <- "A region name is required!"
        
        }
        if (input$customStrandInArea=="") {
            regionValidate["strand"] <- TRUE
            errMsg["strand"] <- "A region strand is required!"
        }
        if (input$customChrInArea=="") {
            regionValidate["chrom"] <- TRUE
            errMsg["chrom"] <- "A region chromosome is required!"
        
        }
        if (input$customStartInArea=="" 
            || as.integer(input$customStartInArea)<=0
            || is.na(suppressWarnings(
                as.numeric(input$customStartInArea)))) {
            regionValidate["start"] <- TRUE
            errMsg["start"] <- 
                "Region start base must be a valid positive integer!"
        }
        if (input$customEndInArea=="" 
            || as.numeric(input$customEndInArea)<=0
            || is.na(suppressWarnings(
                as.integer(input$customEndInArea)))) {
            regionValidate["end"] <- TRUE
            errMsg["end"] <- 
                "Region end base must be a valid positive integer!"
        }
        if (!(regionValidate["start"] || regionValidate["end"])) {
            # Check start before end
            s <- as.integer(input$customStartInArea)
            e <- as.integer(input$customEndInArea)
            if (e-s<=0) {
                regionValidate["startBeforeEnd"] <- TRUE
                errMsg["startBeforeEnd"] <- paste("Region start must be ",
                    "smaller than region end and the region length must ",
                    "be greater than 1!",sep="")
            }
            
            # Check the new region does not already exists
            cc <- customRegionsInArea$chromosome
            ss <- customRegionsInArea$start
            ee <- customRegionsInArea$end
            tt <- customRegionsInArea$strand
            if (input$customNameInArea %in% customRegionsInArea$name) {
                regionValidate["regionExistsName"] <- TRUE
                errMsg["regionExistsName"] <- paste("Custom regions must ",
                    "have each a unique name! ",input$customNameInArea," ",
                    "already exists.",sep="")
            }
            if (length(ss)>0) {
                for (i in 1:length(ss)) {
                    if (input$customChrInArea==cc[i] 
                        && input$customStartInArea==ss[i] 
                        && input$customEndInArea==ee[i] 
                        && input$customStrandInArea==tt[i])
                    regionValidate["regionExistsLocus"] <- TRUE
                    errMsg["regionExistsLocus"] <- 
                        "Custom regions must have unique coordinates!"
                }
            }
        }
        
        if (any(regionValidate)) {
            output$customRegionErrorInArea <- renderUI({
                div(class="error-message",paste(errMsg,collapse="\n"))
            })
            return("")
        }
        else {
            output$customRegionErrorInArea <- renderUI({div()})
            customRegionsInArea$chromosome <- c(customRegionsInArea$chr,
                input$customChrInArea)
            customRegionsInArea$start <- c(customRegionsInArea$start,
                input$customStartInArea)
            customRegionsInArea$end <- c(customRegionsInArea$end,
                input$customEndInArea)
            customRegionsInArea$strand <- c(customRegionsInArea$strand,
                input$customStrandInArea)
            customRegionsInArea$name <- c(customRegionsInArea$name,
                input$customNameInArea)
            
            areaExplorerMessages <- updateMessages(
                areaExplorerMessages,
                type="SUCCESS",
                msg=paste(getTime("SUCCESS"),"Custom region ",
                    input$customNameInArea,
                    " added to the search of custom area to plot! Number",
                    " of custom regions to plot is now ",
                    length(customRegionsInArea$name),sep="")
            )
        }
    })
    
    removeTheCustomRegionInArea <- eventReactive(
        input$removeCustomRegionInArea,{
        output$customRegionErrorInArea <- renderUI({div()})
        j <- as.integer(input$customRegionListInArea_rows_selected)
        if (length(j)>0) {
            ii <- customRegionsInArea$name[j]
            customRegionsInArea$chromosome <- 
                customRegionsInArea$chromosome[-j]
            customRegionsInArea$start <- customRegionsInArea$start[-j]
            customRegionsInArea$end <- customRegionsInArea$end[-j]
            customRegionsInArea$strand <- customRegionsInArea$strand[-j]
            customRegionsInArea$name <- customRegionsInArea$name[-j]
            areaExplorerMessages <- updateMessages(
                areaExplorerMessages,
                type="WARNING",
                msg=paste(getTime("WARNING"),"Custom region ",ii," ",
                    " removed from the search of custom area to plot! ",
                    "Number of custom regions to plot is now ",
                    length(customRegionsInArea$name),sep="")
            )
        }
    })
    
    createAreaProfile <- eventReactive(input$createAreaProfile,{
        if (is.null(currentMetadata$final)) {
            output$areaExplorerError <- renderUI({
                div(class="error-message",paste("You must create a ",
                    "dataset first!",sep=""))
            })
        }
        else {
            if (!is.na(currentMetadata$genome) &&
                is.null(loadedGenomes[[currentMetadata$genome]]$dbExon)) {
                load(file.path("genome",currentMetadata$genome,
                    "summarized_exon.rda"))
                loadedGenomes[[currentMetadata$genome]]$dbExon <<- sexon
                dataSelectorMessages <- updateMessages(
                    dataSelectorMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Genome ",
                        currentMetadata$genome," exons loaded!", sep="")
                )
            }
            isolate(input$areaTypeRadio)
            output$areaExplorerError <- renderUI({div()})
            dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
            dbExon <- loadedGenomes[[currentMetadata$genome]]$dbExon
            if (input$areaTypeRadio=="gene") {
                gn <- as.character(input$areaGeneName)
                fs <- as.integer(input$upstreamFlankA)
                fd <- as.integer(input$downstreamFlankA)
                dbgi <- dbGene[gn]
                refArea <- list(
                    chr=as.character(seqnames(dbgi))[1],
                    start=start(dbgi)-fs,
                    end=end(dbgi)+fd
                )
                customTrans <- NULL
            }
            else {
                refArea <- list(
                    chr=as.character(customArea$chromosome),
                    start=as.integer(customArea$start),
                    end=as.integer(customArea$end)
                )
                if (input$customRegionsFromGeneExplorer=="fromge") {
                    customRegionsInArea <- customRegions
                    ii <- which(customRegions$chromosome==refArea$chr)
                    if (length(ii)>0) {
                        customRegionsInArea$chromosome <- 
                            customRegionsInArea$chromosome[ii]
                        customRegionsInArea$start <- 
                            customRegionsInArea$start[ii]
                        customRegionsInArea$end <- 
                            customRegionsInArea$end[ii]
                        customRegionsInArea$strand <- 
                            customRegionsInArea$strand[ii]
                        customRegionsInArea$name <- 
                            customRegionsInArea$name[ii]
                    }
                    else {
                        customRegionsInArea$chromosome <- character(0)
                        customRegionsInArea$start <- integer(0)
                        customRegionsInArea$end <- integer(0)
                        customRegionsInArea$strand <- character(0)
                        customRegionsInArea$name <- character(0)
                    }
                }
                if (length(customRegionsInArea$chromosome)>0) {
                    customTrans <- as.list(makeGRangesFromDataFrame(
                        df=data.frame(
                            chromosome=
                                as.character(
                                    customRegionsInArea$chromosome),
                            start=as.integer(customRegionsInArea$start),
                            end=as.integer(customRegionsInArea$end),
                            strand=as.character(customRegionsInArea$strand),
                            gene_id=as.character(customRegionsInArea$name)
                        ),
                        keep.extra.columns=TRUE
                    ))
                    names(customTrans) <- as.character(customRegions$name)
                }
                else
                    customTrans <- NULL
            }
            
            progress <- shiny::Progress$new()
            progress$initialize(
                session,
                min=0,
                max=length(unique(as.character(
                    currentMetadata$finalclass)))*2+1
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
            
            ggTrack <- tryCatch(
                getTrack(
                    refArea=refArea,
                    customGene=customTrans,
                    source=as.character(currentMetadata$source),
                    dataset=as.character(currentMetadata$dataset),
                    #class=as.character(currentMetadata$class),
                    class=unique(as.character(currentMetadata$final$class)),
                    sumStat=as.character(currentOpts$sumStat),
                    config=currentMetadata$final,
                    dbGene=dbGene,
                    dbExon=dbExon,
                    trim=as.numeric(currentOpts$trim),
                    classColours=currentOpts$colours,
                    messageContainer=areaExplorerMessages,
                    progressFun=updateProgress,
                    rc=RC
                ),
                error=function(e) {
                    areaExplorerMessages <- updateMessages(
                        areaExplorerMessages,
                        type="WARNING",
                        msg=paste(getTime("WARNING"),"An error occured ",
                            "while calculating profiles from bigWig ",
                            "files. Will try with BAM if present. Please ",
                            "contact the administrator stating the ",
                            "following error: ",e,sep="")
                    )
                    tryCatch(
                        getTrack(
                            refArea=refArea,
                            customGene=customTrans,
                            source=as.character(currentMetadata$source),
                            dataset=as.character(currentMetadata$dataset),
                            #class=as.character(currentMetadata$class),
                            class=unique(as.character(
                                currentMetadata$final$class)),
                            sumStat=as.character(currentOpts$sumStat),
                            config=currentMetadata$final,
                            dbGene=dbGene,
                            dbExon=dbExon,
                            trim=as.numeric(currentOpts$trim),
                            fromBam=TRUE,
                            classColours=currentOpts$colours,
                            messageContainer=areaExplorerMessages,
                            progressFun=updateProgress,
                            rc=RC
                        ),
                        error=function(e) {
                            areaExplorerMessages <- updateMessages(
                                areaExplorerMessages,
                                type="ERROR",
                                msg=paste(getTime("ERROR"),"An error ",
                                    "occured while calculating profiles ",
                                    "from BAM files. Please contact the ",
                                    "administrator stating the following ",
                                    "error: ",e,sep="")
                            )
                        },
                        finally=""
                    )
                },finally=""
            )
            
            ggTrack <- ggTrack + 
                scale_color_manual(values=currentOpts$colours) + 
                scale_fill_manual(values=currentOpts$colours)
            
            areaPlots$areaProfile <- ggTrack
            areaPlots$rendered <- TRUE
        }
    })
    
    return(list(
        updateCurrentArea=updateCurrentArea,
        updateFlanksA=updateFlanksA,
        addTheCustomRegionInArea=addTheCustomRegionInArea,
        removeTheCustomRegionInArea=removeTheCustomRegionInArea,
        createAreaProfile=createAreaProfile
    ))
}

areaSignalTabPanelReactive <- function(input,output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentOpts <- allReactiveVars$currentOpts
    areaExplorerMessages <- allReactiveMsgs$areaExplorerMessages
    
    updateAreaGeneNames <- reactive({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"areaGeneName",
                choices=NULL,
                server=TRUE)
        g <- isolate({input$areaGeneName})
        geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
        i <- grep(paste0("^",g),geneNames,perl=TRUE)
        if (length(i)>0) {
            updateSelectizeInput(session,"areaGeneName",
                choices=geneNames[i],
                selected=input$areaGeneName,
                server=TRUE)
        }
    })
    
    updateAreaColours <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            observeEvent(input[[paste("areaColour_",x,sep="")]],{
                newc <- input[[paste("areaColour_",x,sep="")]]
                if (newc!=currentOpts$colours[x]) {
                    currentOpts$colours[x] <- newc
                    areaExplorerMessages <- updateMessages(
                        areaExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Class ",x," colour ",
                        "changed! New colour: ",newc,sep="")
                    )
                }
            })
        })
    })
    
    return(list(
        updateAreaGeneNames=updateAreaGeneNames,
        updateAreaColours=updateAreaColours
    ))

}

areaSignalTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    customArea <- allReactiveVars$customArea
    customRegionsInArea <- allReactiveVars$customRegionsInArea
    areaPlots <- allReactiveVars$areaPlots
    areaExplorerMessages <- allReactiveMsgs$areaExplorerMessages
    
    output$customRegionListInArea <- DT::renderDataTable(
        data.frame(
            chromosome=customRegionsInArea$chromosome,
            start=customRegionsInArea$start,
            end=customRegionsInArea$end,
            strand=customRegionsInArea$strand,
            name=customRegionsInArea$name
        ),
        class="display compact",
        rownames=FALSE,
        options=list(
            dom='t',
            paging=FALSE
        )
    )
    
    output$setChrsA <- renderUI({
        selectizeInput(
            inputId="customChrA",
            label="", 
            choices=c("",
                getValidChromosomes("hg19")
            ),
            options=list(
                placeholder="Chrom"
            )
        )
    })
    
    output$areaExplorerColours <- renderUI({
        if (!is.null(currentMetadata$final)) {
            #c <- currentMetadata$class
            c <- unique(as.character(currentMetadata$final$class))
            lapply(1:length(c),function(i,c) {
                colourInput(
                    inputId=paste("areaColour_",c[i],sep=""),
                    label=paste("Select colour for",c[i]),
                    value=baseColours[i]
                )
            },c)
        }
    })
    
    output$areaExplorerMessages <- renderUI({
        lapply(areaExplorerMessages$messages,function(x) {
            switch(x$type,
                INFO = {
                    div(class="info-box",x$msg)
                },
                SUCCESS = {
                    div(class="success-box",x$msg)
                },
                WARNING = {
                    div(class="warn-box",x$msg)
                },
                ERROR = {
                    div(class="error-box",x$msg)
                }
            )
        })
    })
    
    output$areaProfile <- renderPlot({
        areaPlots$areaProfile
    })
    
    output$exportAreaPDF <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("area_plot_",tt,".pdf", sep='')
        },
        content=function(con) {
            ggsave(filename=con,plot=areaPlots$areaProfile,
                width=14,height=7)
        }
    )
    
    output$exportAreaPNG <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("area_plot_",tt,".png", sep='')
        },
        content=function(con) {
            ggsave(filename=con,plot=areaPlots$areaProfile,
                width=14,height=7)
        }
    )
    
    output$exportAreaGG2 <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("area_plot_",tt,".rda", sep='')
        },
        content=function(con) {
            gg <- areaPlots$areaProfile
            save(gg,file=con)
        }
    )
}

areaSignalTabPanelObserve <- function(input,output,session,allReactiveVars,
    allReactiveMsgs) {
    customArea <- allReactiveVars$customArea
    customRegions <- allReactiveVars$customRegions
    customRegionsInArea <- allReactiveVars$customRegionsInArea
    currentOpts <- allReactiveVars$currentOpts
    areaPlots <- allReactiveVars$areaPlots
        
    areaSignalTabPanelReactiveEvents <- 
        areaSignalTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    updateCurrentArea <- areaSignalTabPanelReactiveEvents$updateCurrentArea
    updateFlanksA <- areaSignalTabPanelReactiveEvents$updateFlanksA
    addTheCustomRegionInArea <- 
        areaSignalTabPanelReactiveEvents$addTheCustomRegionInArea
    removeTheCustomRegionInArea <- 
        areaSignalTabPanelReactiveEvents$removeTheCustomRegionInArea
    createAreaProfile <- areaSignalTabPanelReactiveEvents$createAreaProfile
    
    areaSignalTabPanelReactiveExprs <- 
        areaSignalTabPanelReactive(input,output,session,allReactiveVars,
            allReactiveMsgs)
    
    updateAreaGeneNames <- areaSignalTabPanelReactiveExprs$updateAreaGeneNames
    updateAreaColours <- areaSignalTabPanelReactiveExprs$updateAreaColours
    
    areaSignalTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        updateTextInput(session,"customChrInArea",label="",
            value=input$customChrA)
    })
    
    observe({
        if (length(customRegionsInArea$name)==0)
            shinyjs::disable("removeCustomRegionInArea")
        else
            shinyjs::enable("removeCustomRegionInArea")
    })
    
    observe({
        addTheCustomRegionInArea()
        removeTheCustomRegionInArea()
    })
    
    # Trimming textbox
    observe({
        if (input$areaSumStatType!="trimmed") {
            shinyjs::disable("areaTrimPct")
            shinyjs::enable("createAreaProfile")
        }
        else {
            shinyjs::enable("areaTrimPct")
            trimp <- as.numeric(input$areaTrimPct)
            if (is.na(trimp) || trimp<0 || trimp>0.5) {
                output$areaExplorerError <- renderUI({
                    div(class="error-message",paste("The trimming ",
                        "must be a number between 0 and 0.5!",sep=""))
                })
                shinyjs::disable("createAreaProfile")
            }
            else {
                output$areaExplorerError <- renderUI({div()})
                if (isEmpty(input$areaGeneName) 
                    && length(customArea$name)==0)
                    shinyjs::disable("createAreaProfile")
                else
                    shinyjs::enable("createAreaProfile")
                currentOpts$trim <- trimp
            }
        }
    })
    
    observe({
        # Engage button status
        customArea$chromosome <- input$customChrA
        customArea$start <- input$customStartA
        customArea$end <- input$customEndA
        if (any(is.null(customArea$chromosome) || customArea$chromosome=="",
            is.null(customArea$start) || customArea$start=="",
            is.null(customArea$end) || customArea$end=="",
            customArea$start>=customArea$end)
            && input$areaTypeRadio=="area")
            shinyjs::disable("createAreaProfile")
        else
            shinyjs::enable("createAreaProfile")
    })
    
    observe({
        tryCatch({
            shinyjs::disable("createAreaProfile")
            createAreaProfile()
        },error=function(e) {
            areaPlots$rendered <- FALSE
            return()
        },
        finally={
            if (isEmpty(input$areaGeneName) 
                && length(customArea$name)==0)
                shinyjs::disable("createAreaProfile")
            else
                shinyjs::enable("createAreaProfile")
        })
    })
    
    observe({
        updateAreaGeneNames()
        updateFlanksA()
        updateCurrentArea()
        updateAreaColours()
    })
}
