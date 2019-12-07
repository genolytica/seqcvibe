geneSignalTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentGenes <- allReactiveVars$currentGenes
    customRegions <- allReactiveVars$customRegions
    currentOpts <- allReactiveVars$currentOpts
    genePlots <- allReactiveVars$genePlots
    geneExplorerMessages <- allReactiveMsgs$geneExplorerMessages

    updateCurrentGenes <- eventReactive(input$geneGeneName,{
        g <- input$geneGeneName
        dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
        if (length(g)>=length(currentGenes$genes)) {
            currentGenes$genes <- g
            if (!isEmpty(g)) {
                mm <- g[length(g)]
                currentGenes$coords$chromosome <- 
                    c(currentGenes$coords$chromosome,
                        as.character(seqnames(dbGene[mm]))[1])
                currentGenes$coords$start <- 
                    c(currentGenes$coords$start,start(dbGene[mm]))
                currentGenes$coords$end <- 
                    c(currentGenes$coords$end,end(dbGene[mm]))
                currentGenes$coords$strand <- 
                    c(currentGenes$coords$strand,
                        as.character(strand(dbGene[mm]))[1])
                currentGenes$coords$name <- 
                    c(currentGenes$coords$name,
                        as.character(dbGene[mm]$gene_name))
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Gene ",g[length(g)],
                        " (",as.character(dbGene[mm]$gene_name),") ",
                        "added to the known gene list! Number of genes ",
                        "to plot is now ",length(currentGenes$genes),sep="")
                )
            }
        }
        else if (length(g)<length(currentGenes$genes) || isEmpty(g)) {
            if (isEmpty(g)) {
                currentGenes$genes <- NULL
                currentGenes$coords$chromosome <- NULL
                currentGenes$coords$start <- NULL
                currentGenes$coords$end <- NULL
                currentGenes$coords$strand <- NULL
                currentGenes$coords$name <- NULL
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="WARNING",
                    msg=paste(getTime("WARNING"),"All genes removed from ",
                        "the known gene list!","sep=")
                )
            }
            else {
                prev <- currentGenes$genes
                removed <- setdiff(prev,g)
                ii <- match(removed,prev)
                gii <- currentGenes$coords$name[ii]
                currentGenes$genes <- g
                currentGenes$coords$chromosome <- 
                        currentGenes$coords$chromosome[-ii]
                currentGenes$coords$start <- currentGenes$coords$start[-ii]
                currentGenes$coords$end <- currentGenes$coords$end[-ii]
                currentGenes$coords$strand <- 
                    currentGenes$coords$strand[-ii]
                currentGenes$coords$name <- currentGenes$coords$name[-ii]
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="WARNING",
                    msg=paste(getTime("WARNING"),"Gene ",removed," (",gii,
                        ") removed from the known gene list! Number of ",
                        "gene to plot is now ",length(currentGenes$genes),
                        sep="")
                )
            }
        }
    })
    
    updateFlanks <- eventReactive({
        input$upstreamFlank
        input$downstreamFlank
    },{
        flankValidate <- c(
            upstream=FALSE,
            downstream=FALSE
        )
        errMsg <- c(
            upstream="",
            downstream=""
        )
        u <- input$upstreamFlank
        d <- input$downstreamFlank
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
            output$regionFlankError <- renderUI({
                div(class="error-message",paste(errMsg,collapse="\n"))
            })
        }
        else {
            output$regionFlankError <- renderUI({div()})
            if (u!=currentOpts$flank[1]) {
                currentOpts$flank[1] <- u
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Upstream flanking region",
                        " changed! New upstream flanking: ",u,sep="")
                )
            }
            if (d!=currentOpts$flank[2]) {
                currentOpts$flank[2] <- d
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Downstream flanking ",
                        "region changed! New downstream flanking: ",d,
                        sep="")
                )
            }
        }
    })
    
    addTheCustomRegion <- eventReactive(input$addCustomRegion,{
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
        if (input$customName=="") {
            regionValidate["name"] <- TRUE
            errMsg["name"] <- "A region name is required!"
        
        }
        if (input$customStrand=="") {
            regionValidate["strand"] <- TRUE
            errMsg["strand"] <- "A region strand is required!"
        }
        # HACK: Temporarily removing this validator
        # This selectize input is a rendered UI so it delays takes some more 
        # time to be rendered on Restore
        # As a result input$customChr won't be set by the time this validator 
        # runs and the app breaks
        #
        
        # if (input$customChr=="") {
        #     regionValidate["chrom"] <- TRUE
        #     errMsg["chrom"] <- "A region chromosome is required!"
        # }
        
        if (input$customStart=="" || as.integer(input$customStart)<=0
            || is.na(suppressWarnings(as.numeric(input$customStart)))) {
            regionValidate["start"] <- TRUE
            errMsg["start"] <- 
                "Region start base must be a valid positive integer!"
        }
        if (input$customEnd=="" || as.numeric(input$customEnd)<=0
            || is.na(suppressWarnings(as.integer(input$customEnd)))) {
            regionValidate["end"] <- TRUE
            errMsg["end"] <- 
                "Region end base must be a valid positive integer!"
        }
        if (!(regionValidate["start"] || regionValidate["end"])) {
            # Check start before end
            s <- as.integer(input$customStart)
            e <- as.integer(input$customEnd)
            if (e-s<=0) {
                regionValidate["startBeforeEnd"] <- TRUE
                errMsg["startBeforeEnd"] <- paste("Region start must be ",
                    "smaller than region end and the region length must ",
                    "be greater than 1!",sep="")
            }
            
            # Check the new region does not already exists
            cc <- customRegions$chromosome
            ss <- customRegions$start
            ee <- customRegions$end
            tt <- customRegions$strand
            if (input$customName %in% customRegions$name) {
                regionValidate["regionExistsName"] <- TRUE
                errMsg["regionExistsName"] <- paste("Custom regions must ",
                    "have each a unique name! ",input$customName," ",
                    "already exists.",sep="")
            }
            if (length(ss)>0) {
                for (i in 1:length(ss)) {
                    if (input$customChr==cc[i] && input$customStart==ss[i] 
                        && input$customEnd==ee[i] 
                        && input$customStrand==tt[i])
                    regionValidate["regionExistsLocus"] <- TRUE
                    errMsg["regionExistsLocus"] <- 
                        "Custom regions must have unique coordinates!"
                }
            }
        }
        
        if (any(regionValidate)) {
            output$customRegionError <- renderUI({
                div(class="error-message",paste(errMsg,collapse="\n"))
            })
            return("")
        }
        else {
            output$customRegionError <- renderUI({div()})
            customRegions$chromosome <- c(customRegions$chromosome,
                input$customChr)
            customRegions$start <- c(customRegions$start,input$customStart)
            customRegions$end <- c(customRegions$end,input$customEnd)
            customRegions$strand <- c(customRegions$strand,
                input$customStrand)
            customRegions$name <- c(customRegions$name,input$customName)
            
            geneExplorerMessages <- updateMessages(
                geneExplorerMessages,
                type="SUCCESS",
                msg=paste(getTime("SUCCESS"),"Custom region ",
                    input$customName,
                    " added to the custom regions to plot! Number of ",
                    "custom regions to plot is now ",
                    length(customRegions$name),sep="")
            )
        }
    })
    
    removeTheCustomRegion <- eventReactive(input$removeCustomRegion,{
        output$customRegionError <- renderUI({div()})
        j <- as.integer(input$customRegionList_rows_selected)
        if (length(j)>0) {
            ii <- customRegions$name[j]
            customRegions$chromosome <- customRegions$chromosome[-j]
            customRegions$start <- customRegions$start[-j]
            customRegions$end <- customRegions$end[-j]
            customRegions$strand <- customRegions$strand[-j]
            customRegions$name <- customRegions$name[-j]
            geneExplorerMessages <- updateMessages(
                geneExplorerMessages,
                type="WARNING",
                msg=paste(getTime("WARNING"),"Custom region ",ii," ",
                    " removed from the custom regions to plot! Number of ",
                    "custom regions to plot is now ",
                    length(customRegions$name),sep="")
            )
        }
    })
    
    createGeneProfile <- eventReactive(input$createGeneProfile,{
        if (is.null(currentMetadata$final)) {
            output$geneExplorerError <- renderUI({
                div(class="error-message",paste("You must create a ",
                    "dataset first!",sep=""))
            })
        }
        else {
            output$geneExplorerError <- renderUI({div()})
            dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
            isolate(input$geneType)
            knownGenes <- as.list(currentGenes$genes)
            names(knownGenes) <- currentGenes$genes
            customGenes <- makeGRangesFromDataFrame(
                df=data.frame(
                    chromosome=as.character(customRegions$chromosome),
                    start=as.numeric(customRegions$start),
                    end=as.numeric(customRegions$end),
                    strand=as.character(customRegions$strand),
                    gene_id=as.character(customRegions$name)
                ),
                keep.extra.columns=TRUE
            )
            customGenes <- as.list(split(customGenes, as.factor(customGenes)))
            names(customGenes) <- as.character(customRegions$name)
            theGenes <- c(knownGenes,customGenes)
            progress <- shiny::Progress$new(
                session,
                min=0,
                max=length(unique(as.character(
                    currentMetadata$finalclass)))*length(theGenes)*2+2
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
            
            ggProf <- tryCatch(
                getProfile(
                    gene=theGenes,
                    flank=as.integer(currentOpts$flank),
                    source=as.character(currentMetadata$source),
                    dataset=as.character(currentMetadata$dataset),
                    #class=as.character(currentMetadata$class),
                    class=unique(as.character(currentMetadata$final$class)),
                    sumStat=as.character(currentOpts$sumStat),
                    config=currentMetadata$final,
                    dbGene=dbGene,
                    trim=as.numeric(currentOpts$trim),
                    pathPrefix=appConfig$paths$data,
                    messageContainer=geneExplorerMessages,
                    progressFun=updateProgress,
                    rc=RC
                ),
                #warning=function(w) {
                #    geneExplorerMessages <- updateMessages(
                #        geneExplorerMessages,
                #        type="WARNING",
                #        msg=paste(getTime("WARNING"),"A warning occured ",
                #            "while calculating profiles from bigWig ",
                #            "files. Please contact the administrator ",
                #            "stating the following warning: ",w,sep="")
                #    )
                #},
                error=function(e) {
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="WARNING",
                        msg=paste(getTime("WARNING"),"An error occured ",
                            "while calculating profiles from bigWig ",
                            "files. Will try with BAM if present. Please ",
                            "contact the administrator stating the ",
                            "following error: ",e,sep="")
                    )
                    tryCatch(
                        getProfile(
                            gene=theGenes,
                            flank=as.integer(currentOpts$flank),
                            source=as.character(currentMetadata$source),
                            dataset=as.character(currentMetadata$dataset),
                            #class=as.character(currentMetadata$class),
                            class=unique(as.character(
                                currentMetadata$final$class)),
                            sumStat=as.character(currentOpts$sumStat),
                            config=currentMetadata$final,
                            dbGene=dbGene,
                            trim=as.numeric(currentOpts$trim),
                            fromBam=TRUE,
                            pathPrefix=appConfig$paths$data,
                            messageContainer=geneExplorerMessages,
                            progressFun=updateProgress,
                            rc=RC
                        ),
                        #warning=function(w) {
                        #    geneExplorerMessages <- updateMessages(
                        #        geneExplorerMessages,
                        #        type="WARNING",
                        #        msg=paste(getTime("WARNING"),"A warning ",
                        #            "occured while calculating profiles ",
                        #            "from BAM files. Please contact the ",
                        #            "administrator stating the following ",
                        #            "warning: ",w,sep="")
                        #    )
                        #},
                        error=function(e) {
                            geneExplorerMessages <- updateMessages(
                                geneExplorerMessages,
                                type="ERROR",
                                msg=paste(getTime("ERROR"),"An error ",
                                    "occured while calculating profiles ",
                                    "from BAM files. Please contact the ",
                                    "administrator stating the following ",
                                    "error: ",e,sep="")
                            )
                        },finally=""
                    )
                },finally=""
            )
            
            ggProf <- ggProf + 
                facet_wrap(~ Locus,scales="free") +
                scale_color_manual(values=currentOpts$colours) + 
                scale_fill_manual(values=currentOpts$colours)
            genePlots$geneProfile <- ggProf
            genePlots$rendered <- TRUE
        }
    })
    
    return(list(
        updateCurrentGenes=updateCurrentGenes,
        updateFlanks=updateFlanks,
        addTheCustomRegion=addTheCustomRegion,
        removeTheCustomRegion=removeTheCustomRegion,
        createGeneProfile=createGeneProfile
    ))

}

geneSignalTabPanelReactive <- function(input,output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentGenes <- allReactiveVars$currentGenes
    customRegions <- allReactiveVars$customRegions
    currentOpts <- allReactiveVars$currentOpts
    genePlots <- allReactiveVars$genePlots
    geneExplorerMessages <- allReactiveMsgs$geneExplorerMessages
    
    updateClassColours <- reactive({
        c <- currentMetadata$class
        currentOpts$colours <- baseColours[c]
    })
    
    updateGeneColours <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            observeEvent(input[[paste("geneColour_",x,sep="")]],{
                newc <- input[[paste("geneColour_",x,sep="")]]
                if (!is.null(currentOpts$colours[x]) 
                    && newc!=currentOpts$colours[x]) {
                    currentOpts$colours[x] <- newc
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Class ",x," colour ",
                        "changed! New colour: ",newc,sep="")
                    )
                }
            })
        })
    })
    
    updateGeneNames <- reactive({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"geneGeneName",
                choices=NULL,
                server=TRUE)
        else {
            g <- isolate({input$geneGeneName})
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            i <- grep(paste0("^",paste(g,collapse="|")),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"geneGeneName",
                    choices=geneNames,
                    selected=g,
                    server=TRUE)
            }
        }
    })
    
    return(list(
        updateClassColours=updateClassColours,
        updateGeneColours=updateGeneColours,
        updateGeneNames=updateGeneNames
    ))

}

geneSignalTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentGenes <- allReactiveVars$currentGenes
    customRegions <- allReactiveVars$customRegions
    genePlots <- allReactiveVars$genePlots
    geneExplorerMessages <- allReactiveMsgs$geneExplorerMessages
    
    # Very dirty hack to update reactive values holding gene coordinates as
    # there is a problem when selectize input get empty, see also this:
    # http://stackoverflow.com/questions/26803536/shiny-how-to-update-a-reactivevalues-object
    dirtyHack <- reactive({
        currentGenes$genes <- NULL
        currentGenes$coords$chromosome <- NULL
        currentGenes$coords$start <- NULL
        currentGenes$coords$end <- NULL
        currentGenes$coords$strand <- NULL
        currentGenes$coords$name <- NULL
        geneExplorerMessages <- updateMessages(
            geneExplorerMessages,
            type="WARNING",
            msg=paste(getTime("WARNING"),"All genes removed from ",
                "the known gene list!",sep="")
        )
        return(data.frame(
            chromosome=character(),
            start=integer(),
            end=integer(),
            strand=character(),
            name=character()
        ))
    })
    
    output$knownGeneList <- DT::renderDataTable(
        if (isEmpty(input$geneGeneName))
            isolate(dirtyHack())
        else
            data.frame(
                chromosome=currentGenes$coords$chromosome,
                start=currentGenes$coords$start,
                end=currentGenes$coords$end,
                strand=currentGenes$coords$strand,
                name=currentGenes$coords$name
            ),
        class="display compact",
        rownames=FALSE,
        options=list(
            dom='t',
            paging=FALSE
        )
    )
    
    output$customRegionList <- DT::renderDataTable(
        data.frame(
            chromosome=customRegions$chromosome,
            start=customRegions$start,
            end=customRegions$end,
            strand=customRegions$strand,
            name=customRegions$name
        ),
        class="display compact",
        rownames=FALSE,
        options=list(
            dom='t',
            paging=FALSE#,
            #columnDefs=list(list(
            #    targets=0,
            #    visible=FALSE,
            #    searchable=FALSE
            #))
        )
    )
        
    output$setChrs <- renderUI({
        selectizeInput(
            inputId="customChr",
            label="", 
            choices=c("",
                getValidChromosomes(currentMetadata$genome)
            ),
            options=list(
                placeholder="Chrom"
            )
        )
    })
    
    output$geneExplorerColours <- renderUI({
        if (!is.null(currentMetadata$final)) {
            c <- unique(as.character(currentMetadata$final$class))
            lapply(1:length(c),function(i,c) {
                colourInput(
                    inputId=paste("geneColour_",c[i],sep=""),
                    label=paste("Select colour for",c[i]),
                    value=baseColours[i]
                )
            },c)
        }
    })
    
    output$geneExplorerMessages <- renderUI({
        lapply(geneExplorerMessages$messages,function(x) {
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
    
    output$geneProfile <- renderPlot({
        genePlots$geneProfile
    })
    
    output$exportGenePDF <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("gene_plot_",tt,".pdf", sep='')
        },
        content=function(con) {
            ggsave(filename=con,plot=genePlots$geneProfile,
                width=10,height=7)
        }
    )
    
    output$exportGenePNG <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("gene_plot_",tt,".png", sep='')
        },
        content=function(con) {
            ggsave(filename=con,plot=genePlots$geneProfile,
                width=10,height=7)
        }
    )
    
    output$exportGeneGG2 <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("gene_plot_",tt,".rda", sep='')
        },
        content=function(con) {
            gg <- genePlots$geneProfile
            save(gg,file=con)
        }
    )
}

geneSignalTabPanelObserve <- function(input,output,session,allReactiveVars,
    allReactiveMsgs) {
    currentOpts <- allReactiveVars$currentOpts
    customRegions <- allReactiveVars$customRegions
    genePlots <- allReactiveVars$genePlots
        
    geneSignalTabPanelReactiveEvents <- 
        geneSignalTabPanelEventReactive(input,output,session,
                allReactiveVars,allReactiveMsgs)
        
    updateCurrentGenes <- 
        geneSignalTabPanelReactiveEvents$updateCurrentGenes
    updateFlanks <- geneSignalTabPanelReactiveEvents$updateFlanks
    addTheCustomRegion <- 
        geneSignalTabPanelReactiveEvents$addTheCustomRegion
    removeTheCustomRegion <- 
        geneSignalTabPanelReactiveEvents$removeTheCustomRegion
    createGeneProfile <- 
        geneSignalTabPanelReactiveEvents$createGeneProfile
    
    geneSignalTabPanelReactiveExprs <- 
        geneSignalTabPanelReactive(input,output,session,allReactiveVars,
            allReactiveMsgs)
    
    updateClassColours <- geneSignalTabPanelReactiveExprs$updateClassColours
    updateGeneColours <- geneSignalTabPanelReactiveExprs$updateGeneColours
    updateGeneNames <- geneSignalTabPanelReactiveExprs$updateGeneNames
    
    geneSignalTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        if (isEmpty(input$geneGeneName) && length(customRegions$name)==0)
            shinyjs::disable("createGeneProfile")
        else
            shinyjs::enable("createGeneProfile")
    })
    
    observe({
        # To elminate danger in certain bookmark restores
        if (!isEmpty(input$geneSumStatType)) { 
            if (input$geneSumStatType!="trimmed") {
                shinyjs::disable("geneTrimPct")
                if (isEmpty(input$geneGeneName) 
                    && length(customRegions$name)==0)
                    shinyjs::disable("createGeneProfile")
                else
                    shinyjs::enable("createGeneProfile")
            }
            else {
                shinyjs::enable("geneTrimPct")
                trimp <- as.numeric(input$geneTrimPct)
                if (is.na(trimp) || trimp<0 || trimp>0.5) {
                    output$geneExplorerError <- renderUI({
                        div(class="error-message",paste("The trimming ",
                            "must be a number between 0 and 0.5!",sep=""))
                    })
                    shinyjs::disable("createGeneProfile")
                }
                else {
                    output$geneExplorerError <- renderUI({div()})
                    if (isEmpty(input$geneGeneName) 
                        && length(customRegions$name)==0)
                        shinyjs::disable("createGeneProfile")
                    else
                        shinyjs::enable("createGeneProfile")
                    currentOpts$trim <- trimp
                }
            }
        }
    })
    
    observe({
        if (length(customRegions$name)==0)
            shinyjs::disable("removeCustomRegion")
        else
            shinyjs::enable("removeCustomRegion")
    })
    
    observe({
        updateGeneNames()
        updateCurrentGenes()
    })
    
    observe({
        addTheCustomRegion()
        removeTheCustomRegion()
    })
    
    observe({
        updateFlanks()
        updateClassColours()
        updateGeneColours()
    })
    
    observe({
        tryCatch({
            shinyjs::disable("createGeneProfile")
            createGeneProfile()
        },error=function(e) {
            genePlots$rendered <- FALSE
            return()
        },
        finally={
            if (isEmpty(input$geneGeneName) 
                && length(customRegions$name)==0)
                shinyjs::disable("createGeneProfile")
            else
                shinyjs::enable("createGeneProfile")
        })
    })
}
