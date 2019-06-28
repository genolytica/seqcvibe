correlationTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    currentCorrelation <- allReactiveVars$currentCorrelation
        
    performRnaCorrelation <- eventReactive(input$performRnaCorrelation,{
        if (is.null(currentMetadata$final)) {
            output$rnaCorrelationSettingsError <- renderUI({
                div(class="error-message",paste("You must create a ",
                    "dataset first!",sep=""))
            })
            return()
        }
        output$rnaCorrelationSettingsError <- renderUI({div()})
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
        D <- loadedData[[s]][[d]]
        genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
        
        # Determine which genes to use
        switch(input$rnaCorrelationGeneList,
            select = {
                g <- input$selectCorrelationGeneName
                if (input$rnaCorrelateWhat=="refgene" 
                    && !isEmpty(input$rnaCorrelationRefGene))
                    g <- unique(c(input$rnaCorrelationRefGene,g))
            },
            expr = {
                if (is.null(currentPipelineOutput$counts)) {
                    output$rnaCorrelationSettingsError <- renderUI({
                        div(class="error-message",paste("You must execute a ",
                            "differential expression analysis first!",sep=""))
                    })
                    return()
                }
                else {
                    output$rnaCorrelationSettingsError <- renderUI({div()})
                    g <- names(which(apply(currentPipelineOutput$flags,1,
                        function(x) all(x==0))))
                }
            },
            custom = {
                g <- input$rnaCorrelationCustomList
                g <- strsplit(g,split="\n")[[1]]
                genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(g,names(genes))
                na <- which(is.na(m))
                if (length(na)>0)
                    m <- m[-na]
                g <- genes[m]
                if (input$rnaCorrelateWhat=="refgene" 
                    && !isEmpty(input$rnaCorrelationRefGene))
                    g <- unique(c(input$rnaCorrelationRefGene,g))
            },
            all = {
                bad <- apply(D$norm,1,function(x) { return(all(x==0)) })
                g <- genes[-which(bad)]
                if (input$rnaCorrelateWhat=="refgene" 
                    && !isEmpty(input$rnaCorrelationRefGene)) {
                        if (!(input$rnaCorrelationRefGene %in% g))
                            g <- unique(c(input$rnaCorrelationRefGene,g))
                    }
            }
        )
        
        # Determine the measurements table
        switch(input$rnaCorrelationMeasureRadio,
            counts = {
                tab <- D$norm[g,samples,drop=FALSE]
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    tab <- rbind(tab,A)
                }
            },
            rpkm = {
                tab <- round(edgeR::rpkm(
                    D$counts[g,samples,drop=FALSE],
                    gene.length=D$length[g],
                    lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(edgeR::rpkm(A,
                        gene.length=currentCustomRnaTables$lengths,
                        lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                    tab <- rbind(tab,A)
                }
            },
            rpgm = {
                tab <- round(D$norm[g,
                    samples,drop=FALSE]/D$length[g],
                        digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(A/currentCustomRnaTables$lengths,
                        digits=6)
                    tab <- rbind(tab,A)
                }
            }
        )
        # We now have the tab, transformations based on other selections
        switch(input$rnaCorrelationScaleRadio,
            natural = {
                tab <- tab
            },
            log2 = {
                tab <- round(log2(tab+1),digits=6)
            }
        )
        if (!is.null(currentMetadata$final$alt_id))
            colnames(tab) <- as.character(currentMetadata$final$alt_id)
        
        currentCorrelation$tooManyGenesAlert <- FALSE
        switch(input$rnaCorrelateWhat,
            samples = {
                currentCorrelation$what <- "samples"
                currentCorrelation$refGene <- NULL
                currentCorrelation$corMatrix <- cor(tab,
                    method=input$rnaCorrelationMethod)
            },
            allgenes = {
                currentCorrelation$what <- "genes"
                currentCorrelation$refGene <- NULL
                if (nrow(tab)>1000)
                    currentCorrelation$tooManyGenesAlert <- TRUE
                else
                    currentCorrelation$corMatrix <- cor(t(tab),
                        method=input$rnaCorrelationMethod)
            },
            refgene = {
                currentCorrelation$what <- "genes"
                currentCorrelation$refGene <- input$rnaCorrelationRefGene
                if (nrow(tab)>1000)
                    currentCorrelation$tooManyGenesAlert <- TRUE
                else
                    currentCorrelation$corMatrix <- cor(t(tab),
                        method=input$rnaCorrelationMethod)
            }
        )
        
        currentCorrelation$datMatrix <- round(tab,digits=3)
        currentCorrelation$opts$method <- input$rnaCorrelationMethod
        currentCorrelation$opts$symm <- input$checkRnaSymmCorHeatmap
        currentCorrelation$opts$colors <- c(
            input$corrColourLow,
            input$corrColourNo,
            input$corrColourHigh
        )
    })
    
    return(list(
        performRnaCorrelation=performRnaCorrelation
    ))
}

correlationTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

correlationTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCorrelation <- allReactiveVars$currentCorrelation
    
    output$correlationOutput <- renderUI({
        if (currentCorrelation$tooManyGenesAlert) {
            output$correlation <- renderPlot({
                ggmessage(msg=paste("Cannot calculate pairwise correlations\n",
                    "for more than 1000 genes asit will\ntake very long.",
                    " Please restrict the number\nof selected\ngenes and try",
                    "again"),type="warning")
            })
            plotOutput("correlation",height="640px")
        }
        else if (is.null(currentCorrelation$corMatrix)) {
            output$correlation <- renderPlot({
                currentCorrelation$entryCor
            })
            plotOutput("correlation",height="640px")
        }
        else {
            require(reshape)
            if (!is.null(currentCorrelation$refGene)) {
                output$correlation <- renderUI(div())
                htmlOutput("correlation",height="640px")
                
                s <- currentMetadata$source
                d <- currentMetadata$dataset
                cc <- currentMetadata$final$class
                D <- currentCorrelation$datMatrix
                r <- currentCorrelation$refGene
                x <- D[r,]
                D <- D[setdiff(rownames(D),r),]
                if (nrow(D)>1000)
                    pwCorrPlot <- ggmessage(paste("Cannot create gene-wise",
                        "correlation plot with more than 1000 genes"),
                        type="error")
                else {
                    gNames <- as.character(loadedGenomes[[
                        currentMetadata$genome]]$dbGene[
                            rownames(D)]$gene_name)
                    rg <- as.character(loadedGenomes[[
                        currentMetadata$genome]]$dbGene[r]$gene_name)
                    rownames(D) <- gNames
                    melted <- melt(t(D))
                    
                    X <- rep(x,nrow(D))
                    Y <- as.numeric(melted$value)
                    X <- X/max(c(X,Y))
                    Y <- Y/max(c(X,Y))
                    
                    pwCorrData <- data.frame(
                        X=X,
                        Y=Y,
                        Gene=melted$X2,
                        Condition=rep(cc,nrow(D))
                    )
                    
                    if (length(levels(melted$X2))>6)
                        pwCorrPlot <- ggplot() +
                            geom_line(data=pwCorrData,aes(x=X,y=Y,
                                colour=Gene)) +
                            geom_point(data=pwCorrData,aes(x=X,y=Y,
                                colour=Gene,shape=Condition),size=2)
                            #,colour="#AEAEAE")
                    else
                        pwCorrPlot <- ggplot() +
                            geom_point(data=pwCorrData,aes(x=X,y=Y,
                                colour=Condition,fill=Condition,shape=Gene),
                                size=2) +
                            geom_line(data=pwCorrData,aes(x=X,y=Y,
                                linetype=Gene),colour="#AEAEAE")
                        
                    pwCorrPlot <- pwCorrPlot +
                        xlab(paste("Reference expression (",rg,")",sep="")) +
                        ylab("Expression of query genes\n") +
                        theme_bw() +
                        theme(
                            axis.title.x=element_text(size=12),
                            axis.title.y=element_text(size=12),
                            axis.text.x=element_text(size=10,face="bold"),
                            axis.text.y=element_text(size=10,face="bold"),
                            legend.position="bottom",
                            legend.title=element_text(size=11,face="bold"),
                            legend.text=element_text(size=10),
                            legend.key=element_blank()
                        )
                    output$correlation <- renderPlot({
                        pwCorrPlot
                    })
                    plotOutput("correlation",height="640px")
                }
            }
            else if (any(is.na(currentCorrelation$corMatrix))) {
                output$correlation <- renderPlot({
                    currentCorrelation$errorCor
                })
                plotOutput("correlation",height="640px")
            }
            else if (all(dim(currentCorrelation$corMatrix)>100)) {
                require(gplots)
                # Switch to static image
                cc <- unique(as.character(currentMetadata$final$class))
                classList <- vector("list",length(cc))
                names(classList) <- cc
                for (cl in cc)
                    classList[[cl]] <- 
                        as.character(
                            currentMetadata$final$sample_id[which(
                                currentMetadata$final$class==cl)])
                names(classList) <- baseColours[1:length(classList)]
                classes <- rep(names(classList),lengths(classList))
                output$correlation <- renderPlot({
                    heatmap.2(
                        currentCorrelation$corMatrix,
                        col=colorRampPalette(currentCorrelation$opts$colors),
                        revC=TRUE,
                        trace="none",
                        symm=currentCorrelation$opts$symm,
                        Colv=TRUE,
                        key=FALSE,
                        labRow="",
                        labCol="",
                        RowSideColors=classes,
                        ColSideColors=classes
                    )
                })
                #div(
                #    class="heatmap-container",
                    plotOutput("correlation",height="640px")
                #)
            }
            else {
                if (currentCorrelation$what=="genes") {
                    gg <- loadedGenomes[[currentMetadata$genome]]$geneNames
                    m <- match(rownames(currentCorrelation$corMatrix),gg)
                    labrow <- labcol <- names(gg[m])
                }
                else
                    labrow <- labcol <- rownames(currentCorrelation$corMatrix)
                output$correlation <- renderD3heatmap({
                    d3heatmap(
                        currentCorrelation$corMatrix,
                        dendrogram="both",
                        Rowv=TRUE,
                        Colv=TRUE,
                        colors=currentCorrelation$opts$colors,
                        revC=TRUE,
                        symm=currentCorrelation$opts$symm,
                        labRow=labrow,
                        labCol=labcol
                    )
                })
                #div(
                #    class="heatmap-container",
                    d3heatmapOutput("correlation",height="640px")
                #)
            }
        }
    })
    
    output$smallMds <- renderPlot({
        if (is.null(currentCorrelation$datMatrix))
            currentCorrelation$entryMds
        else {
            require(ggrepel)
            cc <- unique(as.character(currentMetadata$final$class))
            classList <- vector("list",length(cc))
            names(classList) <- cc
            for (cl in cc)
                classList[[cl]] <- 
                    as.character(
                        currentMetadata$final$sample_id[which(
                            currentMetadata$final$class==cl)])
            classes <- as.factor(rep(names(classList),lengths(classList)))
            x <- currentCorrelation$datMatrix
            d <- as.dist(0.5*(1-cor(x,method=currentCorrelation$opts$method)))
            tryCatch({
                mds.obj <- cmdscale(d,eig=TRUE,k=2)
                nd <- as.dist(0.5*(1-cor(t(mds.obj$points),
                    method=currentCorrelation$opts$method)))
                currentCorrelation$mdsRsq <- cor(c(d),c(nd))^2
                #gofx <- round(100*mds.obj$GOF[1],1)
                #gofy <- round(100*mds.obj$GOF[2],1)
                for.ggplot <- data.frame(
                    x=mds.obj$points[,1],
                    y=mds.obj$points[,2],
                    Condition=classes
                )
                if (nrow(for.ggplot)<=100)
                    rownames(for.ggplot) <- rownames(mds.obj$points)
                mds <- ggplot() +
                    geom_point(data=for.ggplot,mapping=aes(x=x,y=y,
                        colour=Condition,shape=Condition),size=2) +
                    #xlab(paste("\nPrincipal Coordinate 1 (",gofx,
                    #    "% goodness of fit)",sep="")) +
                    #ylab(paste("Principal Coordinate 2 (",gofy,
                    #    "% goodness of fit)\n",sep="")) +
                    xlab("\nPrincipal Coordinate 1") +
                    ylab("Principal Coordinate 2") +
                    theme_bw() +
                    theme(
                        axis.title.x=element_text(size=11),
                        axis.title.y=element_text(size=11),
                        axis.text.x=element_text(size=9,face="bold"),
                        axis.text.y=element_text(size=9,face="bold"),
                        legend.position="bottom",
                        legend.title=element_text(size=10,face="bold"),
                        legend.text=element_text(size=9),
                        legend.key=element_blank()
                    )
                if (nrow(for.ggplot)<=100)
                    mds <- mds +
                        geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
                            label=rownames(for.ggplot)),size=3)
                mds
            },
            warning=function(w) {
                currentCorrelation$mdsRsq <- NULL
                currentCorrelation$warnMds
            },
            error=function(e) {
                currentCorrelation$mdsRsq <- NULL
                currentCorrelation$errorMds
            },
            finally="")            
        }
    })
    
    output$smallMdsRsqDisplay <- renderText({
        if (is.null(currentCorrelation$mdsRsq))
            "R-square goodness of fit:"
        else
            paste("R-square goodness of fit:",
                round(currentCorrelation$mdsRsq,5))
    })
    
    output$rnaCorrelationMatrix <- renderUI({
        if (is.null(currentCorrelation$corMatrix))
            cormat <- data.frame(
                name=character(0),
                correlation=numeric(0)
            )
        else
            cormat <- round(currentCorrelation$corMatrix,digits=3)
        if (currentCorrelation$what=="genes") {
            gg <- loadedGenomes[[currentMetadata$genome]]$geneNames
            m <- match(rownames(currentCorrelation$corMatrix),gg)
            rownames(cormat) <- colnames(cormat) <- names(gg[m])
        }
        output$rnaCorCorMatrix <- 
            DT::renderDataTable(
                cormat,
                class="display compact",
                options=list(
                    searchHighlight=TRUE,
                    pageLength=10,
                    lengthMenu=c(10,20,50,100)
                )
            )
        list(
            div(
                class="small table-container",
                DT::dataTableOutput("rnaCorCorMatrix"),
                br(),
                div(
                    style="display:inline-block; margin:5px;",
                    downloadButton(
                        outputId="exportRnaCorrelationMatrix",
                        label="Export correlation matrix",
                        class="btn-xs"
                    )
                )
            )
        )
    })

    output$rnaCorDataMatrix <- DT::renderDataTable(
        if (is.null(currentCorrelation$datMatrix))
            data.frame(
                name=character(0),
                value=numeric(0)
            )
        else {
            s <- currentMetadata$source
            d <- currentMetadata$dataset
            gNames <- as.character(loadedGenomes[[
                currentMetadata$genome]]$dbGene[rownames(
                currentCorrelation$datMatrix)]$gene_name)
            data.frame(
                gene_name=gNames,
                currentCorrelation$datMatrix
            )
        },
        class="display compact",
        rownames=FALSE,
        options=list(
            searchHighlight=TRUE,
            pageLength=10,
            lengthMenu=c(10,20,50,100)
        )
    )
    
    output$exportRnaCorrelationMatrix <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("corr_matrix_",tt,".txt", sep='')
        },
        content=function(con) {
            cormat <- round(currentCorrelation$corMatrix,digits=5)
            if (currentCorrelation$what=="genes") {
                gg <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(rownames(currentCorrelation$corMatrix),gg)
                rownames(cormat) <- colnames(cormat) <- names(gg[m])
            }
            write.table(cormat,file=con,sep="\t",quote=FALSE,row.names=TRUE,
                col.names=NA)
        }
    )
    
    output$exportCorrelationPDF <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("rna_correlation_",tt,".pdf", sep='')
        },
        content=function(con) {
            require(gplots)
            corMat <- currentCorrelation$corMatrix
            pdf(file=con)
            n <- dim(corMat)[1]
            if (currentCorrelation$what=="genes") {
                gg <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(rownames(currentCorrelation$corMatrix),gg)
                labrow <- labcol <- names(gg[m])
            }
            else
                labrow <- labcol <- rownames(currentCorrelation$corMatrix)
            if (n<=100) {
                labs <- matrix(NA,n,n)
                for (i in 1:n)
                    for (j in 1:n)
                        labs[i,j] <- sprintf("%.2f",corMat[i,j])
                if (n <= 5)
                    notecex <- 1.2
                else if (n > 5 & n < 10)
                    notecex <- 0.9
                else
                    notecex <- 0.7
                heatmap.2(
                    corMat,
                    col=colorRampPalette(currentCorrelation$opts$colors),
                    revC=TRUE,
                    trace="none",
                    symm=currentCorrelation$opts$symm,
                    Colv=TRUE,
                    cellnote=labs,
                    keysize=1,
                    density.info="density",
                    notecex=notecex,
                    cexCol=0.9,
                    cexRow=0.9,
                    font.lab=2,
                    labRow=labrow,
                    labCol=labcol
                )
            }
            else {
                heatmap.2(
                    corMat,
                    col=colorRampPalette(currentCorrelation$opts$colors),
                    revC=TRUE,
                    trace="none",
                    symm=currentCorrelation$opts$symm,
                    Colv=TRUE,
                    keysize=1,
                    density.info="density",
                    cexCol=0.9,
                    cexRow=0.9,
                    font.lab=2
                )
            }
            dev.off()
        }
    )
    
    output$exportCorrelationPNG <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("rna_correlation_",tt,".png", sep='')
        },
        content=function(con) {
            require(gplots)
            corMat <- currentCorrelation$corMatrix
            png(filename=con,width=800,height=800)
            n <- dim(corMat)[1]
            if (currentCorrelation$what=="genes") {
                gg <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(rownames(currentCorrelation$corMatrix),gg)
                labrow <- labcol <- names(gg[m])
            }
            else
                labrow <- labcol <- rownames(currentCorrelation$corMatrix)
            if (n<=100) {
                labs <- matrix(NA,n,n)
                for (i in 1:n)
                    for (j in 1:n)
                        labs[i,j] <- sprintf("%.2f",corMat[i,j])
                if (n <= 5)
                    notecex <- 1.2
                else if (n > 5 & n < 10)
                    notecex <- 0.9
                else
                    notecex <- 0.7
                heatmap.2(
                    corMat,
                    col=colorRampPalette(currentCorrelation$opts$colors),
                    revC=TRUE,
                    trace="none",
                    symm=currentCorrelation$opts$symm,
                    Colv=TRUE,
                    cellnote=labs,
                    keysize=1,
                    density.info="density",
                    notecex=notecex,
                    cexCol=0.9,
                    cexRow=0.9,
                    font.lab=2,
                    labRow=labrow,
                    labCol=labcol
                )
            }
            else {
                heatmap.2(
                    corMat,
                    col=colorRampPalette(currentCorrelation$opts$colors),
                    revC=TRUE,
                    trace="none",
                    symm=currentCorrelation$opts$symm,
                    Colv=TRUE,
                    keysize=1,
                    density.info="density",
                    cexCol=0.9,
                    cexRow=0.9,
                    font.lab=2
                )
            }
            dev.off()
        }
    )
}

correlationTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCorrelation <- allReactiveVars$currentHeatmap
    
    correlationTabPanelReactiveEvents <- 
        correlationTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    performRnaCorrelation <- 
        correlationTabPanelReactiveEvents$performRnaCorrelation
    
    correlationTabPanelReactiveExprs <- 
        correlationTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    correlationTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
        
    observe({
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"selectCorrelationGeneName",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectCorrelationGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectCorrelationGeneName",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"rnaCorrelationRefGene",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectCorrelationGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"rnaCorrelationRefGene",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (input$rnaCorrelateWhat=="refgene" 
            && isEmpty(input$rnaCorrelationRefGene))
            shinyjs::disable("performRnaCorrelation")
        else
            shinyjs::enable("performRnaCorrelation")
    })
    
    observe({
        tryCatch({
            shinyjs::disable("performRnaCorrelation")
            performRnaCorrelation()
        },error=function(e) {
            #print(e)
        },
        finally={
            shinyjs::enable("performRnaCorrelation")
        })
    })
}

