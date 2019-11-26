sessionLoaderTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

sessionLoaderTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

sessionLoaderTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    # Bookmarks table
    currentMetadata <- allReactiveVars$currentMetadata
    output$urlTable = DT::renderDataTable({
        req(currentMetadata$urlDF)
        currentMetadata$urlDF
    }, escape = FALSE)
}

sessionLoaderTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    
    sessionLoaderTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    #if (dbExistsTable(bookmarks,"Bookmarks")) {
    if (dbExistsTable(metadata,"bookmarks")) {
        tmpDf <- dbGetQuery(metadata,paste0("SELECT description,url,
            timestamp,session FROM bookmarks"))
        # Placeholder
        #tmpUrlDF <- dbGetQuery(metadata,paste0("SELECT description,url,
        #   timestamp,session FROM bookmarks WHERE user='",USER,'"))
        tmpUrlDF <- data.table(tmpDf)
        #tmpUrlDF <- data.frame(
        #    description=paste('<a href="',tmpDf$url,'>',tmpDf$description,
        #        '</a>',sep=""),
        #    url=as.POSIXct(timestamp,origin="1970-01-01 00:00")
        #)
        currentMetadata$urlDF <- tmpUrlDF[,timestamp := as.POSIXct(timestamp,
            origin="1970-01-01 00:00")]
        #currentMetadata$urlDF <- tmpUrlDF
    } else
        currentMetadata$urlDF <- NULL

    observeEvent(input$bookmarkBtn, {
        session$doBookmark()
        showNotification(
            paste(
                getTime("SUCCESS:"),
                "Bookmark: '",
                input$description,
                "' created",sep=""
            ),  
            duration = 5, 
            type = 'message')
        #sessionLoaderMessages <- updateMessages(
        #    sessionLoaderMessages,
        #    type="SUCCESS",
        #    msg=paste(getTime("SUCCESS"),
        #        "Bookmark ",input$description," created",sep="")
        #)
    })

    observeEvent(input$deleteBM,{
        if (!is.null(input$urlTable_rows_selected)) {
            showNotification(
                paste(
                    getTime("SUCCESS:"),
                    "Bookmark(s): '",
                    paste(
                        currentMetadata$urlDF$Description[as.numeric(
                            input$urlTable_rows_selected)],
                        collapse=", "
                    ),
                    "' removed",sep=""
                ),  
                duration = 5, 
                type = 'message')
            #sessionLoaderMessages <- updateMessages(
            #     sessionLoaderMessages,
            #     type="SUCCESS",
            #     msg=paste(getTime("SUCCESS"),
            #       "Bookmark(s) ",
            #       paste(currentMetadata$urlDF$Description[as.numeric(
            #           input$urlTable_rows_selected)],collapse=", "),
            #        " removed",sep="")
            #)
            currentMetadata$urlDF <- currentMetadata$urlDF[-as.numeric(
                input$urlTable_rows_selected),]
        }
    })

    observeEvent(input$urlTable_rows_selected,{
        if (!is.null(input$urlTable_rows_selected))
            shinyjs::enable("deleteBM")
        else
            shinyjs::disable("deleteBM")
    })

    session$onSessionEnded(function() {
        tmpUrlDF <- isolate({currentMetadata$urlDF})
        if (!is.null(tmpUrlDF)) {
            ## Delete previous sessions by user in order to massively update
            #dbExecute(metadata,"DELETE FROM bookmarks WHERE user='","'")
            ## Rewrite
            #tmpUrlDF$user <- USER
            #dbWriteTable(metadata,"bookmarks",tmpUrlDF,overwrite=TRUE)
            tmpUrlDF$user <- "user"
            dbWriteTable(metadata,"bookmarks",tmpUrlDF,overwrite=TRUE)
            
            #dbWriteTable(bookmarks,"Bookmarks",tmpUrlDF,overwrite=TRUE)
        }
        #dbDisconnect(metadata)
    })

    onBookmarked(fun=function(url) {
        #if (!url %in% currentMetadata$urlDF$URL) {
        if (!url %in% currentMetadata$urlDF$url) {
            if (is.null(currentMetadata$urlDF)) {
                #currentMetadata$urlDF <- unique(data.table(
                #    Description=input$description,
                #    URL=paste0("<a href='", url, "'>", url, "</a>"),
                #    Timestamp=Sys.time(),
                #    Session=session$token
                #),by="URL")
                currentMetadata$urlDF <- unique(data.table(
                    description=input$description,
                    url=paste0("<a href='", url, "'>", url, "</a>"),
                    timestamp=Sys.time(),
                    session=session$token
                ),by="url")
            } else {
                #currentMetadata$urlDF<-unique(rbindlist(list(
                #    currentMetadata$urlDF,data.table(
                #        description=input$description,
                #        url=paste0("<a href='", url, "'>", url, "</a>"),
                #        timestamp=Sys.time(),
                #        session=session$token
                #))),by="URL")
                currentMetadata$urlDF<-unique(rbindlist(list(
                    currentMetadata$urlDF,data.table(
                        description=input$description,
                        url=paste0("<a href='", url, "'>", url, "</a>"),
                        timestamp=Sys.time(),
                        session=session$token
                ))),by="url")
            }
        }
    })
}
