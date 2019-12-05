sessionLoaderTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    
    bookmarkState <- eventReactive(input$bookmarkBtn,{
        # If no dataset has been created, do not allow session creation
        if (is.null(currentMetadata$final)) {
            showModal(modalDialog("Please create a dataset first! SeqCVIBE ",
                "sessions can be created only if there is an active dataset!",
                title="No dataset!",
                easyClose=TRUE,
                size="s"
            ))
        }
        else {
            session$doBookmark()
            showNotification(
                paste(
                    getTime("SUCCESS:"),
                    "Bookmark: '",
                    input$description,
                    "' created",sep=""
                ),
                duration=5, 
                type='message')
        }
    })
    
    delSt <- reactiveVal()
    delSe <- reactiveVal()
    preDeleteSession <- eventReactive(input$deleteBM, {
        str <- strsplit(input$deleteBM, "_")[[1]]
        delSt(str[2])
        delSe(str[3])
        showModal(modalDialog(
            title="Confirm session delete!",
            "The selected session is about to be deleted. Are you sure? This ",
            "action cannot be undone!",
            easyClose=FALSE,
            footer=tagList(
                modalButton("Cancel",icon=icon("ban")),
                actionButton("confirmSessionDelete","Delete",class="btn-danger",
                    icon=icon("exclamation-triangle"))
            )
        ))
    })
    
    deleteBookmark <- eventReactive(input$confirmSessionDelete,{
        st <- as.character(delSt())
        se <- as.character(delSe())
        nr <- dbExecute(metadata,paste0("DELETE FROM bookmarks ",
            "WHERE state_id='",st,"' AND session='",se,"'"))
        #print(paste0("DELETE FROM bookmarks ",
        #   "WHERE state_id='",st,"' AND session='",se,"'"))
        if (nr > 0)
            showModal(modalDialog(
                title="Session deleted!",
                "The selected session and all related analyses have been ",
                    "succesfully deleted!",
                easyClose=TRUE,
                footer=tagList(
                    modalButton("Dismiss",icon=icon("check"))
                )
            ))
    })
    
    return(list(
        bookmarkState=bookmarkState,
        preDeleteSession=preDeleteSession,
        deleteBookmark=deleteBookmark
    ))
}

sessionLoaderTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    
    sessionInfoTablePoll <- reactivePoll(5000,session,checkFunc=function() {        
        rowcount <- dbGetQuery(metadata,"SELECT Count(*) FROM bookmarks")[1,1]
        return(rowcount)
    },valueFunc=function() {
        uid <- as.numeric(USER_ID())
        if (length(uid) == 0) # USER_ID is NULL
            qq <- paste0("SELECT description, state_id, timestamp, ",
                "session FROM bookmarks WHERE user_id IS NULL")
        else
            qq <- paste0("SELECT description, state_id, timestamp, ",
                "session FROM bookmarks WHERE user_id=",uid)
        urlDF <- dbGetQuery(metadata,qq)
        currentMetadata$urlDF <- urlDF
        return(urlDF)
    })
    
    sessionNameNotEmpty <- reactive({
        if (isEmpty(input$description))
            shinyjs::disable("bookmarkBtn")
        else
            shinyjs::enable("bookmarkBtn")
    })
    
    return(list(
        sessionInfoTablePoll=sessionInfoTablePoll,
        sessionNameNotEmpty=sessionNameNotEmpty
    ))
}

sessionLoaderTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    # Bookmarks table
    currentMetadata <- allReactiveVars$currentMetadata
    
    output$urlTable = renderDT({
        if (dbExistsTable(metadata,"bookmarks")) {
            uid <- as.numeric(USER_ID())
            
            if (length(uid) == 0)
                qq <- paste0("SELECT description, state_id, timestamp, ",
                    "session FROM bookmarks WHERE user_id IS NULL")
            else
                qq <- paste0("SELECT description, state_id, timestamp, ",
                    "session FROM bookmarks WHERE user_id=",uid)
            tmpDf <- dbGetQuery(metadata,qq)
            
            if (nrow(tmpDf) > 0) {
                currentMetadata$urlDF <- tmpDf
                tmpDf$timestamp <- as.POSIXct(tmpDf$timestamp,
                    origin="1970-01-01 00:00")
                 
                base <- paste0(
                    session$clientData$url_protocol,"//",
                    session$clientData$url_hostname,
                    if (!isEmpty(session$clientData$url_port))
                        paste0(":",session$clientData$url_port)
                    else "",
                    session$clientData$url_pathname,"?_state_id_="
                )
                    
                tmpDf$delete <- shinyInput(actionButton,nrow(tmpDf),
                    'deleteSession_',paste0(tmpDf$state_id,"_",tmpDf$session),
                    label="Delete",icon=icon("minus-circle"),
                    class="btn-danger btn-xs",
                    onclick='Shiny.onInputChange(\"deleteBM\",this.id)')
                tmpDf$state_id <- paste0("<a href='",base,tmpDf$state_id,
                    "'>Restore</a>")
                tmpDf <- tmpDf[,c(1:3,5)]
                names(tmpDf)[2:3] <- c("link","date created")
                datatable(tmpDf,
                    rownames=FALSE,
                    class="display",
                    filter="top",
                    escape=FALSE,
                    selection=list(
                        mode="single",
                        target = 'cell'
                    )
                ) %>% formatStyle(1,cursor='alias') %>% 
                formatDate(3,method='toLocaleString')
            }            
        }
    })
}

sessionLoaderTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    
    sessionLoaderTabPanelEventReactiveEvents <- 
        sessionLoaderTabPanelEventReactive(input,output,session,allReactiveVars,
            allReactiveMsgs)
    
    bookmarkState <- sessionLoaderTabPanelEventReactiveEvents$bookmarkState
    preDeleteSession <- 
        sessionLoaderTabPanelEventReactiveEvents$preDeleteSession
    deleteBookmark <- sessionLoaderTabPanelEventReactiveEvents$deleteBookmark
    
    sessionLoaderTabPanelReactiveExprs <- sessionLoaderTabPanelReactive(input,
        output,session,allReactiveVars,allReactiveMsgs)
    
    sessionInfoTablePoll <- 
        sessionLoaderTabPanelReactiveExprs$sessionInfoTablePoll
    sessionNameNotEmpty <- 
        sessionLoaderTabPanelReactiveExprs$sessionNameNotEmpty
    
    sessionLoaderTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
        
    observe({
        bookmarkState()
    })
    observe({
        preDeleteSession()
        deleteBookmark()
    })
    observe({
        sessionNameNotEmpty()
    })
    observe({
        sessionInfoTablePoll()
    })
    
    # We need to store only the _state_id_ and then construct a link 
    # independently - states stored on server only!
    onBookmarked(fun=function(url) {
        stateId <- strsplit(basename(url),"=")[[1]][2]
        
        # Check if stateId already exists
        qCheck <- paste0("SELECT COUNT(1) FROM bookmarks WHERE state_id='",
            stateId,"'")
        n <- dbGetQuery(metadata,qCheck)[1,1]
        
        if (n == 0) {
            uid <- as.numeric(USER_ID())
            if (length(uid) == 0)
                iQuery <- paste0("INSERT INTO bookmarks (description, ",
                    "state_id, timestamp, session) VALUES ('",input$description,
                    "',","'",stateId,"',",as.numeric(Sys.time()),",'",
                    session$token,"')")
            else
                iQuery <- paste0("INSERT INTO bookmarks (description, ",
                    "state_id, timestamp, session, user_id) VALUES ('",
                    input$description,"',","'",stateId,"',",
                    as.numeric(Sys.time()),",'",session$token,"',",uid,")")
            print(iQuery)
            nr <- dbExecute(metadata,iQuery)
        }
        else {
            showModal(modalDialog(
                title="Session exists?",
                "The session you are trying to save was found in the database.",
                " Are you sure you have performed/changed an analysis since ",
                "you last saved the dataset actions?",
                easyClose=TRUE,
                footer=tagList(
                    modalButton("Dismiss",icon=icon("check"))
                )
            ))
        }
    })
}

shinyInput <- function(FUN,len,idPrefix,idSuffix,...) {
    inputs <- character(len)
    if (missing(idSuffix))
        idSuffix <- 1:len
    counter <- 0
    for (i in idSuffix) {
        counter <- counter + 1
        inputs[counter] <- as.character(FUN(paste0(idPrefix,i),...))
    }
        
    return(inputs)
}
