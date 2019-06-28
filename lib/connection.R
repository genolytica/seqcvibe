initSdb <- function() {
    db <- new.env(parent=emptyenv())
    db$connection <- NULL
    db$open <- function() {
        if (!require(RSQLite))
            stop("R package RSQLite is required!")
        con <- tryCatch({
            drv <- dbDriver("SQLite")
			con <- dbConnect(drv,dbname=dbFile)
		},error=function(e) {
			stop(e)
		},finally="")
        db$connection <- con
    }
    db$close <- function() {
        if (!require(RSQLite))
            stop("R package RSQLite is required!")
        invisible(dbDisconnect(db$connection))
    }
    db$insert <- function() {
    }
    db$update <- function() {
    }
    db$find <- function() {
    }
    db$status <- function() {
        if(dbIsValid(db$connection))
            return("Connected!")
        else
            return("Disconnected!")
    }
    db$getConfig <- function() {
    }
    db$getSubConfig <- function(s,d,c) {
    }
    db$getSources <- function() {
    }
    db$getDatasets <- function(s) {
    }
    db$getClasses <- function(s,d) {
    }
    db$getSamples <- function(s,d,c,type=c("prim","alt","both")) {
    }
    db$getSampleTables <- function(s,d,c) {
    }
    db$getNormFactors <- function(s) {
    }
    db$getExpression(s,g,type=c("raw","norm") {
		type=tolower(type[1])
		switch(type,
			raw = {
			},
			norm = {
			}
		)
		SELECT 
	}
    return(db)
}

initMdb <- function(creds) {
    if (is.null(creds$host) || is.na(creds$host) || creds$host == "") {
        #warning("MongoDB host not given! Assuming default (127.0.0.1)...",
        #    immediate.=TRUE)
        creds$host <- "127.0.0.1"
    }
    if (is.null(creds$port) || is.na(creds$port) || creds$port == "") {
        #warning("MongoDB port not given! Assuming default (27017)...",
        #    immediate.=TRUE)
        creds$port <- "27017"
    }
    if (is.null(creds$username) || is.na(creds$username) 
        || creds$username == "") {
        #warning("MongoDB username not given! Assuming no username...",
        #    immediate.=TRUE)
        creds$username <- ""
    }
    if (is.null(creds$password) || is.na(creds$password) 
        || creds$password == "") {
        #warning("MongoDB password not given! Assuming no password...",
        #    immediate.=TRUE)
        creds$password <- ""
    }
    if (is.null(creds$db) || is.na(creds$db) || creds$db == "") {
        #warning("MongoDB database! Assuming default (test)...",
        #    immediate.=TRUE)
        creds$db <- "test"
    }

    mdb <- new.env(parent=emptyenv())
    mdb$connection <- NULL
    mdb$credentials <- creds
    mdb$open <- function() {
        if (!require(rmongodb))
            stop("R package rmongodb is required to connect to MongoDB ",
                "instance!")
        host <- ifelse(is.null(mdb$credentials$host),"127.0.0.1",
            mdb$credentials$host)
        port <- ifelse(is.null(mdb$credentials$port),"27017",
            mdb$credentials$port)
        username <- ifelse(is.null(mdb$credentials$username),"",
            mdb$credentials$username)
        password <- ifelse(is.null(mdb$credentials$password),"",
            mdb$credentials$password)
        db <- mdb$credentials$db
        if (is.null(db) || db == "")
            stop("A valid MongoDB database name must be provided!")
        host <- paste(host,port,sep=":")
        mdb.con <- tryCatch(
            mongo.create(db=db,host=host,username=username,password=password),
            error=function(e) {
                mongo.destroy(mdb.con)
                stop(e)
            },finally="")
        mdb$connection <- mdb.con
    }
    mdb$close <- function() {
        if (!require(rmongodb))
            stop("R package rmongodb is required to disconnect from ",
                "MongoDB!")
        invisible(mongo.destroy(mdb$connection))
    }
    mdb$insert <- function(collection,dat,mod="one") {
        if (!(mod %in% c("one","many")))
            stop("Argument mod must be one of \"one\" or \"many\"")
        if (!is.list(dat) && !is.character(dat))
            stop("Data to import (dat) must be either a list or a JSON ",
                "character vector!")
        if(mongo.is.connected(mdb$connection)) {
            ns <- paste(mdb$credentials$db,collection,sep=".")
            if (mod == "one")
                mongo.insert(mdb$connection,ns,dat)
            else if (mod == "many") {
                for (i in 1:length(dat))
                    mongo.insert(mdb$connection,ns,dat[[i]])
            }
            return(TRUE)
        }
        else {
            message("No connection to MongoDB!")
            return(FALSE)
        }
    }
    mdb$update <- function(collection,crit,dat,mod="one") {
        if (!(mod %in% c("one","many")))
            stop("Argument mod must be one of \"one\" or \"many\"")
        if (!is.list(dat) && !is.character(dat))
            stop("Data to import (dat) must be either a list or a JSON ",
                "character vector!")
        if(mongo.is.connected(mdb$connection)) {
            ns <- paste(mdb$credentials$db,collection,sep=".")
            if (mod == "one")
                mongo.update(mdb$connection,ns,crit,dat)
            else if (mod == "many") {
                for (i in 1:length(dat))
                    mongo.update(mdb$connection,ns,crit,dat[[i]])
            }
            return(TRUE)
        }
        else {
            message("No connection to MongoDB!")
            return(FALSE)
        }
    }
    mdb$find <- function(collection,query=mongo.bson.empty(),
        sort=mongo.bson.empty(),fields=mongo.bson.empty(),limit=0L,skip=0L,
        options=0L,data.frame=FALSE) {
        if (!is.list(query) && !is(query,"mongo.bson"))
            stop("The query must be an R list or a BSON list!")
        if(mongo.is.connected(mdb$connection)) {
            if (data.frame)
                oid2char <- TRUE
            else
                oid2char <- FALSE
            ns <- paste(mdb$credentials$db,collection,sep=".")
            result <- mongo.find.all(mdb$connection,ns,query=query,
                sort=sort,fields=fields,limit=limit,skip=skip,
                options=options,data.frame=data.frame,
                mongo.oid2character=oid2char)
            return(result)
        }
        else {
            message("No connection to MongoDB!")
            return(FALSE)
        }
    }
    mdb$status <- function() {
        if(mongo.is.connected(mdb$connection))
            return("Connected!")
        else
            return("Disconnected!")
    }
    mdb$getConfig <- function() {
        return(mdb$find("metadata",data.frame=TRUE))
    }
    mdb$getSubConfig <- function(s,d,c) {
        query <- list(source=s,dataset=d,class=c)
        return(mdb$find("metadata",query=query,data.frame=TRUE))
    }
    mdb$getSources <- function() {
        if (mongo.is.connected(mdb$connection)) {
            ns <- paste(mdb$credentials$db,"metadata",sep=".")
            return(mongo.distinct(mdb$connection,ns,"source"))
        }
        else {
            message("No connection to MongoDB!")
            return(FALSE)
        }
    }
    mdb$getDatasets <- function(s) {
        query <- list(source=s)
        fields <- list(dataset=1L)
        res <- mdb$find("metadata",query=query,fields=fields)
        return(unique(sapply(res,function(x) return(x$dataset))))
    }
    mdb$getClasses <- function(s,d) {
        query <- list(source=s,dataset=d)
        fields <- list(class=1L)
        res <- mdb$find("metadata",query=query,fields=fields)
        return(unique(sapply(res,function(x) return(x$class))))
    }
    mdb$getSamples <- function(s,d,c,type=c("prim","alt","both")) {
        query <- list(source=s,dataset=d,class=c)
        fields <- list(sample_id=1L,alt_id=1L)
        res <- mdb$find("metadata",query=query,fields=fields)
        sample_ids <- unlist(sapply(res,function(x) return(x$sample_id)))
        alt_ids <- unlist(sapply(res,function(x) {
            if (is.null(x$alt_id))
                return("N/A")
            else
                return(x$alt_id)
        }))
        na <- which(alt_ids=="N/A")
        if (length(na)>0)
            alt_ids[na] <- sample_ids[na]
        names(sample_ids) <- alt_ids
        switch(type,
            prim = {
                return(unlist(sample_ids))
            },
            alt = {
                return(unlist(alt_ids))
            },
            both = {
                the_ids <- unlist(sample_ids)
                names(the_ids) <- paste(names(the_ids)," (",the_ids,")",sep="")
                return(the_ids)
            }
        )
    }
    mdb$getSampleTables <- function(s,d,c) {
        query <- list(source=s,dataset=d,class=c)
        fields <- list(sample_id=1L,alt_id=1L,norm_factor=1L,
            library_strategy=1L,quality=1L)
        return(mdb$find("metadata",query=query,fields=fields,data.frame=TRUE))
    }
    mdb$getLengths <- function(id,chr) {
        query <- list(
            sample_id=id,
            seq=chr
        )
        fields <- list(value=1L)
        res <- mdb$find("runlength",query=query,fields=fields)
        return(as.integer(strsplit(unlist(res[[1]]$value),split="")[[1]]))
    }
    mdb$getValues <- function(id,chr) {
        query <- list(
            sample_id=id,
            seq=chr
        )
        fields <- list(value=1L)
        res <- mdb$find("runvalue",query=query,fields=fields)
        return(as.integer(strsplit(unlist(res[[1]]$value),split="")[[1]]))
    }
    mdb$getNormFactors <- function(s) {
        query <- list(
            "sample_id"=list(
                "$in"=s
            )
        )
        fields <- list(sample_id=1L,norm_factor=1L)
        res <- mdb$find("metadata",query=query,fields=fields)
        normFactors <- unlist(sapply(res,function(x) return(x$norm_factor)))
        names(normFactors) <- unlist(sapply(res,
            function(x) return(x$sample_id)))
        return(normFactors)
    }
    return(mdb)
}


