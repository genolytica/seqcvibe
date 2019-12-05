initDatabase <- function(db) {
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to build the annotation database!")
    if (missing(db))
        stop("A database file must be provided!")
    drv <- dbDriver("SQLite")
    if (file.exists(db))
        # The database has been created at least with the tables defined
        con <- dbConnect(drv,dbname=db)
    else {
        # Create database file and define tables
        con <- dbConnect(drv,dbname=db)
        rs <- .initTables(con)
    }
    return(con)
}

.initTables <- function(con) {
    queries <- .localTblDef()
    rs <- dbSendQuery(con,queries[[1]])
    if (dbHasCompleted(rs))
        dbClearResult(rs)
    for (n in names(queries)) {
        rs <- dbSendStatement(con,queries[[n]])
        if (dbHasCompleted(rs))
            dbClearResult(rs)
    }
}

.localTblDef <- function() {
    return(list(
        enable_fkey="PRAGMA foreign_keys=1;",
        metadata=paste0(
            "CREATE TABLE IF NOT EXISTS metadata (",
            "_id INTEGER PRIMARY KEY AUTOINCREMENT,",
            "sample_id TEXT,",
            "dataset TEXT,",
            "class TEXT,",
            "source TEXT,",
            "sample_dir TEXT,",
            "track_dir TEXT,",
            "alt_id TEXT,",
            "library_strategy TEXT,",
            "quality INTEGER,",
            "norm_factor REAL,",
            "genome TEXT,",
            "is_public INTEGER DEFAULT 0,",
            "user_id INTEGER DEFAULT NULL,",
            "dataset_id INTEGER DEFAULT NULL,",
            "FOREIGN KEY(user_id) REFERENCES users(_id) ",
            "ON DELETE CASCADE,",
            "FOREIGN KEY(dataset_id) REFERENCES datasets(_id) ",
            "ON DELETE CASCADE",
            ");"
        ),
        datasets=paste0(
            "CREATE TABLE IF NOT EXISTS datasets (",
            "_id INTEGER PRIMARY KEY AUTOINCREMENT,",
            "dataset TEXT,",
            "title TEXT,",
            "link TEXT,",
            "short_summary TEXT,",
            "user_id INTEGER DEFAULT NULL,",
            "FOREIGN KEY(user_id) REFERENCES users(_id) ",
            "ON DELETE CASCADE",
            ");"
        ),
        bookmarks=paste(
            "CREATE TABLE IF NOT EXISTS bookmarks (",
            "_id INTEGER PRIMARY KEY AUTOINCREMENT,",
            "description TEXT,",
            "state_id TEXT NOT NULL,",
            "timestamp REAL,",
            "session TEXT,",
            "user_id INTEGER DEFAULT NULL,",
            "FOREIGN KEY(user_id) REFERENCES users(_id) ",
            "ON DELETE CASCADE",
            ");"
        ),
        users=paste(
            "CREATE TABLE IF NOT EXISTS users (",
            "_id INTEGER PRIMARY KEY AUTOINCREMENT,",
            "email TEXT NOT NULL UNIQUE,",
            "name TEXT,",
            "role INTEGER NOT NULL DEFAULT 0",
            ");"
        )
        # On login, the user will be registered in the local db, if not exist
    ))
}

#~ library(RSQLite)
#~ con <- dbConnect(dbDriver("SQLite"),dbname="metadata.sqlite")
#~ tab <- dbGetQuery(con,"SELECT * FROM metadata")

#~ sample_dir <- gsub("/home/makis/elixir-RNAseq/data/","/media/storage/pilot/",
#~  tab$sample_dir)
#~ track_dir <- gsub("/home/makis/elixir-RNAseq/data/","/media/storage/pilot/",
#~  tab$track_dir)

#~ tab$sample_dir <- sample_dir
#~ tab$track_dir <- track_dir

#~ dbWriteTable(con,"metadata",tab,overwrite=TRUE)

#~ dbDisconnect(con)
