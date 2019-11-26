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
            "user TEXT",
            ");"
        ),
        summaries=paste0(
            "CREATE TABLE IF NOT EXISTS summaries (",
            "_id INTEGER PRIMARY KEY AUTOINCREMENT,",
            "dataset TEXT NOT NULL,",
            "title TEXT,",
            "link TEXT,",
            "short_summary TEXT,",
            "content_id INTEGER NOT NULL,",
            "FOREIGN KEY(dataset) REFERENCES metadata(dataset) ",
            "ON DELETE CASCADE",
            ");"
        ),
        bookmarks=paste(
            "CREATE TABLE IF NOT EXISTS bookmarks (",
            "_id INTEGER PRIMARY KEY AUTOINCREMENT,",
            "description TEXT,",
            "url TEXT,",
            "timestamp REAL,",
            "session TEXT,",
            "user TEXT",
            ");"
        )
    ))
}
