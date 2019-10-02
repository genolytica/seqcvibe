docTabPanel <- function() {
    fluidRow(column(12,
        tags$iframe(
            src="manual.html",
            width="100%",height=860,
            frameborder=0
        )
    ))
}
