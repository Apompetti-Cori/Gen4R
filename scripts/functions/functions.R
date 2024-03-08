reorder_gt <- function(input) {

    if(class(input) != "character"){
        stop("input is not a character vector.")
    }

    # Function to sort a character vector containing "0/1" or "1/0" formats.
    reordered_data <- lapply(input, function(item) {
        # Split the string by "/"
        parts <- strsplit(item, split = "/")[[1]]  # Extract the first element (vector)
        # Sort the parts numerically
        sorted_parts <- sort(as.numeric(parts))
        # Join the sorted parts back with "/"
        sorted_parts <- paste(sorted_parts, collapse = "/")
        if(sorted_parts == ""){
            sorted_parts <- NA
        }

        return(sorted_parts)
    })
        return(unlist(reordered_data))
}
