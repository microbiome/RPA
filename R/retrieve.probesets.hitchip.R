#' @title Retrieve probesets 
#' @description List probes for each probeset in taxonomic data.
#' @param tax.table data.frame with oligo - phylotype mapping info
#' @param level phylotype level for probesets
#' @param name specify phylotypes to check (optional)
#' @return A list. Probes for each phylotype.
#' @examples 
#'   #tax.table <- GetPhylogeny('HITChip')
#'   #sets <- retrieve.probesets(tax.table, 'species', 'Weissella confusa')
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
retrieve.probesets <- function(tax.table, level = "species", name = NULL) {

    tax.table <- as.data.frame(tax.table)

    # If name not given, pick all
    if (is.null(name)) {
        name <- unique(as.character(tax.table[[level]]))
    }
    
    phylo <- tax.table[tax.table[[level]] %in% name, ]
    
    if (is.factor(phylo[[level]])) {
        phylo[[level]] <- droplevels(phylo[[level]])
    }
    
    phylo.list <- split(phylo, phylo[[level]])
    probesets <- lapply(phylo.list, function(x) {
        as.character(unique(x$oligoID))
    })
    names(probesets) <- names(phylo.list)
    
    probesets
    
} 

