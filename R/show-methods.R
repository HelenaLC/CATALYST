# ================================================================================
# Show method for dbFrame
# --------------------------------------------------------------------------------
setMethod(f="show",
    signature="dbFrame",
    definition=function(object) {
        dims <- dim(object@exprs)
        n_bcs <- nrow(object@bc_key)
        cat("dbFrame objectect with\n", 
            dims[1], "events,", 
            dims[2], "observables and",
            n_bcs, "barcodes:\n")
        
        cat("\nCurrent assignments:\n")
        cat("     ", sum(object@bc_ids == 0), "events unassigned")

        tmp <- c(sort(table(object@bc_ids[object@bc_ids != 0]), decreasing=TRUE))
        tbl <- t(data.frame("ID"=as.numeric(names(tmp)), "Count"=tmp))
        dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
        print(noquote(tbl))
        
        if (length(object@sep_cutoffs) != 0) {
            cat("\nSeparation cutoffs:")
            tbl <- t(data.frame(
                "ID"=rownames(object@bc_key), 
                "Yield"=object@sep_cutoffs))
            dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
            print(noquote(tbl))
            
            seps <- seq(0, 1, .01)
            inds <- sapply(1:n_bcs, function(x) 
                which(seps == object@sep_cutoffs[x]))
            yields <- object@yields[cbind(1:n_bcs, inds)]
            cat("\nYields upon debarcoding:\n")
            cat("     ", paste0(round(mean(yields), 4) * 100, "%"), "overall yield")
            
            tbl <- t(data.frame(
                "ID"=rownames(object@bc_key), 
                "Yield"=paste0(round(yields, 4) * 100, "%")))
            dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
            print(noquote(tbl))
        }
    })
# --------------------------------------------------------------------------------


