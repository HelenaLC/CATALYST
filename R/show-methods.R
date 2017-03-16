# ==============================================================================
# Show method for dbFrame
# ------------------------------------------------------------------------------

setMethod(f="show",
    signature="dbFrame",
    definition=function(object) {
        dims <- dim(object@exprs)
        cat("dbFrame objectect with\n", 
            dims[1], " events, ", 
            dims[2], " observables and ",
            nrow(object@bc_key), " barcodes:\n", sep="")
        
        cat("\nCurrent assignments:\n")
        cat("     ", sum(object@bc_ids == 0), "event(s) unassigned")
        
        tmp <- c(sort(table(object@bc_ids[object@bc_ids != 0]), 
            decreasing=TRUE))
        tbl <- t(data.frame("ID"=names(tmp), "Count"=tmp))
        dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
        print(noquote(tbl))
        
        if (length(object@sep_cutoffs) != 0) {
            inds <- match(names(tmp), rownames(object@bc_key))
            cat("\nSeparation cutoffs:")
            tbl <- t(data.frame(
                "ID"=names(tmp), 
                "Cutoff"=object@sep_cutoffs[inds]))
            dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
            print(noquote(tbl))
            
            # better way to do this?
            seps <- seq(0, 1, .01)
            ind <- NULL
            for (i in object@sep_cutoffs[inds]) {
                for (j in seps) {
                    if (isTRUE(all.equal(i, j))) 
                        ind <- append(ind, which(seps == j))
                }
            }
            yields <- object@yields[cbind(inds, ind)]
            cat("\nYields upon debarcoding:\n")
            cat("     ", paste0(round(mean(yields), 4) * 100, "%"), 
                "overall yield")
            
            tbl <- t(data.frame(
                "ID"=names(tmp), 
                "Yield"=paste0(round(yields, 4) * 100, "%")))
            dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
            print(noquote(tbl))
        }
    })
