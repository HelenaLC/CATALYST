# ==============================================================================
# Show method for dbFrame
# ------------------------------------------------------------------------------

setMethod(f="show",
    signature="dbFrame",
    definition=function(object) {
        dims <- dim(exprs(object))
        cat("dbFrame objectect with\n", 
            dims[1], " events, ", 
            dims[2], " observables and ",
            nrow(bc_key(object)), " barcodes:\n", sep="")
        
        cat("\nCurrent assignments:\n")
        cat("     ", sum(bc_ids(object) == 0), "event(s) unassigned")
        
        tmp <- c(sort(table(bc_ids(object)[bc_ids(object) != 0]), 
            decreasing=TRUE))
        tbl <- t(data.frame("ID"=names(tmp), "Count"=tmp))
        dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
        print(noquote(tbl))
        
        if (length(sep_cutoffs(object)) != 0) {
            inds <- match(names(tmp), rownames(bc_key(object)))
            cat("\nSeparation cutoffs:")
            tbl <- t(data.frame(
                "ID"=names(tmp), 
                "Cutoff"=sep_cutoffs(object)[inds]))
            dimnames(tbl)[2] <- list(rep("", dim(tbl)[2]))
            print(noquote(tbl))
            
            # better way to do this?
            seps <- seq(0, 1, .01)
            ind <- NULL
            for (i in sep_cutoffs(object)[inds]) {
                for (j in seps) {
                    if (isTRUE(all.equal(i, j))) 
                        ind <- append(ind, which(seps == j))
                }
            }
            yields <- yields(object)[cbind(inds, ind)]
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
