test_that("computeSpillmat() works.", {
    
    get_mass_from_channel <- function(channel){
        return( as.numeric(gsub("[[:punct:][:alpha:]]", "", channel)))
    }
    
    get_metal_from_channel <- function(channel){
        return(gsub("[[:digit:]].*", "", channel))
    }
    # Generate some dummy testdata with ss_exp as a template
    data(ss_exp)
    set.seed(12345)
    ncells = 11
    
    
    # generate a known spillover matrix
    ncol = flowCore::ncol(ss_exp)
    sm = diag(1, ncol, ncol)
    sm_names <- flowCore::colnames(ss_exp)
    colnames(sm) <- rownames(sm) <- sm_names
    # add metal, +-1, +16 specific spillover
    sm_mass = get_mass_from_channel(sm_names)
    sm_met = get_metal_from_channel(sm_names)
    
    for (i in seq_along(sm_mass)){
        curname = sm_names[i]
        curmass = sm_mass[i]
        curmet = sm_met[i]
        pot_spill = sm_names[((sm_met == curmet) | (sm_mass == curmass +1)|
                         (sm_mass == curmass +16) | (sm_mass == curmass-1)
                     ) & (sm_names != curname)]
        sm[curname, pot_spill] <- sm[curname, pot_spill] + runif(length(pot_spill))/20
    }
    
    # generate 100 cells that are positive for each channel
    mat = diag(x=1, nrow=ncol, ncol=ncol)
    mat = do.call(rbind, replicate(ncells, mat, simplify=FALSE)) 
    
    # add some noise to the cells
    mat = mat * (rnorm(ncol*ncells, 1, 0.01)) *100
    
    # add spillover
    mat_sm = mat %*% sm
    
    # create a flowframe from it
    ff_test = flowCore::flowFrame(mat_sm)
    
    # start debarcoding
    bc_ms <- sm_mass
    # debarcode
    re <- assignPrelim(x=ff_test, y=bc_ms, verbose=FALSE)
    
    # check if debarcoding worked flawelessly
    re@bc_ids = rep(as.character(sm_mass), ncells)
    # compute spillover matrix
    spillMat <- computeSpillmat(x=re)

    expect_equal(spillMat, sm)
    
    ## Repeat the test with a constant amount of background
    bg=20
    ff_test = flowCore::flowFrame(mat_sm+bg)
    
    # start debarcoding
    bc_ms <- sm_mass
    # debarcode
    re <- assignPrelim(x=ff_test, y=bc_ms, verbose=FALSE)
    
    # check if debarcoding worked flawelessly
    re@bc_ids = rep(as.character(sm_mass), ncells)
    # compute spillover matrix
    spillMat <- computeSpillmat(x=re)
    
    expect_equal(spillMat, sm)
})
    