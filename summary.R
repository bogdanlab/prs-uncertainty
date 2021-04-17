library(optparse)

parser <- OptionParser()
parser <- add_option(parser, '--out_dir', action='store', type='character', help="output directory of PRS model")
parser <- parse_args(parser)
print(parser)


# go through the files in check corrupted files
missing <- FALSE
for (chr_i in chr_i_list){
    rds_file <- file.path(parser$out_dir, paste0('ldpred2.grid.', chr_i, '.rds'))
    rds <- tryCatch({
        readRDS(rds_file)
    }, error = function(e){
        print(e)
        print(rds_file)
        unlink(rds_file)
        missing <- TRUE
    })
}
if (missing){
    quit()
}


# merge 22 PRS models 
test_posterior <- 0
beta_samplings <- list()

for (chr_i in chr_i_list){
    rds_file <- file.path(parser$pred_dir, paste0('ldpred2.grid.', chr_i, '.rds'))
    rds <- readRDS(rds_file)

    test_posterior <- test_posterior + rds[["pred_test_sample"]]
    beta_samplings[[chr_i]] <- rls$model$best_beta_sample
}


test_prs <- list(posterior_mean = rowMeans(test_posterior), 
                posterior_samples = test_posterior)


saveRDS(test_prs, file.path(parser$out, 'test_prs.rds'))
saveRDS(beta_samplings, file.path(parser$out, 'beta_sampling.rds'))