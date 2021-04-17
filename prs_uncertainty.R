library(optparse)
library(bigsnpr)
library(tibble)
library(readr)
library(dplyr)

parser <- OptionParser()
parser <- add_option(parser, '--train_bfile', action = 'store', type = 'character')
parser <- add_option(parser, '--train_sumstats', action = 'store', type = 'character')

parser <- add_option(parser, '--val_bfile', action = 'store', type = 'character')
parser <- add_option(parser, '--val_pheno', action = 'store', type = 'character')

parser <- add_option(parser, '--test_bfile', action = 'store', type = 'character')

parser <- add_option(parser, '--out_dir', action = 'store', type = 'character')
parser <- add_option(parser, '--cache_dir', action = 'store', type = 'character')

parser <- add_option(parser, '--chr_i', action = 'store', type = 'integer')
parser <- add_option(parser, '--n_cores', action = 'store', type = 'integer')
parser <- add_option(parser, '--output_train_val_num', action = 'store', type = 'integer', default=0, help="How many individuals information should the program output")
parser <- add_option(parser, '--num_burn_in', action = 'store', type = 'integer', default=100)
parser <- add_option(parser, '--num_iter', action = 'store', type = 'integer', default=500)
parser <- parse_args(parser)

print(parser)

# ---------------------------------------------------------------------------------- 
# train the model
# ---------------------------------------------------------------------------------- 

## load training data
train_snps = read_tsv(paste0(parser$train_bfile, '.bim'), 
        col_names=c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'),
        col_types='icdicc')

train_rds_path = snp_readBed2(paste0(parser$train_bfile, '.bed'),
    backingfile = tempfile(),
    ind.col = which(train_snps$CHR == parser$chr_i))

train_bigSNP <- snp_attach(train_rds_path)

# extract data for analyses
train_G  <- snp_fastImputeSimple(train_bigSNP$genotypes, 
        method = 'mean2', ncores = parser$n_cores)
CHR <- train_bigSNP$map$chromosome
POS <- train_bigSNP$map$physical.pos

# sumstats
train_sumstats <- read_tsv(parser$train_sumstats, col_types='icicciddddddc') %>%
    filter(CHR == parser$chr_i) %>% 
    select(chr=CHR, rsid=SNP, pos=BP, a0=A1, a1=A2, beta=BETA, beta_se=SE, n_eff=NMISS, p=P)

map <- train_bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(train_sumstats, map, strand_flip = F)

# calculate LD
POS2 <- snp_asGeneticPos(CHR, POS, dir = parser$cache_dir)
df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

# check whether LD precomputed
ld_file_path <- file.path(parser$cache_dir, sprintf("chr%d.ld.rds", parser$chr_i))

if (file.exists(ld_file_path)){
    corr0 <- readRDS(ld_file_path)
}else{
    corr0 <- snp_cor(train_G,
                    infos.pos = POS2, 
                    size = 3 / 1000,
                    ncores = parser$n_cores)
    saveRDS(corr0, ld_file_path)
}

ld <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

# run ldsc to get h2 estimation and parameters grid for tuning
ldsc <- snp_ldsc2(corr0, df_beta)
h2_est <- ldsc[["h2"]]

# create parameters tuning grid
(h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
(p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

# train model on all parameters
beta_grid <- snp_ldpred2_grid(ld, df_beta, params, ncores = parser$n_cores, burn_in = parser$num_burn_in, num_iter = parser$num_iter)

# ---------------------------------------------------------------------------------- 
# Tune heritability and pcausal paramters on validation data
# ----------------------------------------------------------------------------------

# load validation data
val_snps = read_tsv(paste0(parser$val_bfile, '.bim'), 
        col_names=c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'),
        col_types='icdicc')

val_rds_path = snp_readBed2(paste0(parser$val_bfile, '.bed'),
    backingfile = tempfile(),
    ind.col = which(val_snps$CHR == parser$chr_i))

val_bigSNP <- snp_attach(val_rds_path)

# extract data for analyses
val_G  <- snp_fastImputeSimple(val_bigSNP$genotypes, 
        method = 'mean2', ncores = parser$n_cores)
val_pheno <- as_tibble(read.table(parser$val_pheno, header=TRUE)) %>%
    mutate(FID_IID = paste0(FID, '_', IID)) %>%
    select(FID_IID, pheno=regressed)

# predict on the validation set 
pred_grid_val <- as_tibble(big_prodMat(val_G, beta_grid))
pred_grid_val$FID_IID <- paste0(val_bigSNP$fam$family.ID, '_', val_bigSNP$fam$sample.ID)

# select best parameter
val_df <- inner_join(val_pheno, pred_grid_val, by='FID_IID') %>%
    filter(pheno != -9)

corr <- rep(NA, nrow(params))
for (param_i in 1 : nrow(params)){
    corr[param_i] <- cor(val_df$pheno, val_df[[paste0('V', param_i)]])
}

best_model = which(corr == max(corr))
best_param = params[best_model,]
best_beta = matrix(beta_grid[, best_model], ncol = 1)

# ---------------------------------------------------------------------------------- 
# Train the best model again to obtain posterior samplings 
# ----------------------------------------------------------------------------------

best_beta_sample = snp_ldpred2_sampling(ld, df_beta, best_param, 
            ncores = parser$n_cores, 
            burn_in = parser$num_burn_in, 
            num_iter = parser$num_iter)

best_beta_sample = snp_ldpred2_grid(ld, df_beta, best_param,
            ncores = parser$n_cores, 
            burn_in = parser$num_burn_in, 
            num_iter = parser$num_iter,
            return_sampling_betas = TRUE)

# model info
model <- list(h2_est = h2_est,
                info_snp = info_snp, 
                params = params,
                beta_grid = beta_grid,
                corr = corr,
                best_param = best_param,
                best_beta = best_beta,
                best_beta_sample = best_beta_sample)

# ---------------------------------------------------------------------------------- 
# Estimate PRS and uncertainty on testing data
# ----------------------------------------------------------------------------------

# load testing data
test_snps = read_tsv(paste0(parser$test_bfile, '.bim'), 
        col_names=c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'),
        col_types='icdicc')

test_rds_path = snp_readBed2(paste0(parser$test_bfile, '.bed'),
    backingfile = tempfile(),
    ind.col = which(test_snps$CHR == parser$chr_i))

test_bigSNP <- snp_attach(test_rds_path)

test_G  <- snp_fastImputeSimple(test_bigSNP$genotypes, 
        method = 'mean2', ncores = parser$n_cores)

pred_test = big_prodMat(test_G, best_beta)
pred_test <- tibble(FID_IID = paste0(test_bigSNP$fam$family.ID, '_', test_bigSNP$fam$sample.ID), PREDICTION=as.vector(pred_test))
pred_test_sample <- big_prodMat(test_G, best_beta_sample)

rls <- list(model = model, pred_test=pred_test, pred_test_sample=pred_test_sample)

# also output these for training and testing population, only the top 10k individuals
if (parser$output_train_val_num > 0){
    train_val_print_index <- seq(parser$output_train_val_num)
    pred_train <- big_prodMat(train_G, best_beta, ind.row=train_val_print_index)
    pred_train <- tibble(FID_IID = paste0(train_bigSNP$fam$family.ID[train_val_print_index], '_', train_bigSNP$fam$sample.ID[train_val_print_index]), 
                         PREDICTION=as.vector(pred_train))
    pred_train_sample <- big_prodMat(train_G, best_beta_sample, ind.row=train_val_print_index)

    pred_val <- big_prodMat(val_G, best_beta, ind.row=train_val_print_index)
    pred_val <- tibble(FID_IID = paste0(val_bigSNP$fam$family.ID[train_val_print_index], '_', val_bigSNP$fam$sample.ID[train_val_print_index]), 
                         PREDICTION=as.vector(pred_val))
    pred_val_sample <- big_prodMat(val_G, best_beta_sample, ind.row=train_val_print_index)

    rls <- append(rls, list(pred_train=pred_train, pred_train_sample=pred_train_sample, 
                            pred_val=pred_val, pred_val_sample=pred_val_sample))
}


saveRDS(rls, file = file.path(parser$out_dir, sprintf('ldpred2.grid.%d.rds', parser$chr_i)))
