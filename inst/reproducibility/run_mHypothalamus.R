source("R/funcs_.R")
dir.input <- file.path('data/mHypothalamus/')

filename = paste0(dir.input, 'MERFISH_Animal1_cnts.xlsx')
infoname = paste0(dir.input, 'MERFISH_Animal1_info.xlsx')
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')
results <- list()
aris <- list()

#for (neighborhood_threshold in c(1,2,3)) {
#  for (embedding_dim in c(5,10,15)) {
for (lambda in seq(0, 0.2, by = 0.05)) {
#for (sheet in sheets) {
  sheet <- '-0.24'
  print(paste0("Slice ", sheet))
  cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(cnts) <- cnts[,"...1"]
  cnts <- cnts[ -c(1) ]
  
  xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(xys) <- xys[,"...1"]
  gtlabels <- xys$z
  #gtlabels <- xys$Cell_class
  xys <- xys[-c(1)]
  xys <- xys[-c(-2:-1)]
  xys <- xys[,c(2,1)]
  R <- length(unique(gtlabels))
  
  selected_genes <- "HVG"
  neighborhood_threshold = 3
  embedding_dim = 10
  #lambda = 0.2
  lr = 16
  max_iterations = 20
  loss_tol = 1e-5
  eps = 1e-20
  refinement = FALSE
  refine_k <- 30
  
  results[[as.character(lambda)]] <- DOST(cnts, xys, R, selected_genes = selected_genes,
                  neighborhood_threshold = neighborhood_threshold,
                  embedding_dim = embedding_dim, lambda = lambda,
                  lr = lr, max_iterations = max_iterations,
                  loss_tol = loss_tol, eps = eps,
                  refinement = refinement)
  #aris[[sheet]] <- mclust::adjustedRandIndex(gtlabels, results[[sheet]]$labels)
  #write(unlist(aris), file = paste0("results/mHypothalamus_DOST_", lambda,"_ari.txt"), ncolumns = 1)
  #line <- paste(sheet, neighborhood_threshold, embedding_dim, lambda, aris[[sheet]], sep = ",")
  #write(line, file = "results/mHypo_DOST.txt", append = TRUE)
#}
#line <- paste("mean", neighborhood_threshold, embedding_dim, lambda, mean(unlist(aris)), sep = ",")
#write(line, file = "results/mHypo_DOST.txt", append = TRUE)
}
#  }
#  }

#plotMHypo(xys, results[[sheets[[1]]]]$labels)
#plot_umap(results[[sheets[[1]]]]$Z, gtlabels)
save(results, aris, file = "RData/results_mHypothalamus_DOST_ablation_spatial_domain_-0.24.RData")
