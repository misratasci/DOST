
source("R/funcs_.R")
dir.input <- file.path('data/STARmap_mouse_PFC/')

filename = paste0(dir.input, 'starmap_mpfc_starmap_cnts.xlsx')
infoname = paste0(dir.input, 'starmap_mpfc_starmap_info.xlsx')

sheets <- readxl::excel_sheets(filename)
results <- list()
aris <- list()

#for (neighborhood_threshold in c(1,2,3)) {
#  for (embedding_dim in c(5,10,15)) {
#for (lambda in seq(0.1, 0.9, by = 0.05)) {
for (sheet in sheets) {
  cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(cnts) <- cnts[,"...1"]
  cnts <- cnts[ -c(1) ]
  
  xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(xys) <- xys[,"...1"]
  gtlabels <- xys$z
  xys <- xys[-c(1)]
  xys <- xys[-c(-2:-1)]
  xys <- xys[,c(2,1)]
  R <- length(unique(gtlabels))
  
  selected_genes <- "HVG"
  lr = 16
  max_iterations = 20
  loss_tol = 1e-5
  eps = 1e-20
  refinement = TRUE
  refine_k <- 30
  neighborhood_threshold = 1
  embedding_dim = 10
  lambda = 0.15
  
  results[[sheet]] <- DOST(cnts, xys, R, selected_genes = selected_genes,
                           neighborhood_threshold = neighborhood_threshold,
                           embedding_dim = embedding_dim, lambda = lambda,
                           lr = lr, max_iterations = max_iterations,
                           loss_tol = loss_tol, eps = eps,
                           refinement = refinement)
  aris[[sheet]] <- mclust::adjustedRandIndex(gtlabels, results[[sheet]]$labels)
  #line <- paste(sheet, neighborhood_threshold, embedding_dim, lambda, aris[[sheet]], sep = ",")
  #write(line, file = "results/mPFC_DOST.txt", append = TRUE)
}
#line <- paste("mean", neighborhood_threshold, embedding_dim, lambda, mean(unlist(aris)), sep = ",")
#write(line, file = "results/mPFC_DOST.txt", append = TRUE)
#}
#  }
#  }

save(results, aris, file = "RData/results_mPFC_DOST.RData")
