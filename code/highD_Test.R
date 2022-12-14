###### env setup ######
#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
library(reshape2)


###### generate / load p-values ######
# simulate / load p-values
simu_pval = FALSE

d = 2:20
p = 2
SNRS = c(4, 0.25)
samples = 200
n = sqrt(c(50, 100, 500))

if (simu_pval) {
  numCores = parallel::detectCores()
  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  data_highD = c()

  cl <- makeCluster(numCores)
  registerDoSNOW(cl)

  for (i in 1:length(d)) {
    set.seed(7202)
    if (i == 1) {cat(paste('Dimension d =', d[i], ': \n'))}
    else {cat(paste('\nDimension d =', d[i], ': \n'))}

    means = generateMeans(d[i], p, snr = SNRS)
    simulated = simuSamples(means, n, samples, prl = T)
    gc()
    cat('\nTesting: \n')

    eigvORC = JDTE(means[,1,,], iter = 1000)

    data_i = foreach(est_list = simulated, m =icount(), .inorder = F,
                     .combine = bind_rows, .options.snow = opts,
                     .packages = c('eigTest', 'reshape2', 'tidyverse')) %dopar% {
                       if (m %% max(100, totL %/% 50) == 1) {gc()}
                       #                       print(m)
                       #                       profvis({
                       mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                       SNR = est_list$SNR; CovRate = est_list$CovRate

                       snr = as.numeric(substring(SNR, 7))

                       data_temp = commutatorTest(mu.bar, cn = CovRate, param.out = T)[[1]] %>%
                         list2DF() %>% select(df, pvalue) %>% mutate(covType = 'Oracle')
                       data_temp = commutatorTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                                  param.out = T)[[1]] %>%
                         list2DF() %>% select(df, pvalue) %>% mutate(covType = 'Plugin') %>%
                         bind_rows(data_temp) %>%
                         mutate(SNR = paste0('SNR = ', ifelse(snr == 0, '$\\infty$', round(1/snr))),
                                Dimension = d[i], SampleSize = round(CovRate^2), testType = 'COM')

                       eigvPLG = JDTE(mu.bar)

                       data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                           CV = eigvPLG, poly.sp = F) %>%
                         list2DF() %>% select(df, pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'eigv-PLG')
                       data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                           poly.sp = F, CV = eigvPLG, param.out = T) %>%
                         list2DF() %>% select(df, pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'eigv-PLG') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, param.out = T) %>% list2DF() %>%
                         select(df, pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'poly-PLG') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                           param.out = T) %>%
                         list2DF() %>% select(df, pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'poly-PLG') %>%
                         bind_rows(data_LLR)

                       if (snr == 0) {
                         data_LLR = projTest(mu.bar, refMat = means[,1,,],
                                             cn = CovRate, param.out = T) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Oracle', spaceType = 'poly-ORC') %>%
                           bind_rows(data_LLR)
                         data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                             CV = eigvORC, poly.sp = F) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Oracle', spaceType = 'eigv-ORC') %>%
                           bind_rows(data_LLR)
                         data_LLR = projTest(mu.bar, cn = CovRate, refMat = means[,1,,],
                                             cov.arr = cov.bar, param.out = T) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Plugin', spaceType = 'poly-ORC') %>%
                           bind_rows(data_LLR)
                         data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                             CV = eigvORC, poly.sp = F, param.out = T) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Plugin', spaceType = 'eigv-ORC') %>%
                           bind_rows(data_LLR)
                       }
                       #                       })
                       return(data_LLR %>%
                                mutate(SNR = paste0('SNR = ', ifelse(snr == 0, '$\\infty$', round(1/snr))),
                                       Dimension = d[i],
                                       SampleSize = paste0('Sample size $n = ', round(CovRate^2), '$'),
                                       testType = 'LLR') %>%
                                bind_rows(data_temp)
                       )
                     }
    rm(simulated); gc()
    data_highD = rbind(data_highD, data_i)
  }
  stopCluster(cl)
  save(data_highD, file = 'output/highD_test.RData')
} else {
  load('output/highD_Test.RData')
}


###### Test sizes ######
data_highD_summary = data_highD %>%
  group_by(covType, spaceType, Dimension, SNR, SampleSize, testType) %>%
  summarise(RejRate = mean(pvalue <= 0.05), df = mean(df)) %>%
  mutate(covType = paste(covType, 'Cov'),
         spaceType = gsub('poly-ORC', 'Oracle (A.5)', spaceType),
         spaceType = gsub('poly-PLG', 'Plugin (A.5)', spaceType),
         spaceType = gsub('eigv-ORC', 'Oracle (A.6)', spaceType),
         spaceType = gsub('eigv-PLG', 'Plugin (A.6)', spaceType)) %>%
  ungroup %>% print(n = nrow(.))

