#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
library(tikzDevice)


###### generate / load pvalues ######
# simulate pvalues from scratch
simu_pval = T#FALSE

d = 5
p = 2
k = d
samples = 500
SNRS = c(1000, 10)
n = sqrt(c(50, 250, 1000))

if (simu_pval) {
  set.seed(4000)
  means_l = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)

  numCores = parallel::detectCores()
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  simulated_l = simuSamples(means_l, n, samples, prl = T)

  eigvORC = JDTE(means_l[,1,,], iter = 1000)

  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)

  data_l = foreach(est_list = simulated_l, .inorder = F, .combine = bind_rows,
                   .options.snow = opts, .packages = c('eigTest', 'tidyverse')) %dopar% {
                     mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                     SNR = as.numeric(str_split(est_list$SNR, '=')[[1]][2])
                     CovRate = est_list$CovRate

                     eigvPLG = JDTE(mu.bar)

                     data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                         CV = eigvPLG, poly.sp = F) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(covType = 'Oracle', spaceType = 'Plugin eqref{eqn:eigvP}')
                     data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                         poly.sp = F, CV = eigvPLG, param.out = T) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(covType = 'Plugin', spaceType = 'Plugin eqref{eqn:eigvP}') %>%
                       bind_rows(data_LLR)
                     data_LLR = projTest(mu.bar, cn = CovRate, param.out = T) %>% list2DF() %>%
                       select(pvalue) %>%
                       mutate(covType = 'Oracle', spaceType = 'Plugin eqref{eqn:polyP}') %>%
                       bind_rows(data_LLR)
                     data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                         param.out = T) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(covType = 'Plugin', spaceType = 'Plugin eqref{eqn:polyP}') %>%
                       bind_rows(data_LLR)

                     if (SNR == 0) {
                       data_LLR = projTest(mu.bar, refMat = means_l[,1,,],
                                           cn = CovRate, param.out = T) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'Oracle eqref{eqn:polyP}') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                           CV = eigvORC, poly.sp = F) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'Oracle eqref{eqn:eigvP}') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, refMat = means_l[,1,,],
                                           cov.arr = cov.bar, param.out = T) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'Oracle eqref{eqn:polyP}') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                           CV = eigvORC, poly.sp = F, param.out = T) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'Oracle eqref{eqn:eigvP}') %>%
                         bind_rows(data_LLR)
                     }
                     #                       })
                     return(data_LLR %>%
                              mutate(SNR = paste0('SNR = ', ifelse(SNR==0, '$\\infty$', round(1/SNR))),
                                     SampleSize = paste0('Sample size $n = ', round(CovRate^2), '$'))
                     )
                   }
  stopCluster(cl)
  orders = order(unique(data_l$SNR))
  orders = c(orders[1], rev(orders[2:length(orders)]))
  data_l$SNR = factor(data_l$SNR, levels = unique(data_l$SNR)[orders])
  data_l$SampleSize = factor(data_l$SampleSize, levels = paste0('Sample size $n = ', round(n^2), '$'))
  data_l$spaceType = factor(data_l$spaceType, levels = c('Oracle eqref{eqn:polyP}',
                                                         'Plugin eqref{eqn:polyP}',
                                                         'Oracle eqref{eqn:eigvP}',
                                                         'Plugin eqref{eqn:eigvP}'))
  save(data_l, file = 'output/LLRTest.RData')
} else {
  load('output/LLRTest.RData')
}


###### plot histograms ######
binwidth = 0.05
breaks = seq(0, 1, 0.05)

tikz(file = "output/Plots/PvalueLLR.tikz", standAlone=F, width = 6, height = 7.5)
ggplot(data_l %>% filter(covType == 'Plugin')) +
  geom_histogram(aes(x = pvalue, y = ..density..*binwidth, fill = SNR),
                 breaks = breaks,
                 position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from LLR test') +
  facet_grid(vars(spaceType), vars(SampleSize)) +
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')+
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(size = 1), aspect.ratio = 1,
        panel.spacing.x = unit(4, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        strip.background =element_rect(fill="white", size = 0)) +
  guides(fill=guide_legend(ncol=3)) + xlab('P values') + ylab('Proportion') +
  labs(fill = 'SNR:')
dev.off()


###### Type I/II errors ######
data_l %>% filter(covType == 'Plugin') %>%
  group_by(spaceType, SNR, SampleSize) %>%
  summarise(RejRate = mean(pvalue <= 0.05)) %>% print(n = nrow(.))

