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
  means_c = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)

  numCores = parallel::detectCores()
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  simulated_c = simuSamples(means_c, n, samples, prl = T)

  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)

  data_c = foreach(est_list = simulated_c, .inorder = F, .combine = bind_rows,
                   .options.snow = opts, .packages = c('eigTest', 'tidyverse')) %dopar% {
                     mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                     SNR = as.numeric(str_split(est_list$SNR, '=')[[1]][2])
                     CovRate = est_list$CovRate

                     data_temp = commutatorTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                                param.out = T)[[1]] %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(SNR = paste0('SNR = ', ifelse(SNR==0, '$\\infty$', round(1/SNR))),
                              SampleSize = paste0('Sample size $n = ', round(CovRate^2), '$'))

                     return(data_temp)
                   }
  stopCluster(cl)
  orders = order(unique(data_c$SNR))
  orders = c(orders[1], rev(orders[2:length(orders)]))
  data_c$SNR = factor(data_c$SNR, levels = unique(data_c$SNR)[orders])
  data_c$SampleSize = factor(data_c$SampleSize, levels = paste0('Sample size $n = ', round(n^2), '$'))
  save(data_c, file = 'output/commutatorTest.RData')
} else {
  load('output/commutatorTest.RData')
}


###### plot histograms ######
binwidth = 0.05
breaks = seq(0, 1, 0.05)

tikz(file = "output/Plots/PvalueCommutator.tikz", standAlone=F,width = 6, height = 3)
ggplot(data_c) +
  geom_histogram(aes(x = pvalue, y = ..density..*binwidth, fill = SNR),
                 breaks = breaks,
                 position = position_dodge()) + theme_bw()+
  ggtitle('P-value histogram from Commutator-based test') +
  facet_wrap(~SampleSize, ncol = length(n), dir = 'v') +scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')+
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(size = 1), aspect.ratio = 1,
        panel.spacing.x = unit(4, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        strip.background =element_rect(fill="white", size = 0)) +
  guides(fill=guide_legend(ncol=3)) + xlab('P values') + ylab('Proportion') +
  labs(fill = 'SNR:') #scale_fill_discrete(labels = c('Chi', 'Gam'))+
dev.off()


###### Type I/II errors ######
data_c %>% group_by(SNR, SampleSize) %>%
  summarise(RejRate = mean(pvalue <= 0.05)) %>% print(n = nrow(.))

