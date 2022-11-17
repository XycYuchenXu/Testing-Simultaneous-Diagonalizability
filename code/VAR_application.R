###### env setup ######
library(tidyverse)
library(plotrix)
library(reshape2)
library(tikzDevice)
library(WDI)
library(OECD)
library(x12)
library(pdfetch)
library(jsonlite)
library(latex2exp)

countryNames = c('CHN', 'JPN', 'KOR', 'FRA', 'DEU', 'GBR', 'CAN', 'USA')
indNames = c('GDP', 'M2', 'REER')

###### time series download ######
# GDP from World Bank
GDP = WDI(indicator = 'NYGDPMKTPSACD',
          country = c('CN', 'JP', 'KR', 'FR', 'DE', 'GB', 'CA', 'US'),
          start = '1992Q1', end = '2020Q1') %>%
  select(iso2c, year, NYGDPMKTPSACD) %>%
  rename(country = iso2c, Time = year, gdp = NYGDPMKTPSACD) %>%
  mutate(Time = gsub('Q1', '/03/01', Time),
         Time = gsub('Q2', '/06/01', Time),
         Time = gsub('Q3', '/09/01', Time),
         Time = gsub('Q4', '/12/01', Time),
         Time = as.Date(Time)) %>% arrange(Time) %>%
  pivot_wider(id_cols = Time, names_from = country, values_from = gdp)

# money supply M2 from central banks, FRB & IMF
# USA
M2_USA = read_csv('https://www.federalreserve.gov/datadownload/Output.aspx?rel=H6&series=798e2796917702a5f8423426ba7e6b42&lastobs=&from=01/01/1992&to=03/31/2020&filetype=csv&label=include&layout=seriescolumn',
                  skip = 5) %>% rename(Time = `Time Period`, USA = M2.M) %>%
  mutate(Time = as.Date(paste0(Time, '-01'))) %>%
  filter(Time < '2020-04-01', Time >= '1992-01-01') %>%
  group_by(cut(Time, 'quarter')) %>%
  summarise(Time = last(Time), USA = last(USA)) %>%
  select(Time, USA)

# CAN
M2_CAN = read_csv('https://www.bankofcanada.ca/valet/observations/group/e1_monthly/csv?start_date=1946-01-01',
                  skip = 51) %>% rename(Time = date, CAN = V41552796) %>%
  filter(Time < '2020-04-01', Time >= '1992-01-01') %>%
  group_by(cut(Time, 'quarter')) %>%
  summarise(Time = last(Time), CAN = last(CAN) / 1000) %>%
  select(Time, CAN)

# FRA
M2_FRA = read_csv('https://webstat.banque-france.fr/ws_wsen/ng/exportApi/file?node=5384976&DATASET=BSI1&DATASET=ECOFI&periodSortOrder=DESC&SERIES_KEY=BSI1.M.FR.N.V.M20.A.1.U2.2300.Z01.E&exportType=csv') %>% drop_na() %>% .[-(1:5),]
colnames(M2_FRA) = c('Time', 'FRA')
M2_FRA = M2_FRA %>% mutate(FRA = as.numeric(FRA) / 1000) %>% drop_na() %>%
  mutate(Time = rev(seq(as.Date('1980-01-01'), by = 'month', length = nrow(.)))) %>%
  arrange(Time)
m2 = new('x12Single', ts = ts(M2_FRA$FRA, frequency = 12, start = c(1980, 1)))
M2_FRA = M2_FRA %>% mutate(FRA = x12(m2)@x12Output@d11 %>% as.numeric()) %>%
  filter(Time >= '1992-01-01', Time < '2020-04-01') %>%
  group_by(cut(Time, 'quarter')) %>% summarise(FRA = last(FRA), Time = last(Time)) %>%
  select(Time, FRA)

# DEU
M2_DEU = pdfetch_BUNDESBANK('BBK01.TXI302') %>%
  data.frame(Time = time(.)) %>% as_tibble() %>%
  rename(DEU = 1) %>%
  mutate(DEU = DEU / 1000, Time = Time + 1 - months(1))
m2 = new('x12Single', ts = ts(M2_DEU$DEU, frequency = 12, start = c(1980, 1)))
M2_DEU = M2_DEU %>% mutate(DEU = x12(m2)@x12Output@d11 %>% as.numeric()) %>%
  group_by(cut(Time, 'quarter')) %>% summarise(DEU = last(DEU), Time = last(Time)) %>%
  select(Time, DEU) %>% filter(Time >= '1992-01-01', Time < '2020-04-01')

# GBR
M2_GBR = pdfetch_BOE('LPQVWYW', from = '1992-01-01', to = '2020-04-01') %>%
  data.frame(Time = time(.)) %>% as_tibble() %>%
  rename(GBR = 1) %>%
  mutate(GBR = GBR / 1000, Time = Time + 1 - months(1))

# JPN
M2_JPN = read_csv('https://www.stat-search.boj.or.jp/ssi/mtshtml/csv/md02_m_1_en.csv',
                  skip = 3) %>% select(1, starts_with('M2')) %>% .[-(1:5),]
colnames(M2_JPN) = c('Time', 'Rate', 'JPN')
M2_JPN = M2_JPN %>% mutate(Time = as.Date(paste0(Time, '/01')),
                           Rate = as.numeric(Rate),
                           JPN = as.numeric(JPN)/10) %>%
  mutate(Month = months(Time)) %>% group_by(Month) %>%
  mutate(Rate = cumprod(Rate / 100 + 1),
         first_avai_M2 = first(na.omit(JPN)),
         first_avai_rate = Rate[which(!is.na(JPN))[1]],
         est_M2 = first_avai_M2 * Rate / first_avai_rate,
         JPN = coalesce(JPN, est_M2)) %>% ungroup()
m2 = new('x12Single', ts = ts(M2_JPN$JPN, frequency = 12, start = c(1980, 1)))
M2_JPN = M2_JPN %>%
  mutate(JPN = x12(m2)@x12Output@d11 %>% as.numeric()) %>%
  group_by(cut(Time, 'quarter')) %>% summarise(JPN = mean(JPN), Time = last(Time)) %>%
  select(Time, JPN) %>%
  filter(Time < '2020-04-01', Time >= '1992-01-01')

# KOR
M2_KOR = read_csv('data/M2_KOR_monthly_SA.csv') %>%
  mutate(Time = as.Date(paste0(Time, '/01'))) %>%
  group_by(cut(Time, 'quarter')) %>% summarise(KOR = last(KOR), Time = last(Time)) %>%
  select(Time, KOR)

# CHN
M2_CHN = fromJSON('http://dataservices.imf.org/REST/SDMX_JSON.svc/CompactData/IFS/Q.CN.35L___XDC') %>%
  as.data.frame() %>% as_tibble() %>% select(last_col(1), last_col())
colnames(M2_CHN) = c('Time', 'CHN')
M2_CHN = M2_CHN %>%
  mutate(Time = gsub('Q1', '03-01', Time),
         Time = gsub('Q2', '06-01', Time),
         Time = gsub('Q3', '09-01', Time),
         Time = gsub('Q4', '12-01', Time),
         Time = as.Date(Time), CHN = as.numeric(CHN) / 1000)
m2 = new('x12Single', ts = ts(M2_CHN$CHN, frequency = 4, start = c(1978, 1)))
M2_CHN = M2_CHN %>% mutate(CHN = x12(m2)@x12Output@d11 %>% as.numeric()) %>%
  filter(Time < '2020-04-01', Time >= '1992-01-01')

# remove x12 files from SA
unlink(list.files(pattern = 'Rout'), recursive = T)

# join M2 for all countries
M2_LCU = list(M2_USA, M2_CAN, M2_FRA, M2_DEU, M2_GBR, M2_JPN, M2_KOR, M2_CHN) %>%
  reduce(left_join, by = 'Time') %>%
  pivot_longer(cols = - Time, names_to = 'country', values_to = 'M2')

# exchange rate from World Bank
ER = WDI(indicator = 'DPANUSSPB',
         country = c('US', 'CA', 'FR', 'DE', 'GB', 'JP', 'KR', 'CN'),
         start = '1992M01', end = '2020M03') %>% select(iso2c, year, DPANUSSPB) %>%
  rename(country = iso2c, Time = year, er = DPANUSSPB) %>%
  mutate(Time = gsub('M', '/', Time),
         Time = paste0(Time, '/01'),
         Time = as.Date(Time)) %>% arrange(Time) %>%
  group_by(country, cut(Time, 'quarter')) %>%
  summarise(er = last(er), Time = last(Time)) %>% select(country, er, Time)

# convert money supply to USD
M2 = M2_LCU %>% left_join(ER) %>%
  mutate(M2 = M2 / er) %>% select(-er) %>%
  pivot_wider(id_cols = Time, names_from = country, values_from = M2)

# REER from OECD
REER_nsa = get_dataset('MEI', filter = 'CAN+FRA+DEU+JPN+KOR+GBR+USA+CHN.CCRETT01.IXOB.Q',
                   start_time = 1992, end_time = 2020, pre_formatted = T) %>%
  select(LOCATION, ObsValue, Time) %>%
  mutate(ObsValue = as.numeric(ObsValue),
         Time = gsub('Q1', '03-01', Time),
         Time = gsub('Q2', '06-01', Time),
         Time = gsub('Q3', '09-01', Time),
         Time = gsub('Q4', '12-01', Time),
         Time = as.Date(Time)) %>% arrange(Time) %>% filter(Time < '2020-04-01') %>%
  pivot_wider(id_cols = Time, names_from = LOCATION, values_from = ObsValue)

REER = REER_nsa
for (cc in countryNames) {
  reer = new('x12Single', ts = ts(REER[[cc]], frequency = 4, start = c(1992, 1)))
  REER[[cc]] = x12(reer)@x12Output@d11 %>% as.numeric()
}
unlink(list.files(pattern = 'Rout'), recursive = T)


###### time series plot ######
# merge macro-economic time series
c2p = list(
  GDP %>% mutate_at(vars(-Time), ~./10^6) %>% mutate(key = 'GDP'),
  M2 %>% mutate_at(vars(-Time), ~./10^4) %>% mutate(key = 'M2'),
  mutate(REER, key = 'REER')
) %>% reduce(bind_rows) %>% rename(Quarters = Time) %>%
  pivot_longer(cols = -c(Quarters, key), names_to = 'country', values_to = 'value')

labNames = c('GDP (USD, $\\times 10^{12}$)', 'M2 (USD, $\\times 10^{13}$)', 'REER ($2015=100$)')
names(labNames) = indNames

c2p$country = factor(c2p$country, levels = rev(countryNames))

# generate tikz file
tikz('output/Plots/tikz/ts.tikz', standAlone = F, width = 6, height = 3)
ggplot(data = c2p) +
  geom_path(aes(x = Quarters, y = value, color = country, linetype = country), linewidth = 1.5) +
  facet_wrap(~key, scales = 'free_y', labeller = labeller(key = labNames)) +
  labs(color = 'Country:', linetype = 'Country:') + theme_bw() + ylab('') +
  ggtitle('Quarterly macroeconomic indices') +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank(), legend.position = 'bottom',
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.key.width = unit(1.5, 'cm'), panel.border = element_rect(linewidth = 1),
        legend.spacing.y = unit(0.1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        strip.background =element_rect(fill="white", linewidth = 0)) +
  guides(col=guide_legend(ncol=4))
dev.off()

# generate png file
p1 = ggplot(data = c2p %>% mutate(key = factor(key, levels = indNames, labels = labNames))) +
  geom_path(aes(x = Quarters, y = value, color = country, linetype = country), linewidth = 1.5) +
  facet_wrap(~key, scales = 'free_y', labeller = as_labeller(TeX, default = label_parsed)) +
  labs(color = 'Country:', linetype = 'Country:') + theme_bw() + ylab('') +
  ggtitle('Quarterly macroeconomic indices') +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank(), legend.position = 'bottom',
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.key.width = unit(1.5, 'cm'), panel.border = element_rect(linewidth = 1),
        legend.spacing.y = unit(0.1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        strip.background =element_rect(fill="white", linewidth = 0)) +
  guides(col=guide_legend(ncol=4))
p1
ggsave(filename = 'output/Plots/png/ts.png', p1, width = 6, height = 3, units = 'in')


###### prepare VAR data ######
library(MTS)
USA = cbind(
  GDP_USD_mn = GDP$USA,
  M2_USD_bn = M2$USA,
  REER_ind_2015 = REER$USA
)
CAN = cbind(
  GDP_USD_mn = GDP$CAN,
  M2_USD_bn = M2$CAN,
  REER_ind_2015 = REER$CAN
)
FRA = cbind(
  GDP_USD_mn = GDP$FRA,
  M2_USD_bn = M2$FRA,
  REER_ind_2015 = REER$FRA
)
DEU = cbind(
  GDP_USD_mn = GDP$DEU,
  M2_USD_bn = M2$DEU,
  REER_ind_2015 = REER$DEU
)
GBR = cbind(
  GDP_USD_mn = GDP$GBR,
  M2_USD_bn = M2$GBR,
  REER_ind_2015 = REER$GBR
)
JPN = cbind(
  GDP_USD_mn = GDP$JPN,
  M2_USD_bn = M2$JPN,
  REER_ind_2015 = REER$JPN
)
KOR = cbind(
  GDP_USD_mn = GDP$KOR,
  M2_USD_bn = M2$KOR,
  REER_ind_2015 = REER$KOR
)
CHN = cbind(
  GDP_USD_mn = GDP$CHN,
  M2_USD_bn = M2$CHN,
  REER_ind_2015 = REER$CHN
)

countryMacro = list(
  CHN = CHN, JPN = JPN, KOR = KOR, FRA = FRA, DEU = DEU, GBR = GBR, CAN = CAN, USA = USA
)
save(countryMacro, file = 'output/countryMacro.RData')

# VAR model selection & fit check, on scaled log-difference series
d = ncol(USA)
p = length(countryNames)
n = nrow(USA)
i = 1

bestOrders = matrix(0, p, 3); rownames(bestOrders) = countryNames; colnames(bestOrders) = c('AIC', 'BIC', 'HQ')
coeffAbsEig = matrix(0, p, d); rownames(coeffAbsEig) = countryNames

for (cc in countryNames) {
  indexes = countryMacro[[cc]] %>%
    log() %>% diff() %>% scale() %>% ts(frequency = 4, start = c(1992,1))
  cn = nrow(indexes)
  bestOrder = VARorder(indexes, output =F)
  bestOrders[i,] = c(bestOrder$aicor, bestOrder$bicor, bestOrder$hqor)

  c.var = VAR(indexes, p = 1, output = F, include.mean = T)
  coeff.m = c.var$Phi
  coeffAbsEig[i,] = abs(eigen(coeff.m)$values)
  i = i + 1
}
print(bestOrders); print(coeffAbsEig)


###### test preprocessing ######
#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
m = nrow(countryMacro[[1]]) - 1
d = ncol(countryMacro[[1]])
p = length(countryMacro)

# transpose for test setup
# countryCoeff & countryCovar are also available from package eigTest
countryCoeff_t = countryCoeff
for (i in 1:p) {
  countryCoeff_t[i,,] = t(countryCoeff[i,,])
}
countryCovar_t = countryCovar[,
                              as.vector(matrix(1:d^2, ncol = d, byrow = T)),
                              as.vector(matrix(1:d^2, ncol = d, byrow = T))]


###### global test & estimation ######
# multi-sample test, full-rank, Corollary 4.1
eigTest(countryCoeff_t, cn = sqrt(m), cov.arr = countryCovar_t, testType = 'chi', param.out = T)
eigTest(countryCoeff_t, cn = sqrt(m), cov.arr = countryCovar_t, testType = 'gam', param.out = T)

# partial test, k = 1, Proposition 5.2 & Corollary 5.1
partialTest(countryCoeff_t, cn = sqrt(m), k = 1, cov.arr = countryCovar_t, testType = 'chi', param.out = T)
partialTest(countryCoeff_t, cn = sqrt(m), k = 1, cov.arr = countryCovar_t, testType = 'gam', param.out = T)

# partial test, k = 2, Corollary 5.1
partialTest(countryCoeff_t, cn = sqrt(m), k = 2, cov.arr = countryCovar_t, testType = 'gam', param.out = T)

# estimate common eigenvectors
V = JDTE(countryCoeff_t)
V

# estimate partially common eigenvectors
Q = expmPartSchur(countryCoeff_t, k = 2, warmup = T)
B = array(0, dim = c(p, d, d))
for (i in 1:p) {
  B[i,,] = tcrossprod(crossprod(Q, countryCoeff_t[i,,]), t(Q))
}
V = JDTE(B[,1:2, 1:2])
Q[,1:2] %*% V


###### continent-grouped test ######
# Asia, Corollary 4.1 (Proposition 4.2)
eigTest(countryCoeff_t[1:3,,], cn = sqrt(m), cov.arr = countryCovar_t[1:3,,], testType = 'chi', param.out = T)
eigTest(countryCoeff_t[1:3,,], cn = sqrt(m), cov.arr = countryCovar_t[1:3,,], testType = 'gam', param.out = T)

# Europe, Corollary 4.1 (Proposition 4.2)
eigTest(countryCoeff_t[4:6,,], cn = sqrt(m), cov.arr = countryCovar_t[4:6,,], testType = 'chi', param.out = T)
eigTest(countryCoeff_t[4:6,,], cn = sqrt(m), cov.arr = countryCovar_t[4:6,,], testType = 'gam', param.out = T)

# North America, Corollary 4.1 (Proposition 4.2)
eigTest(countryCoeff_t[7:8,,], cn = sqrt(m), cov.arr = countryCovar_t[7:8,,], testType = 'chi', param.out = T)
eigTest(countryCoeff_t[7:8,,], cn = sqrt(m), cov.arr = countryCovar_t[7:8,,], testType = 'gam', param.out = T)


###### pairwise commutator test ######
# estimate p-value matrix
countries = names(countryMacro)
comm.pair.test = 0.5*diag(length(countries))
for (i in 1:length(countries)) {
  for (j in i:length(countries)) {
    if (i == j) {next()}
    comm.pair.test[i,j] = commutatorTest(countryCoeff[c(i,j),,],
                                         cn = sqrt(m),
                                         cov.arr = countryCovar[c(i,j),,])
  }
}
colnames(comm.pair.test) = countries
rownames(comm.pair.test) = countries
comm.pair.test

# structure for heatmap plot
comm.pair.complete = comm.pair.test + t(comm.pair.test)
pair.test2plot = melt(comm.pair.complete, na.rm = TRUE)
pair.test2plot$value = round(pair.test2plot$value, 3)
text2plot = pair.test2plot

groupMat = matrix(NA, nrow = length(countries), ncol = length(countries))
colnames(groupMat) = countries
rownames(groupMat) = countries
groupMat[1:3,1:3] = 'Asia'
groupMat[4:6,4:6] = 'EU'
groupMat[7:8,7:8] = 'NA'
groupMat = melt(groupMat, na.rm = T)
groupMat$value = as.factor(groupMat$value)

# generate tikz file
tikz('output/Plots/tikz/pvalMat.tikz', standAlone = F, width = 5, height = 4.5)
ggplot(data = pair.test2plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black",
                       limit = c(0,1), space = "Lab",
                       name="p-value:") +
  theme_minimal()+
  geom_text(data = text2plot,
            aes(Var2, Var1, label = sprintf("%0.3f", round(value, digits = 3)), color = value > 0.5),
            size = 3.5)+
  scale_color_manual(guide = 'none', values = c("black", "white")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = 'right',#c(0.5, 0.75),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(0, 0.5),
        legend.direction = 'vertical')+#"horizontal")+
  coord_fixed() + ggtitle('Simultaneous commutator test p-values') +
  scale_y_discrete(position = 'left') +scale_x_discrete(limits = rev(levels(pair.test2plot$Var2))) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 7, title.vjust = 1,
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(margin = margin(0,0,8,0))))+
  geom_tile(data = groupMat, aes(x = Var1, y = Var2), colour = "red", fill = NA, linewidth = 1)
dev.off()

# generate png file
p2 = ggplot(data = pair.test2plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black",
                       limit = c(0,1), space = "Lab",
                       name="p-value:") +
  theme_minimal()+
  geom_text(data = text2plot, aes(Var2, Var1, label = sprintf("%0.3f", round(value, digits = 3)), color = value > 0.5), size = 3.5)+
  scale_color_manual(guide = 'none', values = c("black", "white")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = 'right',#c(0.5, 0.75),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(0, 0.5),
        legend.direction = 'vertical')+#"horizontal")+
  coord_fixed() + ggtitle('Simultaneous commutator test p-values') +
  scale_y_discrete(position = 'left') +scale_x_discrete(limits = rev(levels(pair.test2plot$Var2))) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 7, title.vjust = 1,
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(margin = margin(0,0,8,0))))+
  geom_tile(data = groupMat, aes(x = Var1, y = Var2), colour = "red", fill = NA, linewidth = 1)
p2
ggsave(filename = 'output/Plots/png/pvalMat.png', width = 5, height = 4.5, units = 'in')


