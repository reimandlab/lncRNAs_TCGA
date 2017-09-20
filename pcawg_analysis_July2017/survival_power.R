res.power <- power.stratify(
  n = 70, 
  timeUnit = 1.25, 
  gVec = c(0.5, 0.5),
  PVec = c(0.5, 0.5), 
  HR = 1 / 1.91, 
  lambda0Vec = c(2.303, 1.139),
  power.ini = 0.8, 
  power.low = 0.001, 
  power.upp = 0.999,
  alpha = 0.05, 
  verbose = TRUE
  )


library(powerSurvEpi)

lncsTostudy <- c("LINC00665", "LINC00657")
lncsTostudy <- lnc_rna[,c(which(colnames(lnc_rna) %in% lncsTostudy),216:220)]

#only ovarian
lncsTostudy <- as.data.table(lncsTostudy)
lncsTostudy <- filter(lncsTostudy, canc == "Ovary Serous cystadenocarcinoma")

lncsTostudy$status[lncsTostudy$status=="deceased"] = 1
lncsTostudy$status[lncsTostudy$status=="alive"] = 0
lncsTostudy$time <- as.numeric(lncsTostudy$time)
lncsTostudy$status <- as.numeric(lncsTostudy$status)
lncsTostudy[,2] <- log1p(lncsTostudy[,2])

med1 <- median(lncsTostudy$LINC00665)
lncsTostudy$LINC00665[lncsTostudy$LINC00665 >=med1] = "E"
lncsTostudy$LINC00665[lncsTostudy$LINC00665 <med1] = "C"

lncsTostudy$LINC00665 <- as.factor(lncsTostudy$LINC00665)

res <- powerCT(formula = Surv(time, status) ~ LINC00665, dat = lncsTostudy,
nE = 35, nC = 35, RR = 2, alpha = 0.05)

