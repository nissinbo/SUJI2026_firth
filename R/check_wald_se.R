# Diagnose saved outputs, without refitting models.
z <- readRDS('results/model_fits.rds')
m <- read.csv('data/tables.csv')
m <- m[m$design_id=='E050R1' & m$data_state=='ok', ]
get <- function(key) {q<-z[z$key==key, ]; q[match(m$config_id,q$config_id), ]}
s <- get('sas_proc|fisher_default|Wald')
l <- get('logistf|lf_nr_default|Wald')
b <- get('brglm2|bg_asmean_default|Wald')
stopifnot(identical(s$config_id,l$config_id), identical(s$config_id,b$config_id))
p1 <- (m$a+.5)/26
p0 <- (m$c+.5)/26
beta <- qlogis(p1)-qlogis(p0)
se <- sqrt(1/(25*p1*(1-p1))+1/(25*p0*(1-p0)))
se_aug <- se*sqrt(25/26)
outside <- function(lo,hi) lo>0|hi<0
sd <- outside(s$low,s$high)
ld <- outside(l$low,l$high)
theory_s <- abs(beta)>qnorm(.975)*se
theory_l <- abs(beta)>qnorm(.975)*se_aug
bad <- sd!=ld
out <- data.frame(config_id=m$config_id,exposed_events=m$a,unexposed_events=m$c,
  beta_sas=s$beta,beta_logistf=l$beta,se_sas=s$se,se_logistf=l$se,
  low_sas=s$low,high_sas=s$high,low_logistf=l$low,high_logistf=l$high,
  exact_beta=beta,ordinary_se=se,augmented_se=se_aug,different=bad,
  ordinary_outside=sd,augmented_outside=ld)
stopifnot(sum(bad)==12L,identical(theory_s,sd),identical(theory_l,ld),
          all(!sd[bad]),all(ld[bad]),all(l$se[bad]<s$se[bad]))
# Holding PROC LOGISTIC's point estimate fixed, use only logistf's SE.
swapped <- abs(s$beta)>qnorm(.975)*l$se
stopifnot(identical(swapped,ld))
# In the other direction, keep logistf's beta and use PROC LOGISTIC's SE.
swapped_to_sas <- abs(l$beta)>qnorm(.975)*s$se
stopifnot(identical(swapped_to_sas,sd))
root <- 'results/tables'
write.csv(out,file.path(root,'wald_se.csv'),row.names=FALSE)
record <- c('Saved default Wald fits; individual binary logistic model, exposure only, 25 per group.',
  'At common coefficients h_i=1/25, I_aug=(26/25) I, SE_aug=sqrt(25/26) SE.',
  sprintf('SE ratio: %.12f; reduction percent: %.8f',sqrt(25/26),100*(1-sqrt(25/26))),
  sprintf('Disagreements: %d/674; all ordinary CI includes zero and logistf CI excludes zero.',sum(bad)),
  'Using identical closed-form Firth coefficients and the two SE formulas reproduces all 674 decisions.',
  'Holding PROC LOGISTIC beta fixed and replacing only SE with logistf SE reproduces all 674 logistf decisions.',
  'Holding logistf beta fixed and replacing only SE with PROC LOGISTIC SE reproduces all 674 PROC LOGISTIC decisions.',
  'Aligning SE to either saved implementation gives zero decision disagreements across all 674 tables; no models were refitted.',
  sprintf('Max absolute beta difference among 12 disagreements: %.12g',max(abs(s$beta[bad]-l$beta[bad]))),
  sprintf('Max relative SE formula error, all tables: SAS %.12g; logistf %.12g',max(abs(s$se/se-1)),max(abs(l$se/se_aug-1))))
writeLines(record,'results/checks/wald_se.txt')
cat(paste(record,collapse='\n'),'\n')
print(out[bad,1:11],row.names=FALSE,digits=7)
