#===============================================================================
# LIBRARIES
#===============================================================================
for(pack in c("lme4", "lmerTest")) {
  if(!(pack %in% installed.packages()[,1])) {
    install.packages(pack)  
  }
  library(pack, character.only = T)
}


#===============================================================================
# HELPER FUNCTION
#===============================================================================
# function for displaying results
results <- function(model, effect) {
  coefs <- summary(model)$coefficients      

  ES <- fixef(model)[effect]
  CI <- confint.merMod(model, effect, method = "Wald")

  c(sprintf("%.2f", ES), 
    paste0("[", sprintf("%.2f", CI[1]), ", ",
           sprintf("%.2f", CI[2]), "]")
    )
}


#===============================================================================
# ANALYSIS
#===============================================================================
d <- read.csv("edge1.1.csv", head = T, as.is = T)
remove <- d$UniqueContributors < 2 & d$Type == 2
d <- d[!remove,]
d$Role <- factor(d$Role)
d <- cbind(d, ExtremeStatus = abs(3.5 - d$AcademicHierarchyStrict))

mccall <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))

qd <- d[d$Type == 1,]
qd <- cbind(qd, verbosity = mccall(qd$WC))

cd <- d[d$Type == 2 & d$WC > 50,]
cd <- cbind(cd, dominance = mccall(mccall(cd$we) - mccall(cd$i)))
cd <- cbind(cd, verbosity = mccall(cd$WC))
cd <- cbind(cd, femaleproportion = cd$UniqueFemaleContributors/cd$UniqueContributors)


# Hypothesis 1a
m1a <- lmer(verbosity ~ AcademicHierarchyStrict + (1|ThreadId) + (1|Id), data = qd)
summary(m1a)
  
# Hypothesis 1b
m1b <- lmer(dominance ~ AcademicHierarchyStrict + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m1b)

# Hypothesis 1c
m1c <- lmer(verbosity ~ Male + (1|ThreadId) + (1|Id), data = qd)
summary(m1c)

# Hypothesis 1d
m1d <- lmer(dominance ~ Male + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m1d)

# Hypothesis 2a
m2a <- lmer(verbosity ~ Male*AcademicHierarchyStrict + (1|ThreadId) + (1|Id), data = qd)
summary(m2a)

# Hypothesis 2b
m2b <- lmer(dominance ~ Male*AcademicHierarchyStrict + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m2b)

# Hypothesis 3a
m3a <- lmer(verbosity ~ ExtremeStatus + (1|ThreadId) + (1|Id), data = qd)
summary(m3a)

# Hypothesis 3b
m3b <- lmer(dominance ~ ExtremeStatus + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m3b)

# Hypothesis 4
m4 <- lmer(verbosity ~ Female*femaleproportion + (1|ThreadId) + (1|Id), data = cd)
summary(m4)

# Hypothesis 5a
m5a <- lmer(verbosity ~ Male*Live + (1|ThreadId) + (1|Id), data = cd)
summary(m5a)

# Hypothesis 5b
m5b <- lmer(dominance ~ Male*Live + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m5b)

# Additional Hypothesis 1.1
qd <- cbind(qd, dominance = mccall(mccall(qd$we) - mccall(qd$i)))
m6a <- lmer(dominance ~ AcademicHierarchyStrict + (1|ThreadId) + (1|Id), data = qd)
summary(m6a)

# Additional Hypothesis 1.2
m6b <- lmer(dominance ~ Female + (1|ThreadId) + (1|Id), data = qd)
summary(m6b)

# Additional Hypothesis 1.3
m6c <- lmer(dominance ~ AcademicHierarchyStrict + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m6c)

# Additional Hypothesis 1.4


# Additional Hypothesis 3.1
m8a <- lmer(Exclam ~ AcademicHierarchyStrict + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m8a)

# Additional Hypothesis 3.2
m8b <- lmer(Exclam ~ Female + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m8b)


