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
# HELPER FUNCTIONS
#===============================================================================
# formating p-values
round_p <- function(x) {
  if(x < .001) return("p < .001") 
  if(round(x, 2) == 1) return("p = 1")
  
  if(round(x, 3) >= .01) {
    p <- sprintf("%.2f", x)
  } else {
    p <- sprintf("%.3f", x)
  }
  
  return(paste("p = ", substr(p, 2, nchar(p)), sep = ""))
}

# function for displaying results
results <- function(model, effect, family = "linear") {
  coefs <- summary(model)$coefficients
  if(family == "binomial") {
    z <- coefs[effect, "z value"]
    p <- coefs[effect, "Pr(>|z|)"]
    e.name <- "OR"
    par <- "z"
  } else {
    z <- coefs[effect, "t value"]
    p <- coefs[effect, "Pr(>|t|)"]        
    e.name <- "b"
    par <- paste0("t(", sprintf("%.1f", coefs[effect, "df"]), ")")
  }
  ES <- fixef(model)[effect]
  CI <- confint.merMod(model, effect, method = "Wald")
  if(family == "binomial") {
    ES <- exp(ES)
    CI <- exp(CI)
  }
  cat(par, " = ", sprintf("%.2f", z), ", ", 
      round_p(p), ", ",
      e.name, " = ", sprintf("%.2f", ES), ", ",
      "95% CI = [", sprintf("%.2f", CI[1]), ", ",
      sprintf("%.2f", CI[2]), "]",
      sep = "")
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

dd <- d[sapply(unique(d$Title), function(x) which(d$Title == x)[1]),]
dd <- dd[dd$Type == 2 & dd$UniqueMaleContributors > 0 & dd$UniqueFemaleContributors > 0,]
dd <- cbind(dd, femaleproportion = dd$UniqueFemaleContributors/dd$UniqueContributors)
dd <- cbind(dd, overrepresented = (dd$FemaleParticipation > dd$femaleproportion))


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
m4 <- lmer(verbosity ~ Male*femaleproportion + (1|ThreadId) + (1|Id), data = cd)
summary(m4)

# Hypothesis 5a
m5a <- lmer(verbosity ~ Male*Live + (1|ThreadId) + (1|Id), data = cd)
summary(m5a)

# Hypothesis 5b
m5b <- lmer(dominance ~ Male*Live + Role + (1|ThreadId) + (1|Id), data = cd)
summary(m5b)
