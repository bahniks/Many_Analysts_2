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

mccall <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
qd <- d[d$Type == 1,]
cd <- d[d$Type == 2 & d$Role == 2,]
qd <- cbind(qd, mccall.WC = mccall(qd$WC))


# Hypothesis 1a
m1a.Q <- lmer(mccall.WC ~ AcademicHierarchyStrict + + (1|ThreadId) + (1|Id), 
              data = d[d$Type == 1,])
summary(m1a.Q)
m1a.C <- lmer(Number.Characters ~ AcademicHierarchyStrict + UniqueContributors + (1|Id), 
              data = d[d$Type == 2 & d$Role == 2,])
summary(m1a.C)
  
# Hypothesis 1b
m1b.iQ <- lmer(i ~ AcademicHierarchyStrict + (1|Id), data = d[d$Type == 1,])
summary(m1b.iQ)
m1b.weQ <- lmer(we ~ AcademicHierarchyStrict + (1|Id), data = d[d$Type == 1,])
summary(m1b.weQ)

m1b.iC <- lmer(i ~ AcademicHierarchyStrict*Role + (1|ThreadId) + (1|Id), 
               data = d[d$Type == 2 & d$WC > 10,])
summary(m1b.iC)
m1b.weC <- lmer(we ~ AcademicHierarchyStrict*Role + (1|ThreadId) + (1|Id), 
                data = d[d$Type == 2 & d$WC > 10,])
summary(m1b.weC)


# Hypothesis 1c
m1c.Q <- lmer(Number.Characters ~ Female + (1|Id), data = d[d$Type == 1,])
summary(m1c.Q)
m1c.C <- lmer(Number.Characters ~ Female + UniqueContributors + (1|Id), 
              data = d[d$Type == 2 & d$Role == 2,])
summary(m1c.C)

# Hypothesis 1d
m1b.iQ <- lmer(i ~ Female + (1|Id), data = d[d$Type == 1,])
summary(m1b.iQ)
m1b.weQ <- lmer(we ~ Female + (1|Id), data = d[d$Type == 1,])
summary(m1b.weQ)

m1b.iC <- lmer(i ~ Female*Role + (1|ThreadId) + (1|Id), 
               data = d[d$Type == 2 & d$WC > 10,])
summary(m1b.iC)
m1b.weC <- lmer(we ~ Female*Role + (1|ThreadId) + (1|Id), 
                data = d[d$Type == 2 & d$WC > 10,])
summary(m1b.weC)




