library(lme4)


d <- read.csv("edge1.1.csv", head = T, as.is = T)
simple <- read.csv("edge_simple1.1.csv", head = T, as.is = T)

remove <- d$UniqueContributors < 2 & d$Type == 2
d <- d[!remove,]
d$Role <- factor(d$Role)


# Hypothesis 1a
m1a.Q <- lmer(Number.Characters ~ AcademicHierarchyStrict + (1|Id), data = d[d$Type == 1,])
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




