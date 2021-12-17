if (!"BayesCACE" %in% rownames(installed.packages()))
  install.packages("BayesCACE")
library("BayesCACE")

data("epidural_c", package = "BayesCACE")
epidural_c

data("epidural_ic", package = "BayesCACE")
head(epidural_ic)

plt.noncomp(epidural_c, overall = TRUE)

set.seed(123)
out.study <- cace.study(data = epidural_c, conv.diag = TRUE, mcmc.samples = TRUE, two.step = TRUE)
out.study$CACE
out.study$conv.out[[1]]
out.study$meta

out.meta.c <- cace.meta.c(data = epidural_c, conv.diag = TRUE, mcmc.samples = TRUE, study.specific = TRUE)
out.meta.c$smry
out.meta.c$DIC

out.meta.ic <- cace.meta.ic(data = epidural_ic, conv.diag = TRUE, mcmc.samples = TRUE, study.specific = TRUE)

plt.trace(obj = out.meta.ic)
plt.density(obj = out.meta.ic)
plt.acf(obj = out.meta.ic)

plt.forest(data = epidural_ic, obj = out.meta.ic)
plt.forest(data = epidural_c, obj = out.study, obj2 = out.meta.c)
