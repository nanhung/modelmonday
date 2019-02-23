library(PMXStan)

m <- PMXStanModel(path = "models")
stan_model("models/popPK_2cmpt_ivinfs_clearance_cls.stan")
dat <- prepareInputData(data.file = Theoph, model = m)

PMXStanFit()