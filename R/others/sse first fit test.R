test <- fit_model_sse(default_parameter_list = parameter_list,
                     parms_to_change = list(b1 = 0.07),
                     scenario = "vacc")


inits_test <- c(b1 = 0.05)
tic()
optim(fn=fit_model_sse, par = inits_test, default_parameter_list = parameter_list,
      parms_to_change = NULL, method = "Brent", lower = 0.05, upper = 0.07)
toc()


inits_test <- c(b1 = 0.05, b2 = 0.001, b3 = 0.001)
tic()
optim(fn=fit_model_sse, par = inits_test, default_parameter_list = parameter_list,
      parms_to_change = NULL, method = "Brent", lower = 0.05, upper = 0.07)
toc()



test <- run_model(default_parameter_list = parameter_list,
                      parms_to_change = list(b1 = 0.07),
                      scenario = "vacc")



