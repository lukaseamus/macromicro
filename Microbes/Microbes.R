# 1. Load data ####
# Load output from BMG CLARIOstar Plus microplate reader
require(tidyverse)
require(here)
files <- list.files(path = here("Microbes", "Raw"), pattern = "\\.csv$", full.names = TRUE)

enzyme <- files %>%
  map(~ read.csv(., skip = 9, header = TRUE) %>%
        rename(Fluorescence = "Raw.Data..360.20.450.30.") %>%
        mutate(Fluorescence = as.numeric(Fluorescence))) %>%
  set_names(str_remove(basename(files), "\\..*$") %>% make.names) %>%
  imap(~ .x %>% mutate(Date = str_extract(.y, "\\d+(?=enzyme)") %>% dmy()))

str(enzyme)

# Load ordered metadata
meta <- here("Microbes", "Meta.csv") %>% read.csv()

# Join annotations to respective measurements
require(magrittr)
enzyme %<>% 
  imap(~ bind_cols(.x, meta %>% filter(Plate == str_remove(.y, "^X"))) %>%
         mutate(Date2 = str_extract(Plate, "\\d+(?=enzyme)") %>% dmy()))

str(enzyme)

# Cross-validate that join worked correctly
enzyme %>% 
  imap(~ .x %>% 
         select(Annotation) %>%
         identical(meta %>%
                     filter(Plate == str_remove(.y, "^X")) %>%
                     select(Annotation))
  )
# or
enzyme %>% 
  map(~ .x %$% identical(Date, Date2))

# remove additional Date variable and change Plate to number
enzyme %<>% 
  map(~ .x %>% select(-Date2) %>%
        mutate(Plate = str_extract(Plate, "\\d+(?=[^\\d]*$)") %>%
                 as.numeric()))
str(enzyme)

# The project ran in multiple phases, and I will work through them in order 
# 1. on 19th and 20th October 2021 enzyme kinetics were tested
# 2. on 3rd November I ran several fluorescence tests with ultrapure water and seawater relative to methanol
# 3. from 28th October I ran samples and after 3rd November consistently worked with seawater blanks and controls 

# 2. Enzyme kinetics ####
kinetics <- enzyme %>%
  map(~ .x %>% filter(Date <= "201021" %>% dmy())) %>%
  keep(~ nrow(.x) > 0)

# Split each dataframe in kinetics into its components
kinetics %<>% map(~ list(
  standard = .x %>% # the = instead of <- is important to retain names
    filter(str_detect(Annotation, "^[0-9.]+$")) %>% 
    rename(Concentration = Annotation) %>%
    mutate(Concentration = as.numeric(Concentration)) %>%
    filter(between(Concentration, 0, 100)), # concentrations higher than 100 µM show a decline in fluorescence
  samples = .x %>% filter(!str_detect(Annotation, "^[0-9.]+$"))
))

kinetics$X191021enzyme1$standard # examples
kinetics$X191021enzyme1$samples

# 2.1 Standard curve models ####
# Visualise
require(patchwork)
kinetics %>%
  map(~ .x$standard %>%
        ggplot(aes(Concentration, Fluorescence)) +
          geom_point() +
          geom_line() +
          # coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 2e3)) + # unhash to see intercept
          theme_minimal()
      ) %>%
  wrap_plots()
# In the past I have attempted to approximate this trend with polynomials for simplicity.
# However, polynomials are not easy to work with mathematically, and once the cause of
# non-linearity is understood (inner filter effect), saturation models are the obvious choice.
# There are several non-linear functions to choose from that describe a saturation effect with
# only two parameters: beta (linear slope at low values) and Fmax (the maximal fluorescence) 
# (alternatively these can be parameterised with the half saturation constant K). 

# The simplest are the rectangular hyperbola (i.e. Michaelis-Menten function)
# F = Fmax * beta * c / ( Fmax + beta * c )
# the exponential saturation function (i.e. exponential cumulative distribution function)
# F = Fmax * ( 1 - exp( -beta * c / Fmax ) )
# and the hyperbolic tangent
# F = Fmax * tanh( beta * c / Fmax )

# However, even water has some fluorescence, so an intercept term F0 must be added.
# I instinctively favour the last function because it "sticks" to the linear phase for longest:
Fmax <- 24e4 # example parameter values
beta <- 5e3
F0 <- 500
tibble(c = seq(0, 100)) %>%
  mutate(lm = F0 + beta * c, # linear model for comparison
         rh = Fmax * beta * c / ( Fmax + beta * c ) + F0,
         es = Fmax * ( 1 - exp( -beta * c / Fmax ) ) + F0,
         ht = Fmax * tanh( beta * c / Fmax ) + F0) %>%
  pivot_longer(cols = -c, names_to = "Function", values_to = "F") %>%
  ggplot(aes(c, `F`, colour = Function)) +
    geom_hline(yintercept = 24e4) + # Fmax
    geom_line() +
    theme_minimal() +
    coord_cartesian(ylim = c(0, 24e4))

# 2.1.1 Prior simulation ####
# A sensible prior for Fmax seems to be that chosen above because I never
# recorded fluorescence values above that on the BMG CLARIOstar.
# For beta the sensible value also seems to be that chosen above
# since based on visualisation values range between around 2000 and 8000.
# Baseline fluorescence is usually below 1000 relative to maximal fluorescence,
# so F0 is best centred around 500 as above. All parameters have to be positive,
# so the gamma distribution is best.

kinetics_standard_ht_prior <-
  tibble(n = 1:1e3,
         Fmax = rgamma(n = 1e3, shape = Fmax^2 / 5e4^2, rate = Fmax / 5e4^2), # reparameterised in terms of mean and s.d.
         beta = rgamma(n = 1e3, shape = beta^2 / 3e3^2, rate = beta / 3e3^2),
         F0 = rgamma(n = 1e3, shape = F0^2 / 300^2, rate = F0 / 300^2)) %>%
  expand_grid(c = seq(0, 100)) %>%
  mutate(`F` = Fmax * tanh( beta * c / Fmax ) + F0)

kinetics_standard_ht_prior %>%
  ggplot(aes(c, `F`, group = n)) +
    geom_hline(yintercept = 24e4) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off",
                    # xlim = c(0, 1), ylim = c(0, 2e3) # unhash to check F0
                    ) +
    theme_minimal()
# covers all reasonable possibilities

# 2.1.2 Run model ####
kinetics_standard_ht_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> Fmax;
  real<lower=0> beta;
  real<lower=0> F0;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmax ~ gamma( 24e4^2 / 5e4^2, 24e4 / 5e4^2 );
  beta ~ gamma( 5e3^2 / 3e3^2, 5e3 / 3e3^2 );
  F0 ~ gamma( 500^2 / 300^2, 500 / 300^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmax * tanh( beta * Concentration[i] / Fmax ) + F0;
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

require(cmdstanr)
kinetics_standard_ht_mod <- kinetics_standard_ht_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
kinetics_standard_ht_samples <- kinetics %>%
  map(~ kinetics_standard_ht_mod$sample(data = .x$standard %>%
                                            select(Fluorescence, Concentration) %>%
                                            compose_data(),
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)
      )

# 2.1.3 Model checks ####
# check Rhat, effective sample size and chains
kinetics_standard_ht_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
kinetics_standard_ht_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_rank_overlay()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_ht_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# some correlation between Fmax and beta, indicating some interdependence

# 2.1.4 Prior-posterior comparison ####
source("functions.R")
# sample prior
kinetics_standard_ht_prior <- kinetics %>%
  map(~ prior_samples(model = kinetics_standard_ht_mod,
                      data = .x$standard %>%
                        select(Fluorescence, Concentration) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
kinetics_standard_ht_prior %>%
  map2(kinetics_standard_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA), # no groups so this has to be an empty list or tibble
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "long")) %>%
  map(~ .x %>% prior_posterior_plot()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# some posteriors for F0 broke out of the expected prior probability space

# 2.1.5 Predictions ####
kinetics_standard_ht_predictions <- kinetics_standard_ht_prior %>%
  map2(kinetics_standard_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "short")) %>%
  map2(kinetics, ~ spread_continuous(.x, .y$standard, predictor_name = "Concentration")) %>%
  map(~ .x %>% mutate(mu = Fmax * tanh( beta * Concentration / Fmax ) + F0,
                      obs = rnorm(n(), mu, sigma)))

kinetics_standard_ht_predictions_summary <- kinetics_standard_ht_predictions %>%
  map(~ .x %>% group_by(distribution, Concentration) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_")) 

kinetics_standard_ht_predictions_summary %>%
  map2(kinetics,
       ~ ggplot() +
            geom_point(data = .y$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
                            alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots()
# fit doesn't look as good as expected
# try exponential saturation function 

# 2.1.6 Prior simulation ####
kinetics_standard_es_prior <-
  tibble(n = 1:1e3,
         Fmax = rgamma(n = 1e3, shape = Fmax^2 / 5e4^2, rate = Fmax / 5e4^2),
         beta = rgamma(n = 1e3, shape = beta^2 / 3e3^2, rate = beta / 3e3^2),
         F0 = rgamma(n = 1e3, shape = F0^2 / 300^2, rate = F0 / 300^2)) %>%
  expand_grid(c = seq(0, 100)) %>%
  mutate(`F` = Fmax * ( 1 - exp( -beta * c / Fmax ) ) + F0)

kinetics_standard_es_prior %>%
  ggplot(aes(c, `F`, group = n)) +
    geom_hline(yintercept = 24e4) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off",
                    # xlim = c(0, 1), ylim = c(0, 2e3) # unhash to check F0
    ) +
    theme_minimal()

# 2.1.7 Run model ####
kinetics_standard_es_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> Fmax;
  real<lower=0> beta;
  real<lower=0> F0;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmax ~ gamma( 24e4^2 / 5e4^2, 24e4 / 5e4^2 );
  beta ~ gamma( 5e3^2 / 3e3^2, 5e3 / 3e3^2 );
  F0 ~ gamma( 500^2 / 300^2, 500 / 300^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmax * ( 1 - exp( -beta * Concentration[i] / Fmax ) ) + F0;
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

kinetics_standard_es_mod <- kinetics_standard_es_stan %>%
  write_stan_file() %>%
  cmdstan_model()

kinetics_standard_es_samples <- kinetics %>%
  map(~ kinetics_standard_es_mod$sample(data = .x$standard %>%
                                          select(Fluorescence, Concentration) %>%
                                          compose_data(),
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)
      )

# 2.1.8 Model checks ####
# check Rhat, effective sample size and chains
kinetics_standard_es_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

kinetics_standard_es_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_rank_overlay()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_es_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# some correlation between Fmax and beta, indicating some interdependence

# 2.1.9 Prior-posterior comparison ####
# sample prior
kinetics_standard_es_prior <- kinetics %>%
  map(~ prior_samples(model = kinetics_standard_es_mod,
                      data = .x$standard %>%
                        select(Fluorescence, Concentration) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
kinetics_standard_es_prior %>%
  map2(kinetics_standard_es_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "long")) %>%
  map(~ .x %>% prior_posterior_plot()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# some posteriors for F0 broke out of the expected prior probability space

# 2.1.10 Predictions ####
kinetics_standard_es_predictions <- kinetics_standard_es_prior %>%
  map2(kinetics_standard_es_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "short")) %>%
  map2(kinetics, ~ spread_continuous(.x, .y$standard, predictor_name = "Concentration")) %>%
  map(~ .x %>% mutate(mu = Fmax * ( 1 - exp( -beta * Concentration / Fmax ) ) + F0,
                      obs = rnorm(n(), mu, sigma)))

kinetics_standard_es_predictions_summary <- kinetics_standard_es_predictions %>%
  map(~ .x %>% group_by(distribution, Concentration) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_")) 

kinetics_standard_es_predictions_summary %>%
  map2(kinetics,
       ~ ggplot() +
            geom_point(data = .y$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
                            alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots()
# exponential saturation function fits better but is still not perfect
# try rectangular hyperbola

# 2.1.11 Prior simulation ####
kinetics_standard_rh_prior <-
  tibble(n = 1:1e3,
         Fmax = rgamma(n = 1e3, shape = Fmax^2 / 5e4^2, rate = Fmax / 5e4^2),
         beta = rgamma(n = 1e3, shape = beta^2 / 3e3^2, rate = beta / 3e3^2),
         F0 = rgamma(n = 1e3, shape = F0^2 / 300^2, rate = F0 / 300^2)) %>%
  expand_grid(c = seq(0, 100)) %>%
  mutate(`F` = Fmax * beta * c / ( Fmax + beta * c ) + F0)

kinetics_standard_rh_prior %>%
  ggplot(aes(c, `F`, group = n)) +
    geom_hline(yintercept = 24e4) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off",
                    # xlim = c(0, 1), ylim = c(0, 2e3) # unhash to check F0
    ) +
    theme_minimal()

# 2.1.12 Run model ####
kinetics_standard_rh_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> Fmax;
  real<lower=0> beta;
  real<lower=0> F0;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmax ~ gamma( 24e4^2 / 5e4^2, 24e4 / 5e4^2 );
  beta ~ gamma( 5e3^2 / 3e3^2, 5e3 / 3e3^2 );
  F0 ~ gamma( 500^2 / 300^2, 500 / 300^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmax * beta * Concentration[i] / ( Fmax + beta * Concentration[i] ) + F0;
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

kinetics_standard_rh_mod <- kinetics_standard_rh_stan %>%
  write_stan_file() %>%
  cmdstan_model()

kinetics_standard_rh_samples <- kinetics %>%
  map(~ kinetics_standard_rh_mod$sample(data = .x$standard %>%
                                          select(Fluorescence, Concentration) %>%
                                          compose_data(),
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)
      )

# 2.1.13 Model checks ####
# check Rhat, effective sample size and chains
kinetics_standard_rh_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

kinetics_standard_rh_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_rank_overlay()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_rh_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# some correlation between Fmax and beta, indicating some interdependence

# 2.1.14 Prior-posterior comparison ####
# sample prior
kinetics_standard_rh_prior <- kinetics %>%
  map(~ prior_samples(model = kinetics_standard_rh_mod,
                      data = .x$standard %>%
                        select(Fluorescence, Concentration) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
kinetics_standard_rh_prior %>%
  map2(kinetics_standard_rh_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "long")) %>%
  map(~ .x %>% prior_posterior_plot()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# F0 broke out of the expected prior probability space less

# 2.1.15 Predictions ####
kinetics_standard_rh_predictions <- kinetics_standard_rh_prior %>%
  map2(kinetics_standard_rh_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "short")) %>%
  map2(kinetics, ~ spread_continuous(.x, .y$standard, predictor_name = "Concentration")) %>%
  map(~ .x %>% mutate(mu = Fmax * beta * Concentration / ( Fmax + beta * Concentration ) + F0,
                      obs = rnorm(n(), mu, sigma)))

kinetics_standard_rh_predictions_summary <- kinetics_standard_rh_predictions %>%
  map(~ .x %>% group_by(distribution, Concentration) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

kinetics_standard_rh_predictions_summary %>%
  map2(kinetics,
       ~ ggplot() +
            geom_point(data = .y$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
                            alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots()
# rectangular hyperbola is clearly the best!
# it only slightly overestimates the intercept when beta is low
# and underestimates it when beta is high, but that is less worrisome 
# because I will adjust the intercept based on water autofluorescence

# 2.1.16 Compare functions ####
kinetics_standard_predictions_summary <- 
  kinetics_standard_ht_predictions_summary %>%
  map2(kinetics_standard_es_predictions_summary, ~ bind_rows(.x, .y)) %>%
  map2(kinetics_standard_rh_predictions_summary, ~ bind_rows(.x, .y) %>%
         mutate(Function = c(rep("ht", nrow(.y)), 
                             rep("es", nrow(.y)),
                             rep("rh", nrow(.y))) %>% 
                  fct())
       )

kinetics_standard_predictions_summary %>%
  map2(kinetics,
       ~ ggplot() +
            geom_point(data = .y$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .x %>% filter(distribution == "posterior"), 
                      aes(Concentration, mu_y, colour = Function)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"), 
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width), fill = Function)) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# yes, rectangular hyperbola is clearly best at all scales
# and is chosen as the optimal model

# 2.2 Convert sample fluorescence ####
# solving F = Fmax * beta * c / ( Fmax + beta * c ) + F0 for c
# gives c = Fmax * ( F - F0 ) / ( beta * ( F0 + Fmax - F ) )

kinetics_samples <- kinetics %>%
  map2(kinetics_standard_rh_prior %>%
         map2(kinetics_standard_rh_samples,
              ~ prior_posterior_draws(prior_samples = .x,
                                      posterior_samples = .y,
                                      group = list(NA),
                                      parameters = c("Fmax", "beta", "F0", "sigma"),
                                      format = "short") %>%
                filter(distribution == "posterior")),
       ~ .x$samples %>% 
         mutate(Concentration = str_extract(Annotation, "[0-9]+") %>% as.numeric(),
                Species = str_extract(Annotation, "[a-zA-Z]+") %>% fct()) %>%
         cross_join(.y) %>%
         mutate(Activity = ( Fmax * ( Fluorescence - F0 ) / ( beta * ( F0 + Fmax - Fluorescence ) ) )) 
       )
str(kinetics_samples)

kinetics_samples %>%
  map(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 0) +
          geom_violin(aes(Concentration, Activity, fill = Species, group = Well),
                      colour = NA, alpha = 0.5, position = "identity") +
          coord_cartesian(ylim = c(NA, 1)) +
          theme_minimal() +
          theme(panel.grid = element_blank())
  ) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# These data are exceptionally messy, likely due to bad pipetting, and not worth analysing.
# There was also no correction for the autofluorescence difference between water and methanol,
# so there are impossible negative values. There is no indication of saturation, with the highest 
# enzyme activity always found at the highest substrate concentration. I ended up deciding on 
# 1000 µM as the saturating substrate concentration due to logistical reasons.

# 3. Baseline fluorescence ####
enzyme$X031121enzyme9 %>%
  filter(Annotation %in% c("MQ", "AS", "0", "CMQ", "CAS")) %>%
  mutate(Annotation = if_else(Annotation == "0", "MeOH", Annotation) %>%
           fct_relevel("MQ", "AS", "CMQ")) %>%
  ggplot(aes(Annotation, Fluorescence)) +
    geom_point(shape = 16, alpha = 0.2, size = 3) +
    geom_point(data = . %>% group_by(Annotation) %>% 
                 summarise(Fluorescence = median(Fluorescence)),
               aes(Annotation, Fluorescence), size = 5, shape = 21, fill = NA) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Conclusions: 
# 1. MilliQ water (MQ) and autoclaved seawater (AS) have similar fluorescence properties.
# 2. The enzyme substrate has similar fluorescence properties in MilliQ (CMQ) and seawater (CAS).
# 3. 100% methanol (MeOH) has higher fluorescence than either water with or without substrate.
# It follows that
# 1. The standard curve, which is made up with 100% methanol, needs to be corrected using a water blank
# 2. The water reference can be either MilliQ or autoclaved seawater
# 3. The substrate control can be either MilliQ or autoclaved seawater

rm(list = setdiff(ls(), c("enzyme", "meta")))

# 4. Experimental samples ####
experimental <- enzyme %>%
  map(~ .x %>% filter(Date >= "281021" %>% dmy() & 
                        !(Date == "031121" %>% dmy() & Plate == 9))) %>%
  keep(~ nrow(.x) > 0)

# Split each dataframe in experimental into its components
experimental %<>% map(~ list(
  blank = .x %>% filter(Annotation %in% c("AS", "MQ")),
  control = .x %>% filter(Annotation %in% c("CAS", "CMQ")),
  standard = .x %>%
    filter(str_detect(Annotation, "^[0-9.]+$")) %>% 
    rename(Concentration = Annotation) %>%
    mutate(Concentration = as.numeric(Concentration)) %>%
    filter(between(Concentration, 0, 100)),
  samples = .x %>% filter(!str_detect(Annotation, "^[0-9.]+$") & 
                            !Annotation %in% c("MQ", "AS", "CMQ", "CAS"))
))

str(experimental)

experimental$X031121enzyme1$blank # examples
experimental$X031121enzyme1$control
experimental$X031121enzyme1$standard
experimental$X031121enzyme1$samples

# 4.1 Standard curve models ####
# Visualise
experimental %>%
  map(~ .x$standard %>%
        ggplot(aes(Concentration, Fluorescence)) +
          geom_point() +
          geom_line() +
          theme_minimal() +
          theme(panel.grid = element_blank())
  ) %>%
  wrap_plots() %>% # too large to view in RStudio
  ggsave(filename = "Plots/standard_data.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 4.1.1 Prior simulation ####
experimental_standard_prior <-
  tibble(n = 1:1e3,
         Fmax = rgamma(n = 1e3, shape = Fmax^2 / 5e4^2, rate = Fmax / 5e4^2),
         beta = rgamma(n = 1e3, shape = beta^2 / 3e3^2, rate = beta / 3e3^2),
         F0 = rgamma(n = 1e3, shape = F0^2 / 300^2, rate = F0 / 300^2)) %>%
  expand_grid(c = seq(0, 100)) %>%
  mutate(`F` = Fmax * beta * c / ( Fmax + beta * c ) + F0)

experimental_standard_prior %>%
  ggplot(aes(c, `F`, group = n)) +
    geom_hline(yintercept = 24e4) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off",
                    # xlim = c(0, 1), ylim = c(0, 2e3) # unhash to check F0
    ) +
    theme_minimal()

# 4.1.2 Run model ####
experimental_standard_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> Fmax;
  real<lower=0> beta;
  real<lower=0> F0;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmax ~ gamma( 24e4^2 / 5e4^2, 24e4 / 5e4^2 );
  beta ~ gamma( 5e3^2 / 3e3^2, 5e3 / 3e3^2 );
  F0 ~ gamma( 500^2 / 300^2, 500 / 300^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmax * beta * Concentration[i] / ( Fmax + beta * Concentration[i] ) + F0;
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

experimental_standard_mod <- experimental_standard_stan %>%
  write_stan_file() %>%
  cmdstan_model()

experimental_standard_samples <- experimental %>%
  map(~ experimental_standard_mod$sample(data = .x$standard %>%
                                          select(Fluorescence, Concentration) %>%
                                          compose_data(),
                                         chains = 8,
                                         parallel_chains = parallel::detectCores(),
                                         iter_warmup = 1e4,
                                         iter_sampling = 1e4)
      )

# 4.1.3 Model checks ####
# check Rhat, effective sample size and chains
experimental_standard_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

ggsave( # too large to view in RStudio
  experimental_standard_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "Plots/standard_rank.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

experimental_standard_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots() %>%
  ggsave(filename = "Plots/standard_pairs.pdf", width = 80, height = 80, unit = "cm", device = cairo_pdf)
# some correlation between Fmax and beta, indicating some interdependence

# 4.1.4 Prior-posterior comparison ####
# sample prior
source("functions.R")
experimental_standard_prior <- experimental %>%
  map(~ prior_samples(model = experimental_standard_mod,
                      data = .x$standard %>%
                        select(Fluorescence, Concentration) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  experimental_standard_prior %>%
    map2(experimental_standard_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA),
                                 parameters = c("Fmax", "beta", "F0", "sigma"),
                                 format = "long")) %>%
    map(~ .x %>% prior_posterior_plot()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "Plots/standard_prior_posterior.pdf", width = 80, height = 40,
  unit = "cm", device = cairo_pdf)
# F0 broke out quite a bit, so the model predicts higher intercepts than expected

# 4.1.5 Predictions ####
experimental_standard_prior_posterior <- experimental_standard_prior %>%
  map2(experimental_standard_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "short"))

experimental_standard_predictions <- experimental_standard_prior_posterior %>%
  map2(experimental, ~ spread_continuous(.x, .y$standard, predictor_name = "Concentration",
                                         length = 50)) %>%
  map(~ .x %>% mutate(mu = Fmax * beta * Concentration / ( Fmax + beta * Concentration ) + F0,
                      obs = rnorm(n(), mu, sigma)))

experimental_standard_predictions_summary <- experimental_standard_predictions %>%
  map(~ .x %>% group_by(distribution, Concentration) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

experimental_standard_predictions_summary %>%
  map2(experimental,
       ~ ggplot() +
            geom_point(data = .y$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
                            alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "Plots/standard_predictions.pdf", width = 80, height = 40,
         unit = "cm", device = cairo_pdf)
# great fit except for those deviations in F0

# 4.2 Blank models ####
# To correctly calculate enzyme activities of experimental samples, F0 of the standard curve
# needs to be replaced with a water blank. Blanks were measured in triplicate or more replicates
# so in order to preserve the variability, intercept models must be built. The prediction for mu
# of these intercept models can then be substituted for F0.

# Visualise
experimental %>%
  map(~ .x$blank %>%
        ggplot(aes(Well, Fluorescence)) +
          geom_point(size = 3, shape = 16) +
          geom_hline(data = . %>% summarise(Fluorescence = median(Fluorescence)),
                     aes(yintercept = Fluorescence)) +
          theme_minimal()
  ) %>%
  wrap_plots() %>% # too large to view in RStudio
  ggsave(filename = "Plots/blank_data.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# There are four plates missing blanks, which will have to be substituted with an average.
# This is best done with hierarchical modelling, which requires that the data be combined.
# Also data are clearly right-skewed due to their non-negative nature, so a gamma likelihood is best.

# 4.2.1 Prior simulation ####
# Water autofluorescence is around 200 arbitrary units relative to maximal 100 µM standard fluorescence.
# Since fluorescence cannot be negative, a gamma distribution is best as the likelihood. To properly 
# simulate the outcome we need to transform the normal parameter distribution by exponentiating. While
# the mean can simply be said to the logarithm of the expected mean, the sd must be simulated.

ggplot() +
  geom_density(aes(rnorm(1e5, log(200), 1) %>% exp())) +
  # scale_x_continuous(limits = c(0, 1e3), oob = scales::oob_keep) + # unhash to zoom in on peak
  theme_minimal()
# lots of variability

# here is the underlying distribution
ggplot() +
  geom_density(aes(rnorm(1e5, log(200), 1))) +
  theme_minimal()

# 4.2.2 Run model ####
experimental_blank_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
  array[n] int Group;
  int n_Group;
}

parameters{
  real Fmu_mu; // Hyperparameters
  real<lower=0> Fmu_sigma;

  vector[n_Group] Fmu;

  real<lower=0> sigma;
}

model{
  // Hyperpriors
  Fmu_mu ~ normal( log(200) , 1 );
  Fmu_sigma ~ exponential( 1 );

  // Group-level prior
  Fmu ~ normal( Fmu_mu , Fmu_sigma );

  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( Fmu[Group[i]] ); // Exponential transformation
  }

  // Gamma likelihood (generalsied linear model)
  Fluorescence ~ gamma( mu^2 / sigma^2 , mu / sigma^2 );
}
"

experimental_blank_mod <- experimental_blank_stan %>%
  write_stan_file() %>%
  cmdstan_model()

experimental_blank_samples <- experimental_blank_mod$sample(
  data = experimental %>%
    imap(~ .x$blank %>% # here I am creating the unique plate grouping variable
          mutate(Group = .y %>% fct())) %>%
    keep(~ nrow(.x) > 0) %>% # drop all list items without data
    bind_rows() %>% # here I am combining the blank dataframes to compare across plates
    select(Fluorescence, Group) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4)

# 4.2.3 Model checks ####
# check Rhat, effective sample size and chains
experimental_blank_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

experimental_blank_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

# 4.2.4 Prior-posterior comparison ####
# sample prior
experimental_blank_prior <- prior_samples(
  model = experimental_blank_mod,
  data = experimental %>%
    imap(~ .x$blank %>%
           mutate(Group = .y %>% fct())) %>%
    keep(~ nrow(.x) > 0) %>%
    bind_rows() %>%
    select(Fluorescence, Group) %>%
    compose_data(),
  chains = 8, samples = 1e4)
# prior_samples() doesn't fully support multilevel models but
# for our visualisation purposes the messy samples are enough.

# plot prior-posterior comparison
experimental_blank_prior %>%
    prior_posterior_draws(posterior_samples = experimental_blank_samples,
                          group = experimental %>%
                            imap(~ .x$blank %>%
                                   mutate(Group = .y %>% fct())) %>%
                            keep(~ nrow(.x) > 0) %>%
                            bind_rows() %>%
                            select(Group),
                          parameters = c("Fmu[Group]", "Fmu_mu", "Fmu_sigma", "sigma"),
                          format = "long") %>%
  prior_posterior_plot(group_name = "Group")
# Note the fuzziness of the prior. Stan doesn't seem to be able to smoothly estimate
# hierarchical priors.
# For smooth prior distributions, generate in R:

experimental_blank_prior <-
  tibble(.chain = 1:8 %>% rep(each = 1e4),
         .iteration = 1:1e4 %>% rep(8),
         .draw = 1:8e4,
         Fmu_mu = rnorm(8 * 1e4, log(200) , 1),
         Fmu_sigma = rexp(8 * 1e4, 1),
         Fmu = rnorm(8 * 1e4, Fmu_mu, Fmu_sigma),
         sigma = rexp(8 * 1e4, 1)) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = ".variable", values_to = ".value") %>%
  mutate(rep = if_else(.variable == "Fmu", 24, 1),
         .variable = fct_relevel(.variable, "Fmu")) %>%
  uncount(rep) %>%
  arrange(.variable) %>%
  mutate(Group = experimental %>%
           imap(~ .x$blank %>%
                  mutate(Group = .y %>% fct())) %>%
           keep(~ nrow(.x) > 0) %>%
           bind_rows() %$% levels(Group) %>%
           rep(8 * 1e4) %>%
           c(NA %>% rep(3 * 8 * 1e4)))

experimental_blank_posterior <-
  experimental_blank_samples %>%
  recover_types(experimental %>%
                  imap(~ .x$blank %>%
                         mutate(Group = .y %>% fct())) %>%
                  keep(~ nrow(.x) > 0) %>%
                  bind_rows() %>% select(Group)) %>%
  gather_draws(Fmu[Group], Fmu_mu, Fmu_sigma, sigma) %>%
  ungroup()


experimental_blank_prior %>%
  mutate(distribution = "prior") %>%
  bind_rows(experimental_blank_posterior %>%
              mutate(distribution = "posterior")) %>%
  prior_posterior_plot(group_name = "Group")
# Not fuzzy. Everything fine.

# 4.2.5 Predictions ####
experimental_blank_predictions <-
  experimental_blank_samples %>%
  recover_types(experimental %>%
                  imap(~ .x$blank %>%
                         mutate(Group = .y %>% fct())) %>%
                  keep(~ nrow(.x) > 0) %>%
                  bind_rows() %>% select(Group)) %>%
  spread_draws(Fmu[Group], sigma) %>%
  mutate(Fmu = exp(Fmu),
         obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2)) %>%
  select(-sigma) %>%
  bind_rows(
    experimental_blank_samples %>%
      spread_draws(Fmu_mu, Fmu_sigma, sigma) %>%
      mutate(Fmu = exp( rnorm(n(), Fmu_mu, Fmu_sigma) ),
             obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2),
             Group = "unobserved") %>%
      select(-c(Fmu_mu, Fmu_sigma, sigma))
  ) %>%
  pivot_longer(cols = c("Fmu", "obs"),
               names_to = "level", values_to = "Fluorescence")

experimental_blank_predictions %>%
  ggplot(aes(Fluorescence, alpha = factor(level))) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0, 500), oob = scales::oob_keep) +
    scale_alpha_manual(values = c(0.6, 0.2)) +
    facet_wrap(~ Group, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# It makes most sense to me to use Fmu here rather than predicted observations,
# since we are replacing F0. That is while we want Fmu on the response scale, it is
# to substitute another parameter so does not need to pass through the likelihood.
experimental_blank_predictions %>%
  filter(level == "Fmu") %>%
  ggplot(aes(Fluorescence)) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0, 500), oob = scales::oob_keep) +
    facet_wrap(~ Group, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 4.3 Control models ####
# To correctly calculate enzyme activities of experimental samples, a substrate autogenic
# fluorescence control needs to be subtracted after conversion with the standard curve. If
# the substrate control were subtracted before conversion with the standard curve, water
# autofluorescence which is naturally contained within the substrate control would
# effectively be subtracted twice, once with the control and once during conversion with
# the standard curve which has the water blank as its intercept.

# Visualise
experimental %>%
  imap(~ .x$control %>%
        ggplot(aes(Well, Fluorescence)) +
          geom_point(size = 3, shape = 16) +
          geom_hline(data = . %>% summarise(Fluorescence = median(Fluorescence)),
                     aes(yintercept = Fluorescence)) +
          ggtitle(.y) +
          theme_minimal()
  ) %>%
  wrap_plots() %>%
  ggsave(filename = "Plots/control_data.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# unlike for the blanks, the control data are complete, but there are some outliers that
# are clearly contamination and need removing, e.g. Wells F07 and F08 in Plate X281021enzyme7
# plot without outliers
ggsave(
experimental %>%
  imap(~ .x$control %>%
         mutate(Group = .y %>% fct())) %>%
  bind_rows() %>%
  filter(!(Well %in% c("F07", "F08") & Group == "X281021enzyme7") &
         !(Well %in% c("G10", "G11", "G12") & Group == "X101221enzyme") &
         !(Well == "H08" & Group == "X091121enzyme4")) %>%
        ggplot(aes(Well, Fluorescence)) +
          geom_point(size = 3, shape = 16) +
          geom_hline(data = . %>% group_by(Group) %>% 
                       summarise(Fluorescence = median(Fluorescence)),
                     aes(yintercept = Fluorescence)) +
          facet_wrap(~ Group, scales = "free") +
          theme_minimal(),
filename = "Plots/control_data_nooutliers.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 4.3.1 Prior simulation ####
# Substrate autofluorescence is on average 500 arbitrary units relative to maximal 100 µM standard fluorescence.
# Since fluorescence cannot be negative, a gamma likelihood is best as for the water blank model. However, here
# there is no need for prediction of unobserved controls, because all plates included controls, so no 
# hierarchical model is necessary.

ggplot() +
  geom_density(aes(rnorm(1e5, log(200), 1) %>% exp())) +
  # scale_x_continuous(limits = c(0, 2e3), oob = scales::oob_keep) + # unhash to zoom in on peak
  theme_minimal()
# lots of variability

# here is the underlying distribution
ggplot() +
  geom_density(aes(rnorm(1e5, log(200), 1))) +
  theme_minimal()

# 4.3.2 Run model ####
experimental_control_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
}

parameters{
  real Fmu;
  real<lower=0> sigma;
}

model{
  Fmu ~ normal( log(200) , 2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( Fmu ); // Exponential transformation
  }

  // Gamma likelihood (generalsied linear model)
  Fluorescence ~ gamma( mu^2 / sigma^2 , mu / sigma^2 );
}
"

experimental_control_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
}

parameters{
  real Fmu;
  real<lower=0> sigma;
}

model{
  Fmu ~ normal( 500 , 100 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmu;
  }

  // Normal likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

experimental_control_mod <- experimental_control_stan %>%
  write_stan_file() %>%
  cmdstan_model()

experimental_control_samples <- experimental %>%
  map(~ experimental_control_mod$sample(data = .x$control %>%
                                          select(Fluorescence) %>%
                                          compose_data(),
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)
      )

# 4.3.3 Model checks ####
# check Rhat, effective sample size and chains
experimental_control_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

ggsave(
  experimental_control_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "Plots/control_rank.pdf", width = 80, height = 40, 
  unit = "cm", device = cairo_pdf)
# chains look good

# 4.3.4 Prior-posterior comparison ####
# sample prior
experimental_control_prior <- experimental %>%
  map(~ prior_samples(model = experimental_control_mod,
                      data = .x$control %>%
                        select(Fluorescence) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  experimental_control_prior %>%
    map2(experimental_control_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA),
                                 parameters = c("Fmu", "sigma"),
                                 format = "long")
         # %>% mutate(.value = if_else(.variable == "Fmu", exp(.value), .value))
         ) %>%
    imap(~ .x %>% prior_posterior_plot() + 
           # scale_x_continuous(limits = c(0, 200), oob = scales::oob_keep) + 
           ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "Plots/control_prior_posterior.pdf", width = 80, height = 40,
  unit = "cm", device = cairo_pdf)





####################















experimental_control_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

experimental_control_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

# 4.3.4 Prior-posterior comparison ####
# sample prior
experimental_control_prior <- prior_samples(
  model = experimental_control_mod,
  data = experimental %>%
    imap(~ .x$control %>%
           mutate(Group = .y %>% fct())) %>%
    bind_rows() %>%
    select(Fluorescence, Group) %>%
    compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
experimental_control_prior %>%
    prior_posterior_draws(posterior_samples = experimental_control_samples,
                          group = experimental %>%
                            imap(~ .x$control %>%
                                   mutate(Group = .y %>% fct())) %>%
                            bind_rows() %>%
                            select(Group),
                          parameters = c("Fmu[Group]", "sigma"),
                          format = "long") %>%
  prior_posterior_plot(group_name = "Group")
# Note the fuzziness of the prior. Stan doesn't seem to be able to smoothly estimate
# hierarchical priors.
# For smooth prior distributions, generate in R:

experimental_blank_prior <-
  tibble(.chain = 1:8 %>% rep(each = 1e4),
         .iteration = 1:1e4 %>% rep(8),
         .draw = 1:8e4,
         Fmu_mu = rnorm(8 * 1e4, log(200) , 1),
         Fmu_sigma = rexp(8 * 1e4, 1),
         Fmu = rnorm(8 * 1e4, Fmu_mu, Fmu_sigma),
         sigma = rexp(8 * 1e4, 1)) %>%
  pivot_longer(cols = -starts_with("."),
               names_to = ".variable", values_to = ".value") %>%
  mutate(rep = if_else(.variable == "Fmu", 24, 1),
         .variable = fct_relevel(.variable, "Fmu")) %>%
  uncount(rep) %>%
  arrange(.variable) %>%
  mutate(Group = experimental %>%
           imap(~ .x$blank %>%
                  mutate(Group = .y %>% fct())) %>%
           keep(~ nrow(.x) > 0) %>%
           bind_rows() %$% levels(Group) %>%
           rep(8 * 1e4) %>%
           c(NA %>% rep(3 * 8 * 1e4)))

experimental_blank_posterior <-
  experimental_blank_samples %>%
  recover_types(experimental %>%
                  imap(~ .x$blank %>%
                         mutate(Group = .y %>% fct())) %>%
                  keep(~ nrow(.x) > 0) %>%
                  bind_rows() %>% select(Group)) %>%
  gather_draws(Fmu[Group], Fmu_mu, Fmu_sigma, sigma) %>%
  ungroup()


experimental_blank_prior %>%
  mutate(distribution = "prior") %>%
  bind_rows(experimental_blank_posterior %>%
              mutate(distribution = "posterior")) %>%
  prior_posterior_plot(group_name = "Group")
# Not fuzzy. Everything fine.

# 4.2.5 Predictions ####
experimental_blank_predictions <-
  experimental_blank_samples %>%
  recover_types(experimental %>%
                  imap(~ .x$blank %>%
                         mutate(Group = .y %>% fct())) %>%
                  keep(~ nrow(.x) > 0) %>%
                  bind_rows() %>% select(Group)) %>%
  spread_draws(Fmu[Group], sigma) %>%
  mutate(Fmu = exp(Fmu),
         obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2)) %>%
  select(-sigma) %>%
  bind_rows(
    experimental_blank_samples %>%
      spread_draws(Fmu_mu, Fmu_sigma, sigma) %>%
      mutate(Fmu = exp( rnorm(n(), Fmu_mu, Fmu_sigma) ),
             obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2),
             Group = "unobserved") %>%
      select(-c(Fmu_mu, Fmu_sigma, sigma))
  ) %>%
  pivot_longer(cols = c("Fmu", "obs"),
               names_to = "level", values_to = "Fluorescence")

experimental_blank_predictions %>%
  ggplot(aes(Fluorescence, alpha = factor(level))) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0, 500), oob = scales::oob_keep) +
    scale_alpha_manual(values = c(0.6, 0.2)) +
    facet_wrap(~ Group, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# It makes most sense to me to use Fmu here rather than predicted observations,
# since we are replacing F0. That is while we want Fmu on the response scale, it is
# to substitute another parameter so does not need to pass through the likelihood.
experimental_blank_predictions %>%
  filter(level == "Fmu") %>%
  ggplot(aes(Fluorescence)) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0, 500), oob = scales::oob_keep) +
    facet_wrap(~ Group, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 4.4 Convert sample fluorescence ####








