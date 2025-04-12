# 1. Load data ####
# Load output from BMG CLARIOstar Plus microplate reader
require(tidyverse)
require(here)
files <- here("Microbes", "Raw") %>% list.files(pattern = "\\.csv$", full.names = TRUE)

enzyme <- files %>%
  map(~ read.csv(., skip = 9, header = TRUE) %>%
        rename(Fluorescence = "Raw.Data..360.20.450.30.") %>%
        mutate(Fluorescence = Fluorescence %>% as.numeric(),
               Well = Well %>% fct(),
               Content = Content %>% fct())) %>%
  set_names(str_remove(basename(files), "\\..*$") %>% make.names) %>%
  imap(~ .x %>% mutate(Date = str_extract(.y, "\\d+(?=enzyme)") %>% dmy(),
                       Plate2 = str_remove(.y, "^X") %>% fct()))

str(enzyme)

# Load ordered metadata
meta <- here("Microbes", "Meta.csv") %>% read.csv()

# Join annotations to respective measurements
require(magrittr)
enzyme %<>% 
  imap(~ bind_cols(.x, meta %>% filter(Plate == str_remove(.y, "^X"))) %>%
         mutate(Date2 = str_extract(Plate, "\\d+(?=enzyme)") %>% dmy(),
                Plate = Plate %>% fct()))
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
# or
enzyme %>% 
  map(~ .x %$% identical(Plate, Plate2))

# remove additional Date variable and change Plate to number
enzyme %<>% 
  map(~ .x %>% select(-c(Date2, Plate2)) %>%
        mutate(Plate_number = str_extract(Plate, "\\d+(?=[^\\d]*$)") %>%
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

rm(list = setdiff(ls(), c("enzyme", "meta", "beta", "F0", "Fmax")))

# 4. Experimental samples ####
experimental <- enzyme %>%
  map(~ .x %>% filter(Date >= "281021" %>% dmy() & 
                        !(Date == "031121" %>% dmy() & Plate_number == 9))) %>%
  keep(~ nrow(.x) > 0)

# Split each dataframe in experimental into its components
experimental %<>% map(~ list(
  blank = .x %>% filter(Annotation %in% c("AS", "MQ")),
  # control = .x %>% filter(Annotation %in% c("CAS", "CMQ")), # I chose to include controls in samples
  standard = .x %>%
    filter(str_detect(Annotation, "^[0-9.]+$")) %>% 
    rename(Concentration = Annotation) %>%
    mutate(Concentration = as.numeric(Concentration)) %>%
    filter(between(Concentration, 0, 100)),
  samples = .x %>% filter(!str_detect(Annotation, "^[0-9.]+$") & 
                            !Annotation %in% c("MQ", "AS"))
))

str(experimental)

experimental$X031121enzyme1$blank # examples
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
  ggsave(filename = "standard_data.pdf", path = here("Plots"),
         width = 80, height = 40, unit = "cm", device = cairo_pdf)

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
  filename = "standard_rank.pdf", path = here("Plots"),
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

experimental_standard_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots() %>%
  ggsave(filename = "standard_pairs.pdf", path = here("Plots"),
         width = 80, height = 80, unit = "cm", device = cairo_pdf)
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
  filename = "standard_prior_posterior.pdf", path = here("Plots"),
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# F0 broke out quite a bit, so the model predicts higher intercepts than expected

# 4.1.5 Predictions ####
experimental_standard_predictions <- experimental_standard_prior %>%
  map2(experimental_standard_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "short")) %>%
  map2(experimental, 
       ~ spread_continuous(.x, .y$standard, 
                           predictor_name = "Concentration",
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
  ggsave(filename = "standard_predictions.pdf", path = here("Plots"),
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# great fit except for those deviations in F0

# 4.2 Blank models ####
# To correctly calculate enzyme activities of experimental samples, F0 of the standard curve
# needs to be replaced with a water blank. Blanks were measured in triplicate or more replicates
# so in order to preserve the variability, intercept models must be built. The prediction for mu
# of these intercept models can then be substituted for F0.

# Visualise
experimental %>%
  imap(~ .x$blank %>%
        ggplot(aes(Well, Fluorescence)) +
          geom_point(size = 3, shape = 16) +
          geom_hline(data = . %>% summarise(Fluorescence = median(Fluorescence)),
                     aes(yintercept = Fluorescence)) +
          theme_minimal() +
          ggtitle(.y)
  ) %>%
  wrap_plots() %>% # too large to view in RStudio
  ggsave(filename = "blank_data.pdf", path = here("Plots"),
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# There are four plates missing blanks, which will have to be substituted with an average:
# X281021enzyme1, X281021enzyme3, X281021enzyme5, and X281021enzyme7.
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
  array[n] int Plate;
  int n_Plate;
}

parameters{
  real Fmu_mu; // Hyperparameters
  real<lower=0> Fmu_sigma;

  vector[n_Plate] Fmu;

  real<lower=0> sigma;
}

model{
  // Hyperpriors
  Fmu_mu ~ normal( log(200) , 1 );
  Fmu_sigma ~ exponential( 1 );

  // Plate-level prior
  Fmu ~ normal( Fmu_mu , Fmu_sigma );

  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( Fmu[Plate[i]] ); // Exponential transformation
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
    map(~ .x$blank) %>% # select only blank subset
    keep(~ nrow(.x) > 0) %>% # drop all list items (plates) without data
    bind_rows() %>% # here I am combining the blank dataframes to compare across plates
    select(Fluorescence, Plate) %>%
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
    map(~ .x$blank) %>%
    keep(~ nrow(.x) > 0) %>%
    bind_rows() %>%
    select(Fluorescence, Plate) %>%
    compose_data(),
  chains = 8, samples = 1e4)
# prior_samples() doesn't fully support multilevel models but
# for our visualisation purposes the messy samples are enough.

# plot prior-posterior comparison
experimental_blank_prior %>%
    prior_posterior_draws(posterior_samples = experimental_blank_samples,
                          group = experimental %>%
                            map(~ .x$blank) %>%
                            keep(~ nrow(.x) > 0) %>%
                            bind_rows() %>%
                            select(Plate),
                          parameters = c("Fmu[Plate]", "Fmu_mu", "Fmu_sigma", "sigma"),
                          format = "long") %>%
  prior_posterior_plot(group_name = "Plate")
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
  mutate(rep = if_else(.variable == "Fmu", 24, 1), # there are 24 plates
         .variable = fct_relevel(.variable, "Fmu")) %>%
  uncount(rep) %>%
  arrange(.variable) %>%
  mutate(Plate = experimental %>%
           map(~ .x$blank) %>%
           keep(~ nrow(.x) > 0) %>%
           bind_rows() %$% levels(Plate) %>%
           rep(8 * 1e4) %>%
           c(NA %>% rep(3 * 8 * 1e4)))

experimental_blank_posterior <-
  experimental_blank_samples %>%
  recover_types(experimental %>%
                  map(~ .x$blank) %>%
                  keep(~ nrow(.x) > 0) %>%
                  bind_rows() %>%
                  select(Plate)) %>%
  gather_draws(Fmu[Plate], Fmu_mu, Fmu_sigma, sigma) %>%
  ungroup()

experimental_blank_prior %>%
  mutate(distribution = "prior") %>%
  bind_rows(experimental_blank_posterior %>%
              mutate(distribution = "posterior")) %>%
  prior_posterior_plot(group_name = "Plate")
# Not fuzzy. Everything fine.

# 4.2.5 Predictions ####
experimental_blank_predictions <-
  experimental_blank_samples %>%
  recover_types(experimental %>%
                  map(~ .x$blank) %>%
                  keep(~ nrow(.x) > 0) %>%
                  bind_rows() %>%
                  select(Plate)) %>%
  spread_draws(Fmu[Plate], sigma) %>%
  mutate(Fmu = exp(Fmu),
         obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2)) %>%
  select(-sigma) %>%
  bind_rows(
    experimental_blank_samples %>%
      spread_draws(Fmu_mu, Fmu_sigma, sigma) %>%
      mutate(Fmu = exp( rnorm(n(), Fmu_mu, Fmu_sigma) ),
             obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2)) %>%
      select(-c(Fmu_mu, Fmu_sigma, sigma))
  ) %>%
  pivot_longer(cols = c("Fmu", "obs"),
               names_to = "level", values_to = "Fluorescence")

experimental_blank_predictions %>%
  ggplot(aes(Fluorescence, alpha = factor(level))) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0, 500), oob = scales::oob_keep) +
    scale_alpha_manual(values = c(0.6, 0.2)) +
    facet_wrap(~ Plate, scales = "free") +
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
    facet_wrap(~ Plate, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# # 4.3 Control models ####
# # To correctly calculate enzyme activities of experimental samples, a substrate autogenic
# # fluorescence and contamination control needs to be subtracted after conversion with the 
# # standard curve. If the substrate control were subtracted before conversion with the 
# # standard curve, water autofluorescence which is naturally contained within the substrate 
# # control would effectively be subtracted twice, once with the control and once during 
# # conversion with the standard curve which has the water blank as its intercept.
# 
# # Visualise
# experimental %>%
#   imap(~ .x$control %>%
#         ggplot(aes(Well, Fluorescence)) +
#           geom_point(size = 3, shape = 16) +
#           geom_hline(data = . %>% summarise(Fluorescence = median(Fluorescence)),
#                      aes(yintercept = Fluorescence)) +
#           theme_minimal() +
#           ggtitle(.y)
#   ) %>%
#   wrap_plots() %>%
#   ggsave(filename = "control_data.pdf", path = here("Plots"),
#          width = 80, height = 40, unit = "cm", device = cairo_pdf)
# # unlike for the blanks, the control data are complete
# 
# # 4.3.1 Prior simulation ####
# # Substrate autofluorescence is on average 500 arbitrary units relative to maximal 100 µM standard fluorescence.
# # Since fluorescence cannot be negative, a gamma likelihood is best as for the water blank model. However, here
# # there is no need for prediction of unobserved controls, because all plates included controls, so no 
# # hierarchical model is necessary.
# 
# ggplot() +
#   geom_density(aes(rnorm(1e5, log(500), 1) %>% exp())) +
#   # scale_x_continuous(limits = c(0, 2e3), oob = scales::oob_keep) + # unhash to zoom in on peak
#   theme_minimal()
# # lots of variability
# 
# # here is the underlying distribution
# ggplot() +
#   geom_density(aes(rnorm(1e5, log(500), 1))) +
#   theme_minimal()
# 
# # Additionally, via reparameterisation of the gamma likelihood in terms of mean (mu) and scale (theta),
# # gamma( mu / theta, 1 / theta), I discovered that exponential(1) is too restrictive a prior on the 
# # likelihood uncertainty. It worked with mean and scale, but not with mean and sd (sigma). 
# # I kept my parameterisation but changed the prior to exponential(0.01).
# 
# # 4.3.2 Run model ####
# experimental_control_stan <- "
# data{
#   int n;
#   vector<lower=0>[n] Fluorescence;
#   array[n] int Plate;
#   int n_Plate;
# }
# 
# parameters{
#   real Fmu_mu; // Hyperparameters
#   real<lower=0> Fmu_sigma;
# 
#   vector[n_Plate] Fmu;
# 
#   real<lower=0> sigma;
# }
# 
# model{
#   // Hyperpriors
#   Fmu_mu ~ normal( log(500) , 1 );
#   Fmu_sigma ~ exponential( 1 );
# 
#   // Plate-level prior
#   Fmu ~ normal( Fmu_mu , Fmu_sigma );
# 
#   sigma ~ exponential( 0.01 );
# 
#   // Model
#   vector[n] mu;
#   for ( i in 1:n ) {
#     mu[i] = exp( Fmu[Plate[i]] ); // Exponential transformation
#   }
# 
#   // Gamma likelihood (generalsied linear model)
#   Fluorescence ~ gamma( mu^2 / sigma^2 , mu / sigma^2 );
# }
# "
# 
# experimental_control_mod <- experimental_control_stan %>%
#   write_stan_file() %>%
#   cmdstan_model()
# 
# experimental_control_samples <- experimental_control_mod$sample(
#   data = experimental %>%
#     map(~ .x$control) %>% # select only control subset
#     bind_rows() %>%
#     select(Fluorescence, Plate) %>%
#     compose_data(),
#   chains = 8,
#   parallel_chains = parallel::detectCores(),
#   iter_warmup = 1e4,
#   iter_sampling = 1e4)
# 
# # 4.3.3 Model checks ####
# # check Rhat, effective sample size and chains
# experimental_control_samples$summary() %>%
#   mutate(rhat_check = rhat > 1.001) %>%
#   summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
#             rhat_mean = mean(rhat),
#             rhat_sd = sd(rhat),
#             ess_mean = mean(ess_bulk),
#             ess_sd = sd(ess_bulk))
# # no rhat above 1.001
# # good effective sample size
# 
# experimental_control_samples$draws(format = "df") %>%
#   mcmc_rank_overlay()
# # chains look good
# 
# # 4.3.4 Prior-posterior comparison ####
# # sample prior
# experimental_control_prior <- prior_samples(
#   model = experimental_control_mod,
#   data = experimental %>%
#     map(~ .x$control) %>%
#     bind_rows() %>%
#     select(Fluorescence, Plate) %>%
#     compose_data(),
#   chains = 8, samples = 1e4)
# 
# # plot prior-posterior comparison
# experimental_control_prior %>%
#     prior_posterior_draws(posterior_samples = experimental_control_samples,
#                           group = experimental %>%
#                             map(~ .x$control) %>%
#                             bind_rows() %>%
#                             select(Plate),
#                           parameters = c("Fmu[Plate]", "Fmu_mu", "Fmu_sigma", "sigma"),
#                           format = "long") %>%
#   prior_posterior_plot(group_name = "Plate")
# 
# # prior-posterior comparison in R
# experimental_control_prior <-
#   tibble(.chain = 1:8 %>% rep(each = 1e4),
#          .iteration = 1:1e4 %>% rep(8),
#          .draw = 1:8e4,
#          Fmu_mu = rnorm(8 * 1e4, log(500) , 1),
#          Fmu_sigma = rexp(8 * 1e4, 1),
#          Fmu = rnorm(8 * 1e4, Fmu_mu, Fmu_sigma),
#          sigma = rexp(8 * 1e4, 0.1)) %>%
#   pivot_longer(cols = -starts_with("."),
#                names_to = ".variable", values_to = ".value") %>%
#   mutate(rep = if_else(.variable == "Fmu", 28, 1), # there are 28 plates
#          .variable = fct_relevel(.variable, "Fmu")) %>%
#   uncount(rep) %>%
#   arrange(.variable) %>%
#   mutate(Plate = experimental %>%
#            map(~ .x$control) %>%
#            bind_rows()  %$% levels(Plate) %>%
#            rep(8 * 1e4) %>%
#            c(NA %>% rep(3 * 8 * 1e4)))
# 
# experimental_control_posterior <-
#   experimental_control_samples %>%
#   recover_types(experimental %>%
#                   map(~ .x$control) %>%
#                   bind_rows() %>%
#                   select(Plate)) %>%
#   gather_draws(Fmu[Plate], Fmu_mu, Fmu_sigma, sigma) %>%
#   ungroup()
# 
# experimental_control_prior %>%
#   mutate(distribution = "prior") %>%
#   bind_rows(experimental_control_posterior %>%
#               mutate(distribution = "posterior")) %>%
#   prior_posterior_plot(group_name = "Plate")
# 
# # 4.3.5 Predictions ####
# experimental_control_predictions <-
#   experimental_control_samples %>%
#   recover_types(experimental %>%
#                   map(~ .x$control) %>%
#                   bind_rows() %>%
#                   select(Plate)) %>%
#   spread_draws(Fmu[Plate], sigma) %>%
#   mutate(Fmu = exp(Fmu),
#          obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2)) %>%
#   select(-sigma) %>%
#   bind_rows(
#     experimental_control_samples %>%
#       spread_draws(Fmu_mu, Fmu_sigma, sigma) %>%
#       mutate(Fmu = exp( rnorm(n(), Fmu_mu, Fmu_sigma) ),
#              obs = rgamma(n(), Fmu^2 / sigma^2, Fmu / sigma^2)) %>%
#       select(-c(Fmu_mu, Fmu_sigma, sigma))
#   ) %>%
#   pivot_longer(cols = c("Fmu", "obs"),
#                names_to = "level", values_to = "Fluorescence")
# 
# experimental_control_predictions %>%
#   ggplot(aes(Fluorescence, alpha = factor(level))) +
#     geom_density(colour = NA, fill = "black") +
#     scale_x_continuous(limits = c(0, 2e3), oob = scales::oob_keep) +
#     scale_alpha_manual(values = c(0.6, 0.2)) +
#     facet_wrap(~ Plate, scales = "free") +
#     theme_minimal() +
#     theme(panel.grid = element_blank())
# # It makes most sense to me to use Fmu since that is what I used for the water blank.
# experimental_control_predictions %>%
#   filter(level == "Fmu") %>%
#   ggplot(aes(Fluorescence)) +
#     geom_density(colour = NA, fill = "black") +
#     scale_x_continuous(limits = c(0, 1.2e3), oob = scales::oob_keep) +
#     facet_wrap(~ Plate, scales = "free") +
#     theme_minimal() +
#     theme(panel.grid = element_blank())
# 
# rm(list = setdiff(ls(), c("enzyme", "meta", "experimental", "experimental_standard_samples", 
#                           "experimental_blank_samples", "experimental_standard_samples", 
#                           "experimental_blank_predictions", "experimental_control_predictions")))

# 4.4 Conversion ####
# solving F = Fmax * beta * c / ( Fmax + beta * c ) + F0 for c
# gives c = Fmax * ( F - F0 ) / ( beta * ( F0 + Fmax - F ) )

# 4.4.1 Combine posteriors ####
experimental_standard <- experimental_standard_samples %>%
  imap(~ .x %>% spread_draws(Fmax, beta, F0) %>%
         mutate(Plate = str_remove(.y, "^X") %>% fct())) %>%
  bind_rows()

experimental_blank <- experimental_blank_predictions %>%
  filter(!is.na(Plate)) %>% # filter out predictions for unobserved plates
  bind_rows(
    experimental_blank_predictions %>%
      filter(is.na(Plate)) %>% # replicate prediction for unobserved plates
      slice(rep(1:n(), each = 4)) %>% # and allocate unobserved plate names
      mutate(Plate = c("281021enzyme1", "281021enzyme3",
                       "281021enzyme5", "281021enzyme7") %>%
               rep(length.out = n()) %>% fct())
  ) %>%
  filter(level == "Fmu") %>%
  select(-level) %>%
  rename(blank = Fluorescence)
  
# # Here's how to do the opposite: turn a grouped tibble into a named list
# experimental_blank_predictions %>%
#   filter(level == "Fmu") %>%
#   select(-level) %>%
#   rename(blank = Fluorescence) %>%
#   group_by(Group) %>%
#   group_split() %>%
#   set_names(
#     map(., ~ .x$Group[1])
#   ) %>%
#   map(~ .x %>% select(-Group))
# # But this enables less control over calculations that require group matching

# experimental_control <- experimental_control_predictions %>%
#   filter(!is.na(Plate) & level == "Fmu") %>% # filter out predictions for unobserved plates
#   select(-level) %>% # (no predictions for unobserved plates to be made here)
#   rename(control = Fluorescence)

experimental_parameters <- experimental_standard %>%
  full_join(experimental_blank,
            by = c("Plate", ".chain", ".iteration", ".draw"))


# 4.4.2 Convert sample fluorescence ####
experimental_samples <- experimental %>%
  map(~ .x$samples) %>% # select samples
  bind_rows() %>% # convert list to tibble
  rename(ID = Annotation) %>%
  mutate(ID = if_else(ID %in% c("CMQ", "CAS"), # create plate-specific control name
                      str_c(Plate, "control", sep = "_"), 
                      ID) %>% fct()) %>%
  full_join(experimental_parameters,
            by = "Plate",
            relationship = "many-to-many") %>%
  # mutate(# replace F0 with blank and convert sample and control fluorescence to enzyme activity
  #        Activity_raw = ( Fmax * ( Fluorescence - blank ) / ( beta * ( blank + Fmax - Fluorescence ) ) ),
  #        Activity_control = ( Fmax * ( control - blank ) / ( beta * ( blank + Fmax - control ) ) ),
  #        # correct sample enzyme activity by subtracting the converted control
  #        Activity = Activity_raw - Activity_control, # µM h^-1
  #        ID = ID %>% fct())
  mutate(# replace F0 with blank and convert sample and control fluorescence to enzyme activity
         Activity = ( Fmax * ( Fluorescence - blank ) / ( beta * ( blank + Fmax - Fluorescence ) ) ))
# maybe add sigma?
str(experimental_samples)

# 4.5 Technical replication ####
# In addition to measurement error introduced by the standard curve, control, and standardisation,
# wells vary across microplates due to composition, contamination or pipetting error. This source
# of error was accounted for by including five replicate wells for each biological replicate. The
# most intuitive way to add this source of error to the mix is by estimating an intercept across
# technical replicates for each biological replicate.

# 4.5.1 Summarise data ####
experimental_samples_summary <- experimental_samples %>%
  mutate(Replicate = fct_cross(ID, Well, sep = "_")) %>%
  group_by(Replicate, Well, Content, Fluorescence, Date,
           Plate, ID, Plate_number) %>%
  summarise(Activity_mean = mean(Activity), # this
            Activity_sd = sd(Activity),
            n = length(Activity)) %>%
  ungroup()


################################################

technical_stan <- "
data{
  int n;
  vector<lower=0>[n] Activity;
  array[n] int ID;
  int n_ID;
}

parameters{
  real Fmu_mu; // Hyperparameters
  real<lower=0> Fmu_sigma;

  vector[n_Plate] Fmu;

  real<lower=0> sigma;
}

model{
  // Hyperpriors
  Fmu_mu ~ normal( log(500) , 1 );
  Fmu_sigma ~ exponential( 1 );

  // Plate-level prior
  Fmu ~ normal( Fmu_mu , Fmu_sigma );

  sigma ~ exponential( 0.01 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( Fmu[Plate[i]] ); // Exponential transformation
  }

  // Gamma likelihood (generalsied linear model)
  Fluorescence ~ gamma( mu^2 / sigma^2 , mu / sigma^2 );
}
"

experimental_control_mod <- experimental_control_stan %>%
  write_stan_file() %>%
  cmdstan_model()

experimental_control_samples <- experimental_control_mod$sample(
  data = experimental %>%
    map(~ .x$control) %>% # select only control subset
    bind_rows() %>%
    select(Fluorescence, Plate) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4)







# 4.6 Standardisation ####
# 4.6.1 Load area and mass data ####
area <- here("Microbes", "Disc_Area.csv") %>% read.csv() %>%
  mutate(Species = Species %>% fct())
# Disc area (cm^2) was measured for a few random replicates of each species.
area %>% 
  ggplot(aes(Species, Area)) +
    geom_hline(yintercept = area %$% mean(Area)) +
    geom_point(shape = 16, alpha = 0.2, size = 3) +
    geom_point(data = . %>% group_by(Species) %>% 
                 summarise(Area = mean(Area)),
               aes(Species, Area), size = 5, shape = 21, fill = NA) +
  theme_minimal() +
  theme(panel.grid = element_blank())

mass <- here("Microbes", "Disc_Mass.csv") %>% read.csv() %>%
  mutate(Species = Species %>% fct(),
         ID = ID %>% fct(),
         Date = Date %>% dmy())
# Disc mass (g) was measured for each disc that underwent the enzyme assay.
mass %>% 
  filter(Species != "Sediment") %>% # filter out heavy sediment for clarity
  ggplot(aes(Species, Mass)) +
    geom_hline(yintercept = mass %>% filter(Species != "Sediment") %$% mean(Mass)) +
    geom_point(shape = 16, alpha = 0.2, size = 3) +
    geom_point(data = . %>% 
                 filter(Species != "Sediment") %>%
                 group_by(Species) %>% 
                 summarise(Mass = mean(Mass)),
               aes(Species, Mass), size = 5, shape = 21, fill = NA) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 4.6.2 Match mass ####
# Extracellular enzyme activity can either be reported per cm^2 or per g. I will have a look at both.
# Because disc area was not measured for each disc, it needs to be approximated by modelling. Disc mass
# on the other hand can directly be included in the calculation.

experimental_samples %<>%
  left_join(mass %>%
              rename(Date_weighed = Date), 
            by = "ID")

# 4.6.3 Model area ####
# The grand mean seems to be close to 0.55 cm^2.
ggplot() +
  geom_density(aes(rnorm(1e5, log(0.55), 0.2) %>% exp())) +
  theme_minimal()
# lots of variability

# here is the underlying distribution
ggplot() +
  geom_density(aes(rnorm(1e5, log(0.55), 0.2))) +
  theme_minimal()

# Run model
area_stan <- "
data{
  int n;
  vector<lower=0>[n] Area;
  array[n] int Species;
  int n_Species;
}

parameters{
  real Amu_mu; // Hyperparameters
  real<lower=0> Amu_sigma;

  vector[n_Species] Amu;

  real<lower=0> sigma;
}

model{
  // Hyperpriors
  Amu_mu ~ normal( log(0.55) , 0.2 );
  Amu_sigma ~ exponential( 1 );

  // Plate-level prior
  Amu ~ normal( Amu_mu , Amu_sigma );

  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( Amu[Species[i]] ); // Exponential transformation
  }

  // Gamma likelihood (generalsied linear model)
  Area ~ gamma( mu^2 / sigma^2 , mu / sigma^2 );
}
"

area_mod <- area_stan %>%
  write_stan_file() %>%
  cmdstan_model()

area_samples <- area_mod$sample(
  data = area %>%
    select(Area, Species) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4)

# check Rhat, effective sample size and chains
area_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

area_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

# sample prior
source("functions.R")
area_prior <- prior_samples(
  model = area_mod,
  data = area %>%
    select(Area, Species) %>%
    compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
area_prior %>%
    prior_posterior_draws(posterior_samples = area_samples,
                          group = area %>% select(Species),
                          parameters = c("Amu[Species]", "Amu_mu", "Amu_sigma", "sigma"),
                          format = "long") %>%
  prior_posterior_plot(group_name = "Species")
# priors are fine

# make predictions
area_predictions <-
  area_samples %>%
  recover_types(area %>% select(Species)) %>%
  spread_draws(Amu[Species], sigma) %>%
  mutate(Amu = exp(Amu),
         obs = rgamma(n(), Amu^2 / sigma^2, Amu / sigma^2)) %>%
  select(-sigma) %>%
  bind_rows(
    area_samples %>%
      spread_draws(Amu_mu, Amu_sigma, sigma) %>%
      mutate(Amu = exp( rnorm(n(), Amu_mu, Amu_sigma) ),
             obs = rgamma(n(), Amu^2 / sigma^2, Amu / sigma^2)) %>%
      select(-c(Amu_mu, Amu_sigma, sigma))
  ) %>%
  pivot_longer(cols = c("Amu", "obs"),
               names_to = "level", values_to = "Area")

area_predictions %>%
  ggplot(aes(Area, alpha = factor(level))) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0.3, 0.8), oob = scales::oob_keep) +
    scale_alpha_manual(values = c(0.6, 0.2)) +
    facet_wrap(~ Species, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# I will select only obs because I am interested including the potnetial variation
# caused by cutting a new disc.
area_predictions %>%
  filter(level == "obs") %>%
  ggplot(aes(Area)) +
    geom_density(colour = NA, fill = "black") +
    scale_x_continuous(limits = c(0.3, 0.8), oob = scales::oob_keep) +
    scale_alpha_manual(values = c(0.6, 0.2)) +
    facet_wrap(~ Species, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# There is only one macroalgal species for which I have data but didn't measure
# disc area: Fucus serratus. Hence the area prediction for unobserved species can
# simply be assigned to F. serratus.

area_predictions %<>%
  filter(level == "obs") %>%
  select(-level) %>%
  mutate(Species = if_else(is.na(Species), "Fucus serratus", Species) %>%
           fct())
  
# 4.6.4 Match area ####
experimental_samples %<>%
  left_join(area_predictions, 
            by = c("Species", ".chain", ".iteration", ".draw"))

# 4.6.5 Standardise ####
# Activity needs to be expressed per volume for plankton, per mass for sediment
# and macroalgae and per area for macroalgae. All incubations lasted one hour and
# were done in 1 mL (0.001 L) of sterile seawater. The molar mass of the enzyme 
# substrate is 338.31 g mol^-1. Since activity is presently given as µM h^-1, or
# µmol L^-1 h^-1, it needs to be multiplied by 0.001 L and 338.31 g mol^-1 first.

experimental_samples %<>%
  mutate(Activity_mL = Activity * 0.001 * 338.31, # µg mL^-1 h^-1
         Activity_g = Activity * 0.001 * 338.31 / Mass, # µg g^-1 h^-1
         Activity_cm2 = Activity * 0.001 * 338.31 / Area) # µg cm^-2 h^-1

rm(list = setdiff(ls(), c("experimental_samples")))

# 5. Grazing experiment ####
# 5.1 Prepare data ####
grazing <- experimental_samples %>%
  filter(str_starts(ID, "B") | str_starts(ID, "C")) %>%
  mutate(Well = Well %>% fct_drop(),
         Content = Content %>% fct_drop(),
         ID = ID %>% fct_drop(),
         Species = Species %>% fct_drop(),
         Replicate = fct_cross(ID, Well, sep = "_"))

snails <- here("Snails", "Snails.csv") %>% read.csv() %>%
  mutate(ID = ID %>% fct(), Experiment = Experiment %>% fct(),
         Flask = Flask %>% fct(), Start = Start %>% dmy_hm(),
         End = End %>% dmy_hm(), Species = Species %>% fct(), 
         Treatment = Treatment %>% fct(),
         Days_accurate = Start %--% End %>% as.duration() / ddays())

grazing %<>%
  left_join(snails, by = c("ID", "Species")) %>%
  mutate(ID = ID %>% fct_drop())

grazing_summary <- grazing %>%
  group_by(Experiment, ID, Replicate, Species, Treatment) %>%
  summarise(Activity_g_mean = mean(Activity_g),
            Activity_g_sd = sd(Activity_g),
            Activity_cm2_mean = mean(Activity_cm2),
            Activity_cm2_sd = sd(Activity_cm2),
            n = length(Activity_g)) %>%
  ungroup()

# 5.2 Visualise data ####
grazing_summary %>%
  ggplot() +
    stat_slab(aes(x = Activity_g_mean, y = Species, alpha = Treatment),
              height = 1, colour = "black", size = 0.5, justification = - 0.3) +
    # alternative to justification is Species %>% as.numeric() - 0.3 for geom_point
    # geom_point(aes(x = Activity_g_mean, y = Species, alpha = Treatment),
    #            position = position_jitter(height = 0.2), shape = 16) +
    geom_pointrange(data = grazing_summary %>%
                      group_by(ID, Species, Treatment) %>%
                      summarise(Activity_g_mean_mean = mean(Activity_g_mean),
                                Activity_g_mean_sd = sd(Activity_g_mean)),
                    aes(x = Activity_g_mean_mean, xmin = Activity_g_mean_mean - Activity_g_mean_sd,
                        xmax = Activity_g_mean_mean + Activity_g_mean_sd, y = Species, alpha = Treatment),
                    position = position_jitter(height = 0.2), shape = 16) +
    scale_alpha_manual(values = c(0.6, 0.4, 0.2)) +
    theme_minimal() +
    theme(panel.grid = element_blank())

grazing_summary %>%
  ggplot() +
  stat_slab(aes(x = Activity_cm2_mean, y = Species, alpha = Treatment),
            height = 1, colour = "black", size = 0.5, justification = - 0.3) +
  # alternative to justification is Species %>% as.numeric() - 0.3 for geom_point
  # geom_point(aes(x = Activity_g_mean, y = Species, alpha = Treatment),
  #            position = position_jitter(height = 0.2), shape = 16) +
  geom_pointrange(data = grazing_summary %>%
                    group_by(ID, Species, Treatment) %>%
                    summarise(Activity_cm2_mean_mean = mean(Activity_cm2_mean),
                              Activity_cm2_mean_sd = sd(Activity_cm2_mean)),
                  aes(x = Activity_cm2_mean_mean, xmin = Activity_cm2_mean_mean - Activity_cm2_mean_sd,
                      xmax = Activity_cm2_mean_mean + Activity_cm2_mean_sd, y = Species, alpha = Treatment),
                  position = position_jitter(height = 0.2), shape = 16) +
  scale_alpha_manual(values = c(0.6, 0.4, 0.2)) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# 5.3 Model activity per gram ####
# 5.3.1 Prior simulation ####
grazing_summary %$%
  mean(Activity_g_mean)

ggplot() +
  geom_density(aes(rnorm(1e5, log(75), 0.5) %>% exp())) +
  theme_minimal()

# here is the underlying distribution
ggplot() +
  geom_density(aes(rnorm(1e5, log(75), 0.5))) +
  theme_minimal()

# 5.3.2 Run model ####
grazing_stan <- "
data{
  int n;
  vector<lower=0>[n] Activity_g_mean;
  vector<lower=0>[n] Activity_g_sd;
  array[n] int Species;
  int n_Species;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int ID;
  int n_ID;
}

parameters{
  vector<lower=0>[n] Activity_g;
  real alpha;
  vector[n_Species] beta_S;
  vector[n_Treatment] beta_T;
  matrix[n_Species, n_Treatment] beta_S_T;
  vector[n_ID] beta_ID;
  real<lower=0> sigma;
}

model{
  alpha ~ normal( log(75) , 0.5 );
  beta_S ~ normal( 0 , 0.5 );
  beta_T ~ normal( 0 , 0.5 );
  to_vector(beta_S_T) ~ normal( 0 , 0.5 );
  beta_ID ~ normal( 0 , 0.5 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( alpha + beta_S[Species[i]] + beta_T[Treatment[i]] + 
                 beta_S_T[Species[i], Treatment[i]] + beta_ID[ID[i]] );
  }

  // Gamma likelihood (generalised linear model)
  Activity_g ~ gamma( mu^2 / sigma^2 , mu / sigma^2 );
  // with normal measurement error
  Activity_g_mean ~ normal( Activity_g , Activity_g_sd );
}
"

grazing_mod <- grazing_stan %>%
  write_stan_file() %>%
  cmdstan_model()

grazing_samples <- grazing_mod$sample(
  data = grazing_summary %>%
    select(Activity_g_mean, Activity_g_sd,
           Species, Treatment, ID) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4)

# 5.3.3 Model checks ####
grazing_samples %>%
  recover_types(grazing_summary %>%
                  select(Species, Treatment, ID)) %>%
  spread_draws(alpha, beta_S[Species], beta_T[Treatment], beta_S_T[Species, Treatment], beta_ID[ID]) %>%
  ungroup() %>%
  mutate(sum = alpha + beta_S + beta_T + beta_S_T + beta_ID) %>%
  group_by(Species, Treatment) %>%
  summarise(mean = mean(sum),
            sd = sd(sum))

# with ID
# Species                Treatment               mean    sd
# <fct>                  <fct>                  <dbl> <dbl>
# 1 Laminaria digitata     Steromphala cineraria   3.58  1.46
# 2 Laminaria digitata     Autogenic control       2.85  1.43
# 3 Laminaria digitata     Calliostoma zizyphinum  2.76  1.47
# 4 Laminaria hyperborea   Steromphala cineraria   3.67  1.46
# 5 Laminaria hyperborea   Autogenic control       3.00  1.45
# 6 Laminaria hyperborea   Calliostoma zizyphinum  3.33  1.46
# 7 Saccorhiza polyschides Steromphala cineraria   3.43  1.53
# 8 Saccorhiza polyschides Autogenic control       2.88  1.45
# 9 Saccorhiza polyschides Calliostoma zizyphinum  3.58  1.48
# 10 Saccharina latissima   Steromphala cineraria   3.78  1.47
# 11 Saccharina latissima   Autogenic control       3.32  1.46
# 12 Saccharina latissima   Calliostoma zizyphinum  3.54  1.47

# without ID
# Species                Treatment               mean     sd
# <fct>                  <fct>                  <dbl>  <dbl>
# 1 Laminaria digitata     Steromphala cineraria   4.37 0.107 
# 2 Laminaria digitata     Autogenic control       3.85 0.0890
# 3 Laminaria digitata     Calliostoma zizyphinum  3.83 0.132 
# 4 Laminaria hyperborea   Steromphala cineraria   4.25 0.115 
# 5 Laminaria hyperborea   Autogenic control       3.91 0.0895
# 6 Laminaria hyperborea   Calliostoma zizyphinum  4.03 0.117 
# 7 Saccorhiza polyschides Steromphala cineraria   4.41 0.126 
# 8 Saccorhiza polyschides Autogenic control       4.37 0.114 
# 9 Saccorhiza polyschides Calliostoma zizyphinum  4.42 0.114 
# 10 Saccharina latissima   Steromphala cineraria   4.19 0.120 
# 11 Saccharina latissima   Autogenic control       4.40 0.0899
# 12 Saccharina latissima   Calliostoma zizyphinum  4.23 0.114 

# with ID but not multilevel
# Species                Treatment               mean    sd
# <fct>                  <fct>                  <dbl> <dbl>
# 1 Laminaria digitata     Steromphala cineraria   3.90 0.998
# 2 Laminaria digitata     Autogenic control       2.95 0.987
# 3 Laminaria digitata     Calliostoma zizyphinum  2.51 1.01 
# 4 Laminaria hyperborea   Steromphala cineraria   3.89 0.999
# 5 Laminaria hyperborea   Autogenic control       3.10 1.01 
# 6 Laminaria hyperborea   Calliostoma zizyphinum  3.41 1.00 
# 7 Saccorhiza polyschides Steromphala cineraria   4.05 1.10 
# 8 Saccorhiza polyschides Autogenic control       3.62 1.04 
# 9 Saccorhiza polyschides Calliostoma zizyphinum  4.16 0.998
# 10 Saccharina latissima   Steromphala cineraria   3.96 0.998
# 11 Saccharina latissima   Autogenic control       3.74 1.03 
# 12 Saccharina latissima   Calliostoma zizyphinum  3.69 1.03 


grazing_summary %>%
  group_by(Species, Treatment) %>%
  summarise(mean = mean(Activity_g_mean),
            sd = sd(Activity_g_mean))

# 5.3.4 Prior-posterior comparison ####


# 5.3.5 Predictions ####


# 5.4 Model activity per cm^2 ####


# 6. Decomposition experiment ####


# 7. Remaining samples ####


