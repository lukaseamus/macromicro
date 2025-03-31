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
kinetics_standard_ht_mod <- cmdstan_model(stan_file = write_stan_file(code = kinetics_standard_ht_stan))

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
kinetics_standard_ht_summary <- kinetics_standard_ht_samples %>%
  map(~ .x$summary())

kinetics_standard_ht_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

kinetics_standard_ht_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

kinetics_standard_ht_draws <- kinetics_standard_ht_samples %>%
  map(~ .x$draws(format = "df"))

require(bayesplot)
kinetics_standard_ht_draws %>%
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_ht_draws %>%
  map(~ mcmc_pairs(.x, pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# practically no correlation

# 2.1.4 Prior-posterior comparison ####
kinetics_standard_ht_prior_posterior <- kinetics_standard_ht_samples %>%
  map(~ spread_draws(.x, Fmax, beta, F0, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1),
          F0_prior = rgamma(length(.draw), 500^2 / 300^2, 500 / 300^2),
          beta_prior = rgamma(length(.draw), 5e3^2 / 3e3^2, 5e3 / 3e3^2),
          Fmax_prior = rgamma(length(.draw), 24e4^2 / 5e4^2, 24e4 / 5e4^2)
        )
  )

str(kinetics_standard_ht_prior_posterior)

kinetics_standard_ht_prior_posterior %>%
  map(~ pivot_longer(., cols = c("F0", "F0_prior", "beta", "beta_prior", "Fmax", "Fmax_prior"),
                     names_to = c("Parameter", "Distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message
                     values_to = "Samples") %>%
        mutate(Parameter = fct(Parameter),
               Distribution = fct(ifelse(is.na(Distribution), # here the NAs are dealt with
                                         "posterior", Distribution))) %>%
        ggplot(aes(Samples, fill = Distribution)) +
        facet_wrap(~ Parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# priors are not restrictive

# 2.1.5 Predictions for mu ####
kinetics_standard_ht_mu <- kinetics_standard_ht_prior_posterior %>%
  map2(kinetics,
       ~ .x %>% select(Fmax, beta, F0) %>%
         expand_grid(Concentration = .y$standard %$%
                       seq(min(Concentration), max(Concentration), length.out = 50)) %>%
         mutate(mu = Fmax * tanh( beta * Concentration / Fmax ) + F0) %>%
         select(mu, Concentration)
  )

kinetics_standard_ht_mu_summary <- kinetics_standard_ht_mu %>%
  map(~ .x %>% group_by(Concentration) %>%
        mean_qi(mu, .width = c(.5, .8, .9)))

kinetics %>%
  map2(kinetics_standard_ht_mu_summary,
       ~ ggplot() +
            geom_point(data = .x$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .y, aes(Concentration, mu)) +
            geom_ribbon(data = .y, aes(Concentration, ymin = .lower, ymax = .upper,
                                    alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal()
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

kinetics_standard_es_mod <- cmdstan_model(stan_file = write_stan_file(code = kinetics_standard_es_stan))

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
kinetics_standard_es_summary <- kinetics_standard_es_samples %>%
  map(~ .x$summary())

kinetics_standard_es_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

kinetics_standard_es_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

kinetics_standard_es_draws <- kinetics_standard_es_samples %>%
  map(~ .x$draws(format = "df"))

kinetics_standard_es_draws %>%
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_es_draws %>%
  map(~ mcmc_pairs(.x, pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# practically no correlation

# 2.1.9 Prior-posterior comparison ####
kinetics_standard_es_prior_posterior <- kinetics_standard_es_samples %>%
  map(~ spread_draws(.x, Fmax, beta, F0, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1),
          F0_prior = rgamma(length(.draw), 500^2 / 300^2, 500 / 300^2),
          beta_prior = rgamma(length(.draw), 5e3^2 / 3e3^2, 5e3 / 3e3^2),
          Fmax_prior = rgamma(length(.draw), 24e4^2 / 5e4^2, 24e4 / 5e4^2)
        )
  )

str(kinetics_standard_es_prior_posterior)

kinetics_standard_es_prior_posterior %>%
  map(~ pivot_longer(., cols = c("F0", "F0_prior", "beta", "beta_prior", "Fmax", "Fmax_prior"),
                     names_to = c("Parameter", "Distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message
                     values_to = "Samples") %>%
        mutate(Parameter = fct(Parameter),
               Distribution = fct(ifelse(is.na(Distribution), # here the NAs are dealt with
                                         "posterior", Distribution))) %>%
        ggplot(aes(Samples, fill = Distribution)) +
        facet_wrap(~ Parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# priors are not restrictive

# 2.1.10 Predictions for mu ####
kinetics_standard_es_mu <- kinetics_standard_es_prior_posterior %>%
  map2(kinetics,
       ~ .x %>% select(Fmax, beta, F0) %>%
         expand_grid(Concentration = .y$standard %$%
                       seq(min(Concentration), max(Concentration), length.out = 50)) %>%
         mutate(mu = Fmax * ( 1 - exp( -beta * Concentration / Fmax ) ) + F0) %>%
         select(mu, Concentration)
  )

kinetics_standard_es_mu_summary <- kinetics_standard_es_mu %>%
  map(~ .x %>% group_by(Concentration) %>%
        mean_qi(mu, .width = c(.5, .8, .9)))

kinetics %>%
  map2(kinetics_standard_es_mu_summary,
       ~ ggplot() +
            geom_point(data = .x$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .y, aes(Concentration, mu)) +
            geom_ribbon(data = .y, aes(Concentration, ymin = .lower, ymax = .upper,
                                    alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal()
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

kinetics_standard_rh_mod <- cmdstan_model(stan_file = write_stan_file(code = kinetics_standard_rh_stan))

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
kinetics_standard_rh_summary <- kinetics_standard_rh_samples %>%
  map(~ .x$summary())

kinetics_standard_rh_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

kinetics_standard_rh_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

kinetics_standard_rh_draws <- kinetics_standard_rh_samples %>%
  map(~ .x$draws(format = "df"))

kinetics_standard_rh_draws %>%
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# chains look good

kinetics_standard_rh_draws %>%
  map(~ mcmc_pairs(.x, pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots()
# some correlation, but not bad

# 2.1.14 Prior-posterior comparison ####
kinetics_standard_rh_prior_posterior <- kinetics_standard_rh_samples %>%
  map(~ spread_draws(.x, Fmax, beta, F0, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1),
          F0_prior = rgamma(length(.draw), 500^2 / 300^2, 500 / 300^2),
          beta_prior = rgamma(length(.draw), 5e3^2 / 3e3^2, 5e3 / 3e3^2),
          Fmax_prior = rgamma(length(.draw), 24e4^2 / 5e4^2, 24e4 / 5e4^2)
        )
  )

str(kinetics_standard_rh_prior_posterior)

kinetics_standard_rh_prior_posterior %>%
  map(~ pivot_longer(., cols = c("F0", "F0_prior", "beta", "beta_prior", "Fmax", "Fmax_prior"),
                     names_to = c("Parameter", "Distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message
                     values_to = "Samples") %>%
        mutate(Parameter = fct(Parameter),
               Distribution = fct(ifelse(is.na(Distribution), # here the NAs are dealt with
                                         "posterior", Distribution))) %>%
        ggplot(aes(Samples, fill = Distribution)) +
        facet_wrap(~ Parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# priors are not restrictive

# 2.1.15 Predictions for mu ####
kinetics_standard_rh_mu <- kinetics_standard_rh_prior_posterior %>%
  map2(kinetics,
       ~ .x %>% select(Fmax, beta, F0) %>%
         expand_grid(Concentration = .y$standard %$%
                       seq(min(Concentration), max(Concentration), length.out = 50)) %>%
         mutate(mu = Fmax * beta * Concentration / ( Fmax + beta * Concentration ) + F0) %>%
         select(mu, Concentration)
  )

kinetics_standard_rh_mu_summary <- kinetics_standard_rh_mu %>%
  map(~ .x %>% group_by(Concentration) %>%
        mean_qi(mu, .width = c(.5, .8, .9)))

kinetics %>%
  map2(kinetics_standard_rh_mu_summary,
       ~ ggplot() +
            geom_point(data = .x$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .y, aes(Concentration, mu)) +
            geom_ribbon(data = .y, aes(Concentration, ymin = .lower, ymax = .upper,
                                    alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal()
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots()
# rectangular hyperbola is clearly the best!
# it only slightly overestimates the intercept when beta is low
# and underestimates it when beta is high

# 2.1.16 Compare functions ####
kinetics_standard_mu_summary <- 
  kinetics_standard_ht_mu_summary %>%
  map2(kinetics_standard_es_mu_summary, ~ bind_rows(.x, .y)) %>%
  map2(kinetics_standard_rh_mu_summary, ~ bind_rows(.x, .y) %>%
         mutate(Function = c(rep("ht", nrow(.y)), 
                             rep("es", nrow(.y)),
                             rep("rh", nrow(.y))) %>% 
                  fct())
       )

kinetics %>%
  map2(kinetics_standard_mu_summary,
       ~ ggplot() +
            geom_point(data = .x$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .y, aes(Concentration, mu, colour = Function)) +
            geom_ribbon(data = .y, aes(Concentration, ymin = .lower, ymax = .upper,
                                    alpha = factor(.width), fill = Function)) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal()
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
  map2(kinetics_standard_rh_prior_posterior,
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
          geom_violin(aes(Concentration, Activity, fill = Species, group = Well),
                      colour = NA, alpha = 0.5, position = "identity") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_minimal()
  ) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# These data are exceptionally messy, likely due to bad pipetting, and not worth analysing.
# There is no indication of saturation, with the highest enzyme activity always found at the highest
# substrate concentration. I ended up deciding on 1000 µM as the saturating substrate concentration
# due to logistical reasons.

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
    theme_minimal()
# Conclusions: 
# 1. MilliQ water (MQ) and autoclaved seawater (AS) have similar fluorescence properties.
# 2. The enzyme substrate has similar fluorescence properties in MilliQ (CMQ) and seawater (CAS).
# 3. 100% methanol (MeOH) has higher fluorescence than either water with or without substrate.
# It follows that
# 1. The standard curve, which is made up with 100% methanol, needs to be corrected using a water blank
# 2. The water reference can be either MilliQ or autoclaved seawater
# 3. The control can be either MilliQ or autoclaved seawater

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
          theme_minimal()
  ) %>%
  wrap_plots() %>% # too large to view in RStudio
  ggsave(filename = "standard_data.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)

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

experimental_standard_mod <- cmdstan_model(stan_file = write_stan_file(code = experimental_standard_stan))

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
experimental_standard_summary <- experimental_standard_samples %>%
  map(~ .x$summary())

experimental_standard_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

experimental_standard_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

experimental_standard_draws <- experimental_standard_samples %>%
  map(~ .x$draws(format = "df"))

ggsave( # too large to view in RStudio
experimental_standard_draws %>%
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top"),
filename = "standard_rank.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

experimental_standard_draws %>%
  map(~ mcmc_pairs(.x, pars = c("Fmax", "beta", "F0"))) %>%
  wrap_plots() %>%
  ggsave(filename = "standard_pairs.pdf", width = 80, height = 80, unit = "cm", device = cairo_pdf)
# some correlation between Fmax and beta, indicating some interdependence

# 4.1.4 Prior-posterior comparison ####
experimental_standard_prior_posterior <- experimental_standard_samples %>%
  map(~ spread_draws(.x, Fmax, beta, F0, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1),
          F0_prior = rgamma(length(.draw), 500^2 / 300^2, 500 / 300^2),
          beta_prior = rgamma(length(.draw), 5e3^2 / 3e3^2, 5e3 / 3e3^2),
          Fmax_prior = rgamma(length(.draw), 24e4^2 / 5e4^2, 24e4 / 5e4^2)
        )
  )

str(experimental_standard_prior_posterior)

ggsave(
experimental_standard_prior_posterior %>%
  map(~ pivot_longer(., cols = c("F0", "F0_prior", "beta", "beta_prior", "Fmax", "Fmax_prior"),
                     names_to = c("Parameter", "Distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message
                     values_to = "Samples") %>%
        mutate(Parameter = fct(Parameter),
               Distribution = fct(ifelse(is.na(Distribution), # here the NAs are dealt with
                                         "posterior", Distribution))) %>%
        ggplot(aes(Samples, fill = Distribution)) +
        facet_wrap(~ Parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top"),
filename = "standard_prior_posterior.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# priors are not restrictive

# 4.1.5 Predictions for mu ####
experimental_standard_mu <- experimental_standard_prior_posterior %>%
  map2(experimental,
       ~ .x %>% select(Fmax, beta, F0) %>%
         expand_grid(Concentration = .y$standard %$%
                       seq(min(Concentration), max(Concentration), length.out = 50)) %>%
         mutate(mu = Fmax * beta * Concentration / ( Fmax + beta * Concentration ) + F0) %>%
         select(mu, Concentration)
  )

experimental_standard_mu_summary <- experimental_standard_mu %>%
  map(~ .x %>% group_by(Concentration) %>%
        mean_qi(mu, .width = c(.5, .8, .9)))

experimental %>%
  map2(experimental_standard_mu_summary,
       ~ ggplot() +
            geom_point(data = .x$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .y, aes(Concentration, mu)) +
            geom_ribbon(data = .y, aes(Concentration, ymin = .lower, ymax = .upper,
                                    alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal()
            # + coord_cartesian(xlim = c(0, 10), ylim = c(0, 100e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "standard_mu.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# great fit

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
  ggsave(filename = "blank_data.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# there are four plates missing blanks, which will have to be substituted with an average

# 4.2.1 Prior simulation ####
# Water autofluorescence is on average 150 arbitrary units relative to maximal 100 µM standard fluorescence.
# Since fluorescence cannot be negative, a gamma distribution is best

ggplot() +
  geom_density(aes(rgamma(1e5, 150^2 / 50^2, 150 / 50^2))) +
  theme_minimal()
# seems reasonable

# 4.2.2 Run model ####
experimental_blank_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
}

parameters{
  real<lower=0> Fmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmu ~ gamma( 150^2 / 50^2, 150 / 50^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmu;
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

experimental_blank_mod <- cmdstan_model(stan_file = write_stan_file(code = experimental_blank_stan))

experimental_blank_samples <- experimental %>%
  map(~ experimental_blank_mod$sample(data = .x$blank %>%
                                            select(Fluorescence) %>%
                                            compose_data(),
                                      chains = 8,
                                      parallel_chains = parallel::detectCores(),
                                      iter_warmup = 1e4,
                                      iter_sampling = 1e4)
      )

# 4.2.3 Model checks ####
# check Rhat, effective sample size and chains
experimental_blank_summary <- experimental_blank_samples %>%
  map(~ .x$summary())

experimental_blank_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

experimental_blank_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

experimental_blank_draws <- experimental_blank_samples %>%
  map(~ .x$draws(format = "df"))

ggsave(
experimental_blank_draws %>%
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top"),
filename = "blank_rank.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 4.2.4 Prior-posterior comparison ####
experimental_blank_prior_posterior <- experimental_blank_samples %>%
  map(~ spread_draws(.x, Fmu, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1),
          Fmu_prior = rgamma(length(.draw), 150^2 / 50^2, 150 / 50^2)
        )
  )

str(experimental_blank_prior_posterior)

ggsave(
experimental_blank_prior_posterior %>%
  map(~ pivot_longer(., cols = c("Fmu", "Fmu_prior", "sigma", "sigma_prior"),
                     names_to = c("Parameter", "Distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message
                     values_to = "Samples") %>%
        mutate(Parameter = fct(Parameter),
               Distribution = fct(ifelse(is.na(Distribution), # here the NAs are dealt with
                                         "posterior", Distribution))) %>%
        ggplot(aes(Samples, fill = Distribution)) +
        facet_wrap(~ Parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top"),
filename = "blank_prior_posterior.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# priors are not restrictive

# 4.2.5 Predictions for mu ####
experimental_blank_mu <- experimental_blank_prior_posterior %>%
  map(~ .x %>% select(.chain, .iteration, .draw, Fmu))

# 4.3 Control models ####
# To correctly calculate enzyme activities of experimental samples, a substrate autogenic
# fluorescence control needs to be subtracted. In this case, it probably does not matter 
# if it is subtracted before or after conversion with the standard curve, but it is most 
# intuitive to subtract the raw fluorescence values because autogenic fluorescence does
# not correspond to fluorescence produced as a result of enzyme activity.

# Visualise
experimental %>%
  map(~ .x$control %>%
        ggplot(aes(Well, Fluorescence)) +
          geom_point(size = 3, shape = 16) +
          geom_hline(data = . %>% summarise(Fluorescence = median(Fluorescence)),
                     aes(yintercept = Fluorescence)) +
          theme_minimal()
  ) %>%
  wrap_plots() %>%
  ggsave(filename = "control_data.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# there are four plates missing blanks, which will have to be substituted with an average

# 4.3.1 Prior simulation ####
# Water autofluorescence is on average 500 arbitrary units relative to maximal 100 µM standard fluorescence.
# Since fluorescence cannot be negative, a gamma distribution is best

ggplot() +
  geom_density(aes(rgamma(1e5, 500^2 / 200^2, 500 / 200^2))) +
  theme_minimal()
# seems reasonable

# 4.3.2 Run model ####
experimental_control_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
}

parameters{
  real<lower=0> Fmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmu ~ gamma( 500^2 / 200^2, 500 / 200^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmu;
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

experimental_control_mod <- cmdstan_model(stan_file = write_stan_file(code = experimental_control_stan))

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
experimental_control_summary <- experimental_control_samples %>%
  map(~ .x$summary())

experimental_control_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

experimental_control_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

experimental_control_draws <- experimental_control_samples %>%
  map(~ .x$draws(format = "df"))

ggsave(
experimental_control_draws %>%
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top"),
filename = "control_rank.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 4.3.4 Prior-posterior comparison ####
experimental_control_prior_posterior <- experimental_control_samples %>%
  map(~ spread_draws(.x, Fmu, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1),
          Fmu_prior = rgamma(length(.draw), 500^2 / 200^2, 500 / 200^2)
        )
  )

str(experimental_control_prior_posterior)

ggsave(
  experimental_control_prior_posterior %>%
  map(~ pivot_longer(., cols = c("Fmu", "Fmu_prior", "sigma", "sigma_prior"),
                     names_to = c("Parameter", "Distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message
                     values_to = "Samples") %>%
        mutate(Parameter = fct(Parameter),
               Distribution = fct(ifelse(is.na(Distribution), # here the NAs are dealt with
                                         "posterior", Distribution))) %>%
        ggplot(aes(Samples, fill = Distribution)) +
        facet_wrap(~ Parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top"),
filename = "control_prior_posterior.pdf", width = 80, height = 40, unit = "cm", device = cairo_pdf)
# one very high posterior for Fmu, likely due to contamination of CAS
# here it is reasonable to remove the outlier






# 4.3.5 Predictions for mu ####
experimental_control_mu <- experimental_control_prior_posterior %>%
  map(~ .x %>% select(.chain, .iteration, .draw, Fmu))

# 4.4 Convert sample fluorescence ####








