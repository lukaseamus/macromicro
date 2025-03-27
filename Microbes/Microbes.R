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

# remove additional date variable
enzyme %<>% 
  map(~ .x %>% select(-Date2))
str(enzyme)

# The project ran in multiple phases, and I will work through them in order 
# 1. on 19th and 20th October 2021 enzyme kinetics were tested
# 2. on 28th October and 3rd November I ran several tests with ultrapure water and seawater,
#    with and without samples, some of which can be excluded
# 3. after 3rd November I consistently ran samples and worked with seawater blanks and controls 

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
# I instinctively favour the last function because it "sticks" to the linear phase for longest:
Fmax <- 24e4 # example parameter values
beta <- 5e3
tibble(c = seq(0, 100)) %>%
  mutate(lm = beta * c, # linear model for comparison
         rh = Fmax * beta * c / ( Fmax + beta * c ),
         es = Fmax * ( 1 - exp( -beta * c / Fmax ) ),
         ht = Fmax * tanh( beta * c / Fmax )) %>%
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
# Both parameters have to be positive, so the gamma distribution is best.

kinetics_standard_prior <- 
  tibble(n = 1:1e3,
         Fmax = rgamma(n = 1e3, shape = Fmax^2 / 5e4^2, rate = Fmax / 5e4^2), # reparameterised in terms of mean and s.d.
         beta = rgamma(n = 1e3, shape = beta^2 / 2e3^2, rate = beta / 2e3^2)) %>% 
  expand_grid(c = seq(0, 100)) %>%
  mutate(`F` = Fmax * tanh( beta * c / Fmax ))

kinetics_standard_prior %>%
  ggplot(aes(c, `F`, group = n)) +
    geom_hline(yintercept = 24e4) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()
# covers all reasonable possibilities

# 2.1.2 Run model ####
kinetics_standard_stan <- "
data{
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> Fmax;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  Fmax ~ gamma( 24e4^2 / 5e4^2, 24e4 / 5e4^2 );
  beta ~ gamma( 5e3^2 / 2e3^2, 5e3 / 2e3^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Fmax * tanh( beta * Concentration[i] / Fmax );
  }

  // Likelihood
  Fluorescence ~ normal( mu , sigma );
}
"

require(cmdstanr)
kinetics_standard_mod <- cmdstan_model(stan_file = write_stan_file(code = kinetics_standard_stan))

require(tidybayes)
kinetics_standard_samples <- kinetics %>%
  map(~ kinetics_standard_mod$sample(data = .x$standard %>%
                                       select(Fluorescence, Concentration) %>%
                                       compose_data(),
                                     chains = 8,
                                     parallel_chains = parallel::detectCores(),
                                     iter_warmup = 1e4,
                                     iter_sampling = 1e4)
      )

# 2.1.3 Model checks ####
# check Rhat, effective sample size and chains
kinetics_standard_summary <- kinetics_standard_samples %>%
  map(~ .x$summary())

kinetics_standard_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

kinetics_standard_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# great Rhat and effective sample size

kinetics_standard_draws <- kinetics_standard_samples %>%
  map(~ .x$draws(format = "df"))

require(bayesplot)
kinetics_standard_draws %>% 
  map(~ mcmc_rank_overlay(.x)) %>%
  wrap_plots() + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")
# chains look good

kinetics_standard_draws %>% 
  map(~ mcmc_pairs(.x, pars = c("Fmax", "beta"))) %>%
  wrap_plots()
# practically no correlation

# 2.1.4 Prior-posterior comparison ####
kinetics_standard_prior_posterior <- kinetics_standard_draws %>%
  map(~ spread_draws(.x, Fmax, beta, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1), 
          beta_prior = rgamma(length(.draw), 5e3^2 / 2e3^2, 5e3 / 2e3^2),
          Fmax_prior = rgamma(length(.draw), 24e4^2 / 5e4^2, 24e4 / 5e4^2)
        )
  )

str(kinetics_standard_prior_posterior)

kinetics_standard_prior_posterior %>%
  map(~ pivot_longer(., cols = c("beta", "beta_prior", "Fmax", "Fmax_prior"),
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

kinetics_standard_mu <- kinetics_standard_prior_posterior %>% 
  map2(kinetics,
       ~ select(.x, Fmax, beta) %>%
         expand_grid(Concentration = .y$standard %$% 
                       seq(min(Concentration), max(Concentration), length.out = 50)) %>%
         mutate(mu = Fmax * tanh( beta * Concentration / Fmax )) %>%
         select(mu, Concentration)
  )

kinetics_standard_mu_summary <- kinetics_standard_mu %>%
  map(~ group_by(.x, Concentration) %>%
        mean_qi(mu, .width = c(.5, .8, .9)))

kinetics %>%
  map2(kinetics_standard_mu_summary,
       ~ ggplot() +
            geom_point(data = .x$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .y, aes(Concentration, mu)) +
            geom_ribbon(data = .y, aes(Concentration, ymin = .lower, ymax = .upper,
                                    alpha = factor(.width))) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal()
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots()
# fit looks good











# Split each dataframe in enzyme into its components
enzyme %<>% map(~ list(
  water_blank = .x %>% filter(Annotation == "MQ"), # the = instead of <- is important to retain names
  seawater_blank = .x %>% filter(Annotation == "AS"),
  standard = .x %>%
    filter(str_detect(Annotation, "^[0-9.]+$")) %>% 
    rename(Concentration = Annotation) %>%
    mutate(Concentration = as.numeric(Concentration)) %>%
    filter(between(Concentration, 0, 100)), # concentrations higher than 100 µM showed a decline in fluorescence
  samples = .x %>% filter(!(Annotation %in% c("MQ", "AS", "CMQ", "CAS")) & !str_detect(Annotation, "^[0-9.]+$")),
  water_control = .x %>% filter(Annotation == "CMQ"),
  seawater_control = .x %>% filter(Annotation == "CAS")
))

str(enzyme)
enzyme$X031121enzyme1$samples # example




enzyme %>%
  map(~ .x$standard %>% pull(Concentration) %>% range())
# 031121enzyme9 does not fit in because it's range is only 0 to 1 µM













enzyme %>%
  map(~ .x$seawater_blank %>% pull(Fluorescence))

# there are NAs where AS and MQ were not measured
# these should be populated with the median value of the others,
# excluding 031121enzyme9 because it's range is not the same
water_median <- enzyme %>%
  bind_rows() %>%
  filter(Annotation %in% c("AS", "MQ"),
         Plate != "031121enzyme9") %>% # exclude plate
  mutate(Fluorescence = as.numeric(Fluorescence)) %>%
  pull(Fluorescence) %>% median

water %<>%
  map(~ ifelse(is.na(.), water_median, .)) # replace NAs

standard %<>%
  map2(water, ~ .x %>%
         mutate(Fluorescence_c = Fluorescence - (median(Fluorescence[.x$Concentration == 0]) - .y)))

# Visualise standard curves
standard %>%
  bind_rows() %>%
  ggplot(., aes(x = Concentration, y = Fluorescence_c)) +
  geom_point() +
  facet_wrap(~Plate, scales = "free") +
  labs(title = "Standard curve",
       x = "Standard concentration (µM)",
       y = "Fluorescence") +
  theme_minimal() +
  theme(axis.text.y = element_blank())

# Build standard curve models
# the Baly equation with an intercept term adequately describes the initially linear increase 
# of fluorescence with fluorophore concentration, followed by the plateauing of
# fluorescence starting at concentrations beyond the limit of linearity at ~10 µM,
# which is caused by the inner filter effect:
# F = Fmax * b * C / (Fmax + b * C) + F0,
# where F is fluorescence, Fmax is the maximal fluorescence, b is the initial
# slope F/C, C is the fluorophore concentration (µM) and F0 is the fluorescence
# at a fluorophore concentration of  0 µM (i.e. autofluorescence of the solvent)

# Prior selection
# the target value for automatic gain adjustment is 90% of the maximal value,
# which is 26000 arbitrary units, so Fmax is expected to be 0.9 * 260000
Fmax.mu <- 0.9 * 260000
Fmax.sigma <- 0.2 * Fmax.mu # 20% standard deviation
# there is little prior information about the linear slope but it can be estimated to be 
# about 2 * Fmax.mu / 100, because the plateau underestimates the linear slope
b.mu <- 2 * Fmax.mu / 100
b.sigma <- 0.2 * b.mu # 20% standard deviation
# F0 is expected to be near water_median but is necessarily > 0 because there is 
# no negative fluorescence, so log parameters for a log-normal distribution are required
F0.mu <- log(water_median)
F0.sigma <- 0.2 * F0.mu # 20% standard deviation

tibble(n = 1:1e4,
       Fmax = rnorm(n = 1e4, mean = Fmax.mu, sd = Fmax.sigma),
       b = rnorm(n = 1e4, mean = b.mu, sd = b.sigma),
       F0 = rlnorm(n = 1e4, meanlog = F0.mu, sdlog = F0.sigma)) %>% # log-normal because > 0
  expand_grid(Concentration = seq(0, 100)) %>%
  mutate(Fluorescence = Fmax * b * Concentration / (Fmax + b * Concentration) + F0) %>%
  ggplot(aes(x = Concentration, y = Fluorescence, group = n)) +
  geom_hline(yintercept = standard %>% bind_rows() %$% range(Fluorescence_c)) +
  geom_line(alpha = 0.05) +
  labs(x = "Standard concentration (µM)", y = "Fluorescence") +
  theme_minimal()
# prior simulation looks good


# Hamiltonian Monte Carlo model
require(rethinking)






enzyme_standard %>%
  map(~ lm(Fluorescence ~ Annotation, data = .x))










