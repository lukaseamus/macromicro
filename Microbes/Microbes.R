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
#    some of which can be excluded
# 3. after 3rd November I consistently worked with seawater

# Enzyme kinetics
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










