#### Load and explore data ####
biomass <- read.csv("~/Desktop/Projects/QUB/Data/Decomposition/Biomass.csv")

require(tidyverse)
biomass <- biomass %>% mutate(BB = Blotted/Buoyant) # compute blotted to buoyant ratio
biomass <- biomass %>% mutate(LB = Lyophilised/Blotted) # compute lyophilised to blotted ratio

ggplot(biomass, aes(Buoyant, Blotted)) + 
  geom_point() + 
  facet_grid(~Species) # looks linear

ggplot(biomass, aes(Buoyant, BB)) + 
  geom_point() + 
  facet_grid(~Species) # one logically impossible technical outlier, where blotted-buoyant ratio > 1
biomass %>% filter(BB > 1)
# remove outlier for S. latissima
require(magrittr)
biomass %<>% filter(BB < 1 | is.na(BB))


ggplot(biomass, aes(Blotted, Lyophilised)) + 
  geom_point() + 
  facet_grid(~Species) # looks non-linear

ggplot(biomass, aes(Blotted, LB)) + 
  geom_point() + 
  facet_grid(~Species) # no technical outliers

#### Buoyant vs. Blotted ####
#### Select data ####
bb <- biomass %>% 
  filter(!is.na(Buoyant)) %>%
  dplyr::select(-c(Lyophilised, LB))

#### Prior selection ####
# this model only requires a prior for the slope (b) as we know the intercept is 0
# we also know that b must be positive, since any other relationship between different
# measurements of mass is unimaginable
# we also know that b cannot exceed 1 because it is not imaginable that for 1 unit
# increase in buoyant mass, blotted mass increases by more than 1 (see removal of outlier above)
# a lognormal prior would fulfill the first requirement but not the second, so a uniform prior
# with equal probability between 0 and 1 is best

b.prior <- runif(n = 1e5, min = 0, max = 1)
ggplot() +
  geom_density(aes(x = b.prior)) +
  theme_minimal()

predictor <- bb %$% seq(min(Buoyant), max(Buoyant), length.out = 100)
plot(NULL, xlim = c(0, 40), ylim = c(0, 20), # base plot is better with loops
     xlab = "buoyant mass", ylab = "blotted mass") # empty plot
bb %$% abline(h = min(Blotted), lty = 2) # data range
bb %$% abline(h = max(Blotted), lty = 2)
abline(a = 0, b = 1, lty = 2) # 1:1, slope that cannot be exceeded
for(i in 1:1e3) curve(b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = alpha("black", 0.1))

#### Build model ####
require(rethinking)
# prepare data for ulam
bbl <- bb %>% 
  mutate(Species = as.integer(fct(Species))) %>%
  dplyr::select(Species, Buoyant, Blotted) %>%
  as.list()

# usually partial pooling should be used to be able to predict for new species
# such a model would look like this:
# bbm <- ulam(
#   alist(
#     # likelihood
#     Blotted ~ dnorm(mean = mu, sd = sigma),
# 
#     # linear model
#     mu <- b[Species] * Buoyant,
# 
#     # higher level priors for partial pooling
#     b[Species] ~ dnorm(mean = b_, sd = sb),
# 
#     # lower level priors
#     b_ ~ dunif(min = 0, max = 0),
#     c(sigma, sb) ~ dexp(rate = 1)
#   ),
#   data = bbl, chains = 8, cores = parallel::detectCores(), iter = 1e4,
# )
# however, this causes the logical limits 0 and 1 to be disobeyed because the normal
# higher level priors are not truncated

# for that reason partial pooling will not be used in this case, which is fine because
# the question does not relate to prediction for other species, but is restricted to
# the species in the dataset

bbm <- ulam(
  alist(
    # likelihood
    Blotted ~ dnorm(mean = mu, sd = sigma),
    
    # linear model
    mu <- b[Species] * Buoyant,
    
    # priors
    vector[4]:b ~ dunif(min = 0, max = 1), # the typical syntax "b[Species] ~" does not work with dunif
    sigma ~ dexp(rate = 1)
  ),
  data = bbl, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

trankplot(bbm)
precis(bbm, depth = 2)
# looks good

#### Prior-posterior comparison ####
bbm.posterior.l <- extract.samples(bbm)
bbm.prior.l <- extract.prior(bbm, n = 2e4)

bbm.posterior <- data.frame(bbm.posterior.l) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

bbm.prior <- data.frame(bbm.prior.l) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

bbm.prior.posterior <- data.frame(rbind(bbm.prior, bbm.posterior),
                                  distribution = rep(c("prior", "posterior"), each = 50000)) %>%
  mutate(distribution = factor(distribution))

ggplot(data = bbm.prior.posterior) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_wrap(~ coefficient + level, scales = "free") +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")
# looks good

#### Posterior prediction ####
#### Pairwise comparison ####
bbm.diff <- bbm.posterior %>%
  filter(coefficient == "b") %>%
  pivot_wider(names_from = c(coefficient, level), values_from = sample) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4) %>%
  dplyr::select(-c(b_1, b_2, b_3, b_4)) %>%
  pivot_longer(cols = starts_with("b"), names_to = "level", values_to = "diff") %>%
  arrange(level)

bbm.comp <- bbm.diff %>%
  group_by(level) %>%
  summarise(mean = mean(diff),
            sd = sd(diff),
            p = sum(diff < 0)/n(),
            invp = 1 - p)

#### Coefficients ####
bbm.coef <- bbm.posterior %>%
  filter(coefficient == "b") %>%
  group_by(level) %>%
  summarise(mean = mean(sample),
            sd = sd(sample),
            mean.r = signif(mean, digits = 2),
            sd.r = signif(sd, digits = 2))
bbm.coef
# Laminaria digitata: 0.87 ± 0.0066
# Laminaria hyperborea: 0.85 ± 0.005
# Saccharina latissima: 0.91 ± 0.0071
# Saccorhiza polyschides: 0.91 ± 0.0072

#### Lines and intervals ####
# create plot annotation dataframe
bb.annotation <- bb %>%
  group_by(Species) %>%
  summarise(n = n(),
            min = min(Buoyant),
            max = max(Buoyant))
bb.annotation
# n = 70 for all species

bb.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 70"),
                            expression(italic("n ")*"= 70"),
                            expression(italic("n ")*"= 70"),
                            expression(italic("n ")*"= 70"))),
         equation = as.character(c(expression(italic("y ")*"= 0.87"*italic("x")),
                                   expression(italic("y ")*"= 0.85"*italic("x")),
                                   expression(italic("y ")*"= 0.91"*italic("x")),
                                   expression(italic("y ")*"= 0.91"*italic("x")))))

# create prediction dataframe
predictor <- bb.annotation %$%
  data.frame(Buoyant = c(seq(min[1], max[1], length.out = 200),
                         seq(min[2], max[2], length.out = 200),
                         seq(min[3], max[3], length.out = 200),
                         seq(min[4], max[4], length.out = 200)),
             Species = rep(1:4, each = 200))

# extract posterior across prediction range
bbm.mu <- link(bbm, data = predictor)

# extract simulated prediction interval across prediction range
bbm.sigma <- sim(bbm, data = predictor)

# summarise posterior probability distributions
predictor %<>%
  mutate(bbm.mu.mean = apply(bbm.mu, 2, mean),
         bbm.mu.pi.lwr = t(apply(bbm.mu, 2, PI, prob = 0.99))[,1],
         bbm.mu.pi.upr = t(apply(bbm.mu, 2, PI, prob = 0.99))[,2],
         bbm.sigma.pi.lwr = t(apply(bbm.sigma, 2, PI, prob = 0.99))[,1],
         bbm.sigma.pi.upr = t(apply(bbm.sigma, 2, PI, prob = 0.99))[,2],
         Species = ifelse(Species == 1, "Laminaria digitata",
                          ifelse(Species == 2, "Laminaria hyperborea",
                                 ifelse(Species == 3, "Saccharina latissima",
                                        ifelse(Species == 4, "Saccorhiza polyschides", NA)))))

#### Plot ####
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_text(size = 12, face = "bold"),
                 text = element_text(family = "Helvetica Neue"))

bbp <- ggplot() +
  geom_abline(intercept = 0, slope = 1, colour = "#c9d2d7") +
  geom_line(data = predictor, aes(Buoyant, bbm.mu.mean, colour = Species)) +
  geom_ribbon(data = predictor, aes(Buoyant, ymin = bbm.mu.pi.lwr, 
                                    ymax = bbm.mu.pi.upr, fill = Species), 
              alpha = 0.5) +
  # geom_ribbon(data = predictor, aes(Buoyant, ymin = bbm.sigma.pi.lwr,
  #                                   ymax = bbm.sigma.pi.upr, fill = Species),
  #             alpha = 0.4) + # prediction interval produces illogical values
  geom_point(data = bb, aes(Buoyant, Blotted, colour = Species),
             size = 2, shape = 16, alpha = 0.3) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      guide = "none") +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    guide = "none") +
  facet_grid(~Species) +
  geom_text(data = bb.annotation, aes(0.58, 15.5, label = n),
            family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
  geom_text(data = bb.annotation, aes(0.58, 14, label = equation),
            family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
  labs(y = "Blotted mass (g)",
       x = "Buoyant mass (g)") +
  scale_x_continuous(breaks = seq(0, 18, by = 6)) +
  coord_cartesian(xlim = c(0, 18), ylim = c(0, 16), expand = FALSE) +
  mytheme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, hjust = 0,
                                  face = "italic"),
        panel.spacing = unit(.5, "cm"))

bbp # dimensions: 3 x 9 in



### could incorporate Date to add seasonal/individual variation

