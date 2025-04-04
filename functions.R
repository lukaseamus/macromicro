prior_samples <- function(model, data, chains, samples) {
  
  empty_data <- data %>% 
    purrr::modify_if(~ typeof(.x) == "integer", ~ 0L, 
                     .else = ~ numeric(0))
  
  require(cmdstanr)
  prior_samples <- model$sample(
    data = empty_data,
    chains = chains,
    parallel_chains = parallel::detectCores(),
    iter_warmup = samples,
    iter_sampling = samples)
  
  return(prior_samples)
  
}

prior_posterior_draws <- function(prior_samples, posterior_samples, 
                                  group, parameters, format) {
  
  parameters_list <- parameters %>% 
    purrr::map(~ .x %>% rlang::parse_expr())
  
  if(format %in% c("long", "longer", "gather", "gathered") &
     !is.na(group)) {
    prior_posterior_draws <- prior_samples %>%
      tidybayes::recover_types(group) %>%
      tidybayes::gather_draws(!!!parameters_list) %>%
      ungroup() %>%
      mutate(distribution = "prior") %>%
      bind_rows(
        posterior_samples %>%
          tidybayes::recover_types(group) %>%
          tidybayes::gather_draws(!!!parameters_list) %>%
          ungroup() %>%
          mutate(distribution = "posterior")
      )
    } else if(format %in% c("short", "shorter", "spread") &
              !is.na(group)) {
    prior_posterior_draws <- prior_samples %>%
      tidybayes::recover_types(group) %>%
      tidybayes::spread_draws(!!!parameters_list) %>%
      ungroup() %>%
      mutate(distribution = "prior") %>%
      bind_rows(
        posterior_samples %>%
          tidybayes::recover_types(group) %>%
          tidybayes::spread_draws(!!!parameters_list) %>%
          ungroup() %>%
          mutate(distribution = "posterior")
      )
    } else if(format %in% c("long", "longer", "gather", "gathered") &
              is.na(group)) {
      prior_posterior_draws <- prior_samples %>%
        tidybayes::gather_draws(!!!parameters_list) %>%
        ungroup() %>%
        mutate(distribution = "prior") %>%
        bind_rows(
          posterior_samples %>%
            tidybayes::recover_types(group) %>%
            tidybayes::gather_draws(!!!parameters_list) %>%
            ungroup() %>%
            mutate(distribution = "posterior")
        )
    } else if(format %in% c("short", "shorter", "spread") &
              is.na(group)) {
      prior_posterior_draws <- prior_samples %>%
        tidybayes::spread_draws(!!!parameters_list) %>%
        ungroup() %>%
        mutate(distribution = "prior") %>%
        bind_rows(
          posterior_samples %>%
            tidybayes::recover_types(group) %>%
            tidybayes::spread_draws(!!!parameters_list) %>%
            ungroup() %>%
            mutate(distribution = "posterior")
        )
    } else {
    prior_posterior_draws <- "Not a valid format."
    }
  
    return(prior_posterior_draws)
  
}

prior_posterior_plot <- function(prior_posterior_draws, group) {

  if(is.na(group)) {
    prior_posterior_plot <- 
      ggplot2::ggplot(data = prior_posterior_draws,
                      aes(x = .value, alpha = distribution)) +
      geom_density(colour = NA, fill = "black") +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ .variable, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank())
  } else if(!is.na(group)) {
    prior_posterior_plot <- 
      ggplot2::ggplot(data = prior_posterior_draws,
                      aes(x = .value, y = group, alpha = distribution)) +
      ggdist::stat_slab(height = 2) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ .variable, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank())
  } else {
    prior_posterior_plot <- "Not a valid data format or group."
  }
  
  return(prior_posterior_plot)
  
}