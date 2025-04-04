prior_samples <- function(model, data, n_name = "n", chains = 8, samples = 1e4) {
  
  empty_data <- data %>%
    purrr::modify_if(.p = ~ typeof(.x) == "double", ~ double(0)) %>%
    purrr::modify_at(.at = ~ str_match(.x, n_name), ~ 0L)
  
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
                                  group = NULL, parameters, format) {
  
  parameters_list <- parameters %>% 
    purrr::map(~ .x %>% rlang::parse_expr())
  
  if(format %in% c("long", "longer", "narrow", "narrower", "gather", "gathered") &
     !is.null(group)) {
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
    } else if(format %in% c("short", "shorter", "wide", "wider", "spread") &
              !is.null(group)) {
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
              is.null(group)) {
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
              is.null(group)) {
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

prior_posterior_plot <- function(prior_posterior_draws_long, 
                                 group_name = NULL, ridges = TRUE) {
  
  if(is.null(group_name)) {
    prior_posterior_plot <- 
      ggplot2::ggplot(data = prior_posterior_draws_long,
                      aes(x = .value, alpha = distribution)) +
                geom_density(colour = NA, fill = "black") +
                scale_alpha_manual(values = c(0.6, 0.2)) +
                facet_wrap(~ .variable, scales = "free") +
                theme_minimal() +
                theme(panel.grid = element_blank())
  } else if(!is.null(group_name) & ridges == TRUE) {
    group_name_parsed <- group_name %>% rlang::parse_expr()
    prior_posterior_plot <- 
      ggplot2::ggplot(data = prior_posterior_draws_long,
                      aes(x = .value, y = !!group_name_parsed, alpha = distribution)) +
                ggdist::stat_slab(height = 2, fill = "black") +
                scale_alpha_manual(values = c(0.6, 0.2)) +
                facet_wrap(~ .variable, scales = "free") +
                theme_minimal() +
                theme(panel.grid = element_blank())
  } else if(!is.null(group_name) & ridges == FALSE) {
    group_name_parsed <- group_name %>% rlang::parse_expr()
    prior_posterior_plot <- 
      ggplot2::ggplot(data = prior_posterior_draws_long,
                      aes(x = .value, alpha = distribution)) +
                geom_density(colour = NA, fill = "black") +
                scale_alpha_manual(values = c(0.6, 0.2)) +
                ggh4x::facet_nested_wrap(facets = vars(.variable, !!group_name_parsed), 
                                         scales = "free", nest_line = TRUE) +
                theme_minimal() +
                theme(panel.grid = element_blank())
  } else {
    prior_posterior_plot <- "Not a valid data format or group."
  }
  
  return(prior_posterior_plot)
}

spread_continuous <- function(prior_posterior_draws_short, data, predictor_name, 
                              group_name = NULL, length = 100) {
  
  predictor_name_parsed <- predictor_name %>% rlang::parse_expr()
  
  if(is.null(group_name)) {
    require(magrittr)
    spread_continuous <- prior_posterior_draws_short %>%
      expand_grid(!!predictor_name := data %$% 
                  seq(min(!!predictor_name_parsed),
                      max(!!predictor_name_parsed),
                      length.out = length))
  } else if(!is.null(group_name)) {
    group_name_parsed <- group_name %>% rlang::parse_expr()
    spread_continuous <- prior_posterior_draws_short %>%
      left_join(data %>%
                  group_by(!!group_name_parsed) %>%
                  summarise(min = min(!!predictor_name_parsed),
                            max = max(!!predictor_name_parsed)),
                by = group_name) %>%
      rowwise() %>%
      mutate(!!predictor_name := list( seq(min, max, length.out = length) )) %>%
      select(-c(min, max)) %>%
      unnest(!!predictor_name_parsed)
  } else {
    spread_continuous <- "Not a valid data format, predictor or group."
  }
  
  return(spread_continuous)
}