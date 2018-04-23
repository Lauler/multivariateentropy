#' 4 categories.
#'
#' @param cutoffs Discrete bins.
#' @param data Numeric vector. The degrees of entropy to calculate.
#'
#' @return Data frame of discretized variables.
#' @export
#'
#' @examples \code(categorize(c(0, 0, 0.2, 0.4), df))
categorize <- function(cutoffs, data){

  vars <- names(data)

  data <- data %>%
    dplyr::mutate_at(.vars = vars(vars),
                     .funs = function(x) as.factor(dplyr::case_when(
                        x >= cutoffs[1] ~ "once a week or more",
                        x >= cutoffs[2] ~ "several times a week",
                        x > cutoffs[3] ~ "once to twice a week",
                        x <= cutoffs[4] ~ "less than once a week",
                        TRUE ~ "other"
                      )))

  return(data)
}


#' Calculate n-variate entropy
#'
#' @param data Data frame containing only discrete data.
#' @param n_variate Numeric vector. The degrees of entropy to calculate.
#'
#' @return Returns a list of lists with each variables' n-variate entropies.
#' @export
#'
#' @examples multi_entropy(df, n_variate = 1:3)
multi_entropy <- function(data, n_variate = 2){
  n_variate <- as.list(n_variate)
  vars <- names(data)

  data <- data %>%
    dplyr::group_by_at(vars) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    dplyr::ungroup()

  var_combinations <- lapply(n_variate,
                             function(x) as_tibble(utils::combn(x = names(data %>% dplyr::select(-n)), m = x)))

  entropies <- list()
  i <- 1

  for (combination in var_combinations){
    names(combination) <- lapply(combination, function(x) paste(x, collapse = ":")) %>% unlist()

    entropies[[i]] <- lapply(combination, function(x){
      data %>%
        dplyr::group_by_at(x) %>%
        dplyr::summarise(sum = sum(n)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(prop = sum/sum(sum),
                      entropies = ifelse(prop == 0,
                                         yes = 0,
                                         no = prop*log(1/prop, base = 2))) %>%
        dplyr::summarise(entropy = sum(entropies))

    }
    )

    i <- i + 1
  }

  entropies <- lapply(X = entropies,
                      FUN = function(x) dplyr::bind_rows(x, .id = "id"))

  return(entropies)
}



#' Calculate the sum of entropies from all variables
#'
#' @param cutoffs Initial values for cutoffs.
#' @param data Data frame containing only discrete data.
#'
#' @return Numeric. Sum of all variables' entropies.
#' @export
sum_entropy <- function(cutoffs, data){
  data <- categorize(cutoffs = cutoffs, data = data)

  entropies <- multi_entropy(data, n_variate = 1)

  sum_entropy <- entropies[[1]] %>%
    dplyr::summarise(sum = sum(entropy))

  return(sum_entropy$sum)
}


#' Find bin sizes which optimize the entropy
#'
#' @param data
#'
#' @return Returns an optim() object
#' @export
#'
#' @examples optimize_entropy(df)
optimize_entropy <- function(data){
  best_params <- stats::optim(par = c(0, 0, 0, 0),
                              fn = sum_entropy,
                              data = data,
                              lower = c(0, 0, 0, 0),
                              upper = c(5, 5, 5, 5),
                              method = "L-BFGS-B",
                              control = list(fnscale = -1))

  return(best_params)
}


#' Title
#'
#' @param data Data frame with discrete variables.
#'
#' @return Merges the list of lists from the \code(multi_entrop()) function.
#' @export
merge_entropy_table <- function(data){

  i <- 0
  nr_vars_start <- strsplit(x = data[[1]]$id[1], split = ":") %>% unlist() %>% length()

  for (entropy in data){
    var_names <- paste0(rep("var", nr_vars_start+i), 1:(nr_vars_start+i))

    data[[i+1]] <- entropy %>%
      tidyr::separate(id, into = var_names, sep = ":") %>%
      dplyr::arrange_at(.vars = var_names)

    i <- i + 1
  }

  return(data)
}


