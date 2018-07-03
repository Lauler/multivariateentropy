#' Calculate n-variate entropy
#'
#' @importFrom magrittr %>%
#' @param data Data frame containing only discrete data.
#' @param n_variate Numeric vector. The degrees of entropy to calculate.
#'
#' @return Returns a list of lists with each variables' n-variate entropies.
#' @export
#'
#' @examples
#' \dontrun{
#' multi_entropy(df, n_variate = 1:3)}
multi_entropy <- function(data, n_variate = 2){
  n_variate <- as.list(n_variate)
  vars <- names(data)

  data <- data %>%
    dplyr::group_by_at(vars) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    dplyr::ungroup()

  var_combinations <- lapply(n_variate,
                             function(x) dplyr::as_tibble(utils::combn(x = names(data %>% dplyr::select(-n)), m = x)))

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
#' @param data Data frame containing only discrete data.
#'
#' @return Numeric. Sum of all variables' entropies.
#' @export
sum_entropy <- function(data){

  entropies <- multi_entropy(data, n_variate = 1)

  sum_entropy <- entropies[[1]] %>%
    dplyr::summarise(sum = sum(entropy))

  return(sum_entropy$sum)
}


#' Title
#'
#' @param data Data frame with discrete variables.
#'
#' @return Merges the list of lists from the \code{multi_entropy()} function.
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


#' A table of entropy measures
#'
#' @param entropy_list Takes the output of the \code{multi_entropy()} function as input a list of data frames
#' containing univariate and bivariate entropies.
#'
#'
#' @return Returns a data frame of entropy measures: univariate entropies,
#' bivariate entropies, mutual information and dependence measures.
#' @export
#'
#' @examples
#' \dontrun{
#' df_list <- multi_entropy(df, n_variate = 1:2)
#' df_entropies <- entropy_table(df_list)}
entropy_table <- function(entropy_list){
  entropy_list <- merge_entropy_table(entropy_list)

  df <- entropy_list[[2]] %>%
    dplyr::left_join(entropy_list[[1]], by = c("var1" = "var1")) %>%
    dplyr::left_join(entropy_list[[1]], by = c("var2" = "var1")) %>%
    dplyr::rename(entropy = entropy.x,
                  entropy.x = entropy.y,
                  entropy.y = entropy) %>%
    dplyr::mutate(mutual_info = (entropy.x + entropy.y - entropy),
                  dependence_var2_var1 = (mutual_info)/entropy.x,
                  dependence_var1_var2 = (mutual_info)/entropy.y)

  return(df)
}
