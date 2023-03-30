
#------------------------------------------------
#' @title Define parameters of transmission simulation model
#'
#' @description Define the parameters that will be used when simulating data
#'  from the inbuilt SIMPLEGEN transmission model.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. \code{mu = -log(p)} unless
#'  specified otherwise.
#' @param u intrinsic incubation period. The number of days from initial
#'   inoculation to blood-stage infection in a human host.
#' @param v extrinsic incubation period. The number of days from initial
#'   inoculation to becoming infectious in a mosquito.
#' @param g lag time between human blood-stage infection and production of
#'   gametocytes.
#' @param prob_infection probability that a human becomes infected after being
#'   bitten by an infected mosquito. If a vector, then each subsequent value
#'   applies to each subsequent infectious bite.
#' @param prob_acute probability that an infection passes to the acute stage
#'   rather than going directly to the chronic stage. If a vector, then each
#'   subsequent value applies to each subsequent infection.
#' @param prob_AC probability that an infection in the acute stage passes to the
#'   chronic stage rather than recovering directly. If a vector, then each
#'   subsequent value applies to each subsequent infection.
#' @param duration_acute,duration_chronic specifies the probability distribution
#'   of duration (in days) of acute or chronic disease. Can be a vector or a
#'   list. If a list then the first element specifies the probability
#'   distribution for the first incident of disease, the second element for the
#'   second incident and so on (the final element is used for all remaining
#'   incidents). If a single vector then the same probability distribution is
#'   used for all incidents of disease. Values are automatically normalised to a
#'   proper probability mass distribution internally (i.e. a distribution that
#'   sums to 1).
#' @param detectability_microscopy_acute,detectability_microscopy_chronic
#'   specifies the probability of parasite detection via microscopy at each day
#'   since entering the acute or chronic stages. Can be a vector or a list. If a
#'   list then the first element specifies the probability distribution for the
#'   first incident of disease, the second element for the second incident and
#'   so on (the final element is used for all remaining incidents). If a single
#'   vector then the same probability distribution is used for all incidents of
#'   disease. This is not a probability distribution, meaning values each day
#'   can be between 0 and 1.
#' @param detectability_PCR_acute,detectability_PCR_chronic equivalent to
#'   `detectability_microscopy_acute` and `detectability_microscopy_chronic`,
#'   but for detection via polymerase chain reaction (PCR).
#' @param time_treatment_acute,time_treatment_chronic the probability
#'   distribution of time (in days since entering the acute or chronic stages)
#'   until requiring treatment. At this time, the host seeks treatment according
#'   to their host-specific treatment seeking parameter (see
#'   \code{treatment_seeking_alpha} and \code{treatment_seeking_beta}). Can be a
#'   vector or a list. If a list then the first element specifies the
#'   probability distribution for the first incident of disease, the second
#'   element for the second incident and so on (the final element is used for
#'   all remaining incidents). If a single vector then the same probability
#'   distribution is used for all incidents of disease. Values are automatically
#'   normalised to a proper probability mass distribution internally (i.e. a
#'   distribution that sums to 1).
#' @param treatment_seeking_mean,treatment_seeking_sd treatment seeking is
#'   modelled in two stages; first, the time at which hosts require treatment is
#'   drawn from the probability distribution in \code{time_treatment_acute} or
#'   \code{time_treatment_chronic}, and second, whether they go ahead with
#'   treatment seeking at this point depends on their personal host-specific
#'   treatment seeking probability. These latter probabilities are drawn
#'   independently for each host at birth from a Beta distribution with mean
#'   `treatment_seeking_mean` and standard deviation `treatment_seeking_sd`.
#'   Hence, these parameters describe the average propensity to seek treatment
#'   in the population, and the degree of inequality in treatment seeking. Note
#'   that `treatment_seeking_sd`` must be less than
#'   `sqrt(treatment_seeking_mean*(1 - treatment_seeking_mean))``, otherwise the
#'   Beta distribution is undefined.
#' @param duration_prophylactic vector defining the probability distribution of
#'   duration (in days) of prophylaxis following treatment.
#' @param infectivity_acute,infectivity_chronic the probability of onward
#'   infectivity to mosquitoes at each day since entering the acute or chronic
#'   states. Can be a vector or a list. If a list then the first element
#'   specifies the infectivity for the first incident of disease, the second
#'   element for the second incident and so on (the final element is used for
#'   all remaining incidents). If a single vector then the same infectivity is
#'   used for all episodes of disease.
#' @param max_inoculations maximum number of inoculations that an individual
#'   can carry simultaneously.
#' @param H vector specifying human population size in each deme.
#' @param seed_infections vector specifying the initial number of infected human
#'   hosts in each deme. Infected hosts are assumed to have just been bitten,
#'   meaning they have just entered the latent phase.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#' @param mig_mat migration matrix specifying the daily probability of human
#'   migration. The source deme is given in rows, the destination deme in
#'   columns.
#' @param life_table vector specifying the probability of death of a human host
#'   in each one-year age group. Defaults to a life table taken from Mali - see
#'   \code{?life_table_Mali} for details.
#'
#' @importFrom stats dgeom dpois
#' @export

define_epi_model_parameters <- function(project,
                                        a = 0.3,
                                        p = 0.868,
                                        mu = -log(p),
                                        u = 12,
                                        v = 10,
                                        g = 12,
                                        prob_infection = seq(0.8, 0.5, l = 50),
                                        prob_acute = seq(1, 0.1, l = 50),
                                        prob_AC = 1.0,
                                        duration_acute = dpois(1:100, 10),
                                        duration_chronic = dpois(1:500, 100),
                                        detectability_microscopy_acute = 0.8,
                                        detectability_microscopy_chronic = 0.8 * 100 * dgeom(0:500, 1 / 100),
                                        detectability_PCR_acute = 1,
                                        detectability_PCR_chronic = 1,
                                        time_treatment_acute = dgeom(1:25, 1 / 5),
                                        time_treatment_chronic = dgeom(1:100, 1/20),
                                        treatment_seeking_mean = 0.5,
                                        treatment_seeking_sd = 0.1,
                                        duration_prophylactic = dgeom(1:25, 1/5),
                                        infectivity_acute = 0.115,
                                        infectivity_chronic = 0.078,
                                        max_inoculations = 5,
                                        H = 1000,
                                        seed_infections = 100,
                                        M = 1000,
                                        mig_mat = diag(1),
                                        life_table = life_table_Mali()) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone.
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # get list of all input values, including those set by default
  all_args <- within(as.list(environment()), rm(project))
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values
  if (is.null(project$epi_model_parameters)) {
    project$epi_model_parameters <- all_args
  }
  
  # get list of only input values defined by user
  user_arg_names <- names(as.list(match.call()))
  user_arg_names <- setdiff(user_arg_names, c("", "project"))
  user_args <- all_args[user_arg_names]
  
  # overwrite parameters defined by user
  project$epi_model_parameters[user_arg_names] <- user_args
  
  # perform checks on parameters
  check_epi_model_params(project)
  
  # standardise and process parameters
  project <- process_epi_model_params(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define outputs produced from the transmission model
#'
#' @description Loads one or more dataframes into the SIMPLEGEN project that
#'   specify which outputs will be returned from the epidemiological model.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param daily a dataframe of daily outputs.
#' @param sweeps a dataframe of outputs at specific time points.
#' @param surveys a dataframe specifying surveys to be conducted on the
#'   population.
#'
#' @export

define_epi_sampling_parameters <- function(project,
                                           daily = NULL,
                                           sweeps = NULL,
                                           surveys = NULL) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # at least one input must be non-NULL
  if (is.null(daily) & is.null(sweeps) & is.null(surveys)) {
    stop("must define at least one output type")
  }
  
  # get list of all input values, including those set by default
  all_args <- within(as.list(environment()), rm(project))
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values
  if (is.null(project$epi_sampling_parameters)) {
    project$epi_sampling_parameters <- all_args
  }
  
  # get list of only input values defined by user
  user_arg_names <- names(as.list(match.call()))
  user_arg_names <- setdiff(user_arg_names, c("", "project"))
  user_args <- all_args[user_arg_names]
  
  # overwrite parameters defined by user
  project$epi_sampling_parameters[user_arg_names] <- user_args
  
  # TODO - perform checks on final parameters
  #check_epi_sampling_params(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from the inbuilt SIMPLEGEN transmission model. Model
#'   parameters are taken from the \code{epi_model_parameters} slot of the
#'   project, and the sampling strategy is taken from the
#'   \code{epi_sampling_parameters} slot. Outputs are written to the
#'   \code{epi_output} slot.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param max_time run simulation for this many days.
#' @param save_transmission_record whether to write the transmission record to
#'   file.
#' @param transmission_record_location the file path to the transmission record.
#' @param overwrite_transmission_record if \code{TRUE} the transmission record
#'   will overwrite any existing file by the same name. \code{FALSE} by default.
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100% to avoid large amounts of output
#'   being printed to markdown files.
#' @param silent whether to suppress written messages to the console.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats setNames
#' @export

sim_epi <- function(project,
                    max_time = 365,
                    save_transmission_record = FALSE,
                    transmission_record_location = "",
                    overwrite_transmission_record = FALSE,
                    pb_markdown = FALSE,
                    silent = FALSE) {
  
  # avoid "no visible binding" warning
  #numer <- denom <- measure <- NULL
  
  # ---------- check inputs ----------
  
  assert_class(project, "simplegen_project")
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_single_logical(save_transmission_record)
  assert_string(transmission_record_location)
  assert_single_logical(overwrite_transmission_record)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # optionally return warning if will overwrite transmission record file
  if (save_transmission_record & !overwrite_transmission_record) {
    if (file.exists(transmission_record_location)) {
      stop(sprintf("file already exists at %s. Change target location, or use argument `overwrite_transmission_record = TRUE` to manually override this warning", transmission_record_location))
    }
  }
  
  # ensure that parameters are loaded and pass checks
  check_epi_model_params(project)
  # TODO - check_epi_sampling_params(project)
  
  
  # ---------- define arguments  ----------
  
  # create argument list by concatenating project parameters
  args <- c(project$epi_model_parameters,
            project$epi_sampling_parameters)
  
  # append values that were arguments to this function
  args <- c(args,
            list(max_time = max_time,
                 save_transmission_record = save_transmission_record,
                 transmission_record_location = transmission_record_location,
                 pb_markdown = pb_markdown,
                 silent = silent))
  
  # split migration matrix into list
  args$mig_mat <- split(args$mig_mat, f = 1:nrow(args$mig_mat))
  
  # get complete demography distribution from life table
  demog <- get_demography(project$epi_model_parameters$life_table)
  args <- c(args, list(age_death = demog$age_death,
                       age_stable = demog$age_stable))
  
  # establish which sampling outputs are required. This is because C++ cannot
  # interpret a NULL input value, so we need to flag which input elements to
  # avoid
  args$any_daily_outputs <- !is.null(args$daily)
  args$any_sweep_outputs <- !is.null(args$sweeps)
  args$any_survey_outputs <- !is.null(args$surveys)
  
  # process sampling objects to make compatible with cpp11 format
  args$daily_name <- args$daily$name
  args$daily <- process_sweep(args$daily)
  
  # # get sampling strategy indices into 0-indexed numerical format
  # sampling_to_cpp_format <- function(x) {
  #   if (is.null(x)) {
  #     return(x)
  #   }
  #   w <- which(x$deme != -1)
  #   x$deme[w] <- x$deme[w] - 1
  #   x$measure <- match(x$measure, c("count", "prevalence", "incidence", "EIR")) - 1
  #   if ("state" %in% names(x)) {
  #     x$state <- match(x$state, c("S", "E", "A", "C", "P", "H", "Sv", "Ev", "Iv", "M")) - 1
  #   }
  #   if ("diagnostic" %in% names(x)) {
  #     x$diagnostic <- match(x$diagnostic, c("true", "microscopy", "PCR")) - 1
  #   }
  #   if ("sampling" %in% names(x)) {
  #     x$sampling <- match(x$sampling, c(NA, "ACD", "PCD")) - 1
  #   }
  #   return(x)
  # }
  # args$daily <- sampling_to_cpp_format(args$daily)
  # args$sweeps <- sampling_to_cpp_format(args$sweeps)
  # surveys_raw <- args$surveys
  # args$surveys <- sampling_to_cpp_format(args$surveys)
  # 
  # # get unique times at which sweeps happen
  # args$sweep_time_ordered <- NULL
  # if (args$any_sweep_outputs) {
  #   args$sweep_time_ordered <- sort(unique(args$sweeps$time))
  # }
  # 
  # # add columns to surveys
  # if (args$any_survey_outputs) {
  #   args$surveys$n_days <- args$surveys$t_end - args$surveys$t_start + 1
  # }
  # 
  # # get expanded version of surveys dataframe
  # args$surveys_expanded <- expand_surveys(args$surveys)
  # 
  # # functions
  # args_functions <- list(update_progress = update_progress)
  # 
  # make progress bars
  pb_sim <- txtProgressBar(min = 0, max = max_time, initial = NA, style = 3)
  args_progress <- list(pb_sim = pb_sim)
  
  
  # ---------- run simulation ----------
  
  #return(args)
  
  # run efficient C++ function
  output_raw <- sim_indiv_deploy(args, args_progress)
  
  #return(output_raw)
  
  # ---------- process output ----------
  
  # wrangle daily output
  daily_output <- NULL
  if (args$any_daily_outputs) {
    daily_output <- output_raw$daily %>%
      as.data.frame() %>%
      setNames(args$daily_name) %>%
      dplyr::mutate(time = seq_len(max_time), .before = 1)
  }

  # # wrangle sweep output
  sweeps_output <- NULL
  # if (args$any_sweep_outputs) {
  #   
  #   # make output dataframe with raw values
  #   sweeps_output <- cbind(project$epi_sampling_parameters$sweeps,
  #                          value = unlist(output_raw$sweep_output),
  #                          row.names = NULL)
  # }
  # 
  # # wrangle surveys output
  sample_details <- NULL
  if (args$any_survey_outputs) {
    sample_details <- data.frame(host_ID = output_raw$host_ID) %>%
      mutate(inoc_ID = output_raw$inoc_ID)
  }
  #   
  #   # make individual-level output dataframe
  #   sample_details <- mapply(as.data.frame, output_raw$survey_output, SIMPLIFY = FALSE) %>%
  #     dplyr::bind_rows()
  #   sample_details$infection_IDs <- output_raw$survey_output_infection_IDs
  #   sample_details <- dplyr::arrange(sample_details, .data$study_ID)
  #   
  #   # get summary output from individual-level
  #   surveys_output <- get_survey_summary(sample_details,
  #                                        surveys_raw,
  #                                        args$surveys_expanded)
  #   
  # }
  # 
  # append to project
  project$epi_output <- list(daily = daily_output,
                             sweeps = sweeps_output,
                             surveys = NULL)
  
  # append sample details
  project$sample_details <- sample_details
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Prune the transmission record
#'
#' @description Reads in a saved transmission record from file. Combines this
#'   with sampling information in the project to produce a pruned version of the
#'   transmission record that contains only the events relevant to the final
#'   sample. This pruned record is saved to the project.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param transmission_record_location the file path to a transmission record
#'   already written to file.
#' @param silent whether to suppress written messages to the console.
#' 
#' @importFrom  utils read.csv
#' @export

prune_transmission_record <- function(project,
                                      transmission_record_location = "",
                                      silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_string(transmission_record_location)
  assert_single_logical(silent)
  
  # check that project contains individual-level info
  sample_details <- project$sample_details
  #indlevel_error <- "project must contain a dataframe in project$sample_details slot, and this dataframe must contain the columns 'study_ID' and 'infection_IDs'"
  #assert_dataframe(sample_details, message = indlevel_error)
  #assert_in(c("study_ID", "infection_IDs"), names(sample_details), message = indlevel_error)
  #assert_vector_pos_int(sample_details$study_ID, name = "study_ID")
  #assert_vector_pos_int(unlist(sample_details$infection_IDs), name = "infection_IDs")
  
  # check transmission record exists
  if (!file.exists(transmission_record_location)) {
    stop(sprintf("could not find file at %s", transmission_record_location))
  }
  
  # check that headers formatted correctly
  first_row <- read.csv(transmission_record_location, nrows = 1)
  required_names <- c("time", "event", "human_ID", "mosquito_ID", "child_inoc_ID", "parent_inoc_ID")
  assert_in(required_names, names(first_row), message = sprintf("transmission record must contain the following column headers: {%s}", paste(required_names, collapse = ", ")))
  
  # extract starting inoc IDs from sample info
  inoc_IDs <- unlist(project$sample_details$inoc_ID)
  if (length(inoc_IDs) == 0) {
    stop("no inoc_ID values found in sample details. Check that there is at least one malaria positive host")
  }
  
  # define argument list
  args <- list(transmission_record_location = transmission_record_location,
               inoc_IDs = inoc_IDs,
               silent = silent)
  
  # run efficient C++ code
  output_raw <- prune_transmission_record_deploy(args)
  
  #return(output_raw)
  
  # process output into dataframe. parent_infection_ID column needs to be dealt
  # with separately because each element can be a vector
  output_processed <- output_raw[-which(names(output_raw) == "parent_infection_ID")] %>%
    as.data.frame() %>%
    dplyr::mutate(parent_infection_ID = output_raw$parent_infection_ID, .after = "child_infection_ID")
  
  # reverse order so same as original transmission record
  output_processed <- output_processed[rev(seq_len(nrow(output_processed))),]
  row.names(output_processed) <- NULL
  
  # add to project
  project$pruned_record <- output_processed
  
  # return project
  invisible(project)
}

#------------------------------------------------
#' @title Define genetic parameters
#'
#' @description Defines all parameters relating to the genetic model. This
#'   includes parameters related to recombination etc. that are used when
#'   simulating relatedness between lineages, as well as parameters related to
#'   mutation and sequencing errors etc. that are used at a later stage when
#'   converting relatedness into actual genotypes.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param oocyst_distribution vector specifying the probability distribution of
#'   each number of oocysts within the mosquito midgut.
#' @param hepatocyte_distribution vector specifying the probability distribution
#'   of the number of infected hepatocytes in a human host. More broadly, this
#'   defines the number of independent draws from the oocyst products that make
#'   it into the host bloodstream upon a bite from an infectious mosquito.
#' @param alpha parameter dictating the skew of lineage densities. Small
#'   values of \code{alpha} create a large skew, and hence make it likely that
#'   an oocyst will be produced from the same parents.
#' @param r the rate of recombination. The expected number of base pairs in a
#'   single recombinant block is 1/r.
#' @param contig_lengths vector of lengths (in bp) of each contig.
#' 
#' @importFrom stats dpois
#' @export

define_genetic_params <- function(project,
                                  r = 1e-6,
                                  alpha = 1.0,
                                  oocyst_distribution = dpois(1:10, lambda = 2),
                                  hepatocyte_distribution = dpois(1:10, lambda = 5),
                                  contig_lengths = c(643292, 947102, 1060087, 1204112, 1343552,
                                                     1418244, 1501717, 1419563, 1541723, 1687655,
                                                     2038337, 2271478, 2895605, 3291871)) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # get list of all input values, including those set by default
  all_args <- within(as.list(environment()), rm(project))
  
  # get list of only input values defined by user
  user_arg_names <- names(as.list(match.call()))
  user_arg_names <- setdiff(user_arg_names, c("", "project"))
  user_args <- all_args[user_arg_names]
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values
  if (is.null(project$genetic_parameters)) {
    project$genetic_parameters <- all_args
  }
  
  # otherwise overwrite parameters defined by user
  project$genetic_parameters[user_arg_names] <- user_args
  
  # perform checks on parameters
  check_genetic_params(project$genetic_parameters)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Simulate haplotype tree
#'
#' @description Reads in the pruned transmission record from file and uses this,
#'   along with a specified genetic model, to simulate a tree relating
#'   haplotypes to one another.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

sim_haplotype_tree <- function(project,
                               silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_single_logical(silent)
  
  # check pruned record exists
  assert_non_null(project$pruned_record,
                  message = "no pruned transmission record found in project")
  
  # check sample details exist
  assert_non_null(project$sample_details,
                  message = "no sample details dataframe found in project")
  
  # check for defined genetic params
  assert_non_null(project$genetic_parameters,
                  message = "no genetic parameters defined. See ?define_genetic_params")
  
  # subset sample to those with infections
  sample_positive <- project$sample_details %>%
    dplyr::filter(mapply(length, inoc_ID) > 0)
  
  # determine whether optional extra columns have been used
  defined_densities <- ("parent_infection_density" %in% names(project$pruned_record))
  defined_deme <- ("deme" %in% names(project$pruned_record))
  
  # define arguments
  args <- append(project$genetic_parameters,
                 list(sample_human_IDs = sample_positive$host_ID,
                      sample_infection_IDs = sample_positive$inoc_ID,
                      pruned_record = project$pruned_record,
                      defined_densities = defined_densities,
                      defined_deme = defined_deme,
                      silent = silent))
  
  # run efficient C++ code
  output_raw <- sim_haplotype_tree_cpp(args)
  
  return(output_raw)
  
  # get haplotype tree into dataframe by ading columns to pruned transmission
  # record
  haplotype_tree <- project$pruned_record[output_raw$record_row + 1, ] %>%
    dplyr::mutate(child_haplo_ID = output_raw$child_haplo_ID,
                  parent_haplo_ID = output_raw$parent_haplo_ID)
  row.names(haplotype_tree) <- NULL
  
  # append haplotype IDs and densities to sample_details dataframe
  sample_details <- project$sample_details
  w <- which(sample_details$true_positive)
  haplo_ID_list <- replicate(nrow(sample_details), integer())
  haplo_ID_list[w] <- output_raw$sample_haplo_IDs
  sample_details$haplo_IDs <- haplo_ID_list
  
  # save objects back into project
  project$haplotype_tree <- haplotype_tree
  project$sample_details <- sample_details
  
  # return project
  invisible(project)
}
