#------------------------------------------------
#' Measures of annual EIR and prevalence (Beier et al.,1999)
#' 
#' Data from Beier et al.(1999) consisting of measures of the prevalence of Plasmodium
#' falciparum detected by microscopy in children under five years old and annual
#' Entomological Innoculation Rate (EIR) for the same year, from 31 sites in Africa. Data were read
#' from Figure 1 in the original article using WebPlotDigitizer [https://github.com/ankitrohatgi/WebPlotDigitizer].
#'
#' @docType data
#'
#' @usage data(EIR_prev_beier1999)
#'
#' @format A data frame of 28 rows and two columns. Column one, "annual_EIR" is
#'   the annual Entomological Inoculation Rate calculated for Plasmodium
#'   falciparum and column 2, "prevalence" is the parasite prevalence of
#'   plasmodium falciparum measured at the same site
#'
#' @keywords datasets
#'
#' @references Beier, John C., Gerry F. Killeen, and John I. Githure.
#'   Entomologic inoculation rates and Plasmodium falciparum malaria prevalence
#'   in Africa. The American journal of tropical medicine and hygiene 61.1
#'   (1999): 109-113. (\href{https://www.ncbi.nlm.nih.gov/pubmed/10432066}{PubMed})
#'
#' @source
#'   \href{https://www.researchgate.net/profile/John_Beier/publication/12867829_Short_report_Entomologic_inoculation_rates_and_Plasmodium_falciparum_malaria_prevalence_in_Africa/links/004635142eb595e31a000000/Short-report-Entomologic-inoculation-rates-and-Plasmodium-falciparum-malaria-prevalence-in-Africa.pdf}{Research
#'   Gate} \href{https://github.com/ankitrohatgi/WebPlotDigitizer}{Data extractor hhg}
#' 
#' @examples
#' data(EIR_prev_beier1999)
#' plot( EIR_prev_beier1999$annual_EIR, EIR_prev_beier1999$prevalence,
#' xlab = "Annual EIR (Pf)",
#' ylab = "Parasite prevalence (Pf)",
#' main = "Annual EIR vs Prevalence, Beier et al., 1999")
#' 
"EIR_prev_beier1999"


#------------------------------------------------
#' Measures of annual EIR and prevalence within the same site
#'
#' Data from Hay et al. (2005) of measures of the prevalence and annual Entomological
#' Innoculation Rate (APfEIR) of P. falciparum in Africa between 1984 and 2004 . Data were read from Supplement S2.
#' APfEIR was defined as falciparum infected bites per adult per night indoors, using human biting rates that were averaged over 
#' one year and standardized to human bait catch equivalents on adults.
#'
#' @docType data
#'
#' @usage data(EIR_prev_hay2005)
#'
#' @format A data frame of 130 rows and two columns. Column one, "annual_EIR" is
#'   the annual Entomological Inoculation Rate calculated for Plasmodium
#'   falciparum and column 2, "prevalence" is the parasite prevalence of
#'   plasmodium falciparum measured at the same site
#'
#' @keywords datasets
#'
#' @references Hay, S., Guerra, C., Tatem, A. et al. Urbanization, malaria
#'   transmission and disease burden in Africa. Nat Rev Microbiol 3, 81–90
#'   (2005) doi:10.1038/nrmicro1069 
#'   (\href{https://www.ncbi.nlm.nih.gov/pubmed/15608702}{PubMed}
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3130901/bin/NIHMS35895-supplement-S2.pdf}{Supplement S2}
#'
#' @examples
#' data(EIR_prev_hay2005)
#' plot( EIR_prev_hay2005$annual_EIR, EIR_prev_hay2005$prevalence,
#' xlab = "Annual EIR (Pf)",
#' ylab = "Parasite prevalence (Pf)",
#' main = "Annual EIR vs Prevalence, Hay et al., 2005")
#' 
#' 
"EIR_prev_hay2005"

#------------------------------------------------
#' Measures of prevalence by PCR and Microscopy 
#'
#' Data on Plasmodium falciparum prevalence by microscopy vs. by PCR.
#' Taken from the systematic review by Whittaker et al. (Lancet Microbe, 2021),
#' downloaded from the associated
#' \href{https://github.com/cwhittaker1000/Sub_Patent_Malaria_Analysis/}{GitHub
#' repository}. This dataset contains prevalence in \emph{all ages} - see
#' \code{?PCR_micro_age_whittaker2021} for values broken down by age. Note, the
#' original repository contains more information about where these values come
#' from, and also more columns that were not considered relevant for our purposes.
#' 
#' @docType data
#'
#' @usage data(PCR_micro_full_whittaker2021)
#'
#' @format A data frame of 387 rows and 19 columns. Columns give the name and
#'   year of the study, the study location at multiple geographic levels,
#'   measures of historical and current transmission intensity in terms of
#'   regional prevalence, the numerator, denominator and resulting prevalence by
#'   both PCR and microscopy, details of the PCR and microscopy methods, and
#'   finally the sampling season. Some fields are given in terms of "Raw" values
#'   in addition to processed values. Raw values have more flexibility and are
#'   closer to the description in the original papers, while processed values
#'   are simplified to a smaller number of possible levels.
#'
#' @keywords datasets
#'
#' @references Charles Whittaker, Hannah Slater, Rebecca Nash, Teun Bousema,
#'   Chris Drakeley, Azra C. Ghani, and Lucy C. Okell. Global patterns of
#'   submicroscopic Plasmodium falciparum malaria infection: insights from a
#'   systematic review and meta-analysis of population surveys. Lancet Microbe
#'   (2021) https://doi.org/10.1016/S2666-5247(21)00055-0
#'   (\href{https://www.sciencedirect.com/science/article/pii/S2666524721000550}{ScienceDirect})
#'
#' @source \href{https://github.com/cwhittaker1000/submicroscopic_malaria/blob/master/Data/SI_Systematic_Review_Results_R_Import.csv}{Systematic Review Results}
#' 
#' 
"PCR_micro_full_whittaker2021"

#------------------------------------------------
#' Measures of prevalence by PCR and Microscopy, broken down by age
#'
#' Age-dissagregated data on Plasmodium falciparum prevalence by microscopy vs.
#' by PCR. Taken from the systematic review by Whittaker et al. (Lancet Microbe,
#' 2021), downloaded from the associated
#' \href{https://github.com/cwhittaker1000/Sub_Patent_Malaria_Analysis/}{GitHub
#' repository}. See \code{?PCR_micro_full_whittaker2021} for more information
#' about this and related datasets.
#' 
#' @docType data
#'
#' @usage data(PCR_micro_age_whittaker2021)
#'
#' @format A data frame of 164 rows and 21 columns. Columns are identical to the
#'   \code{PCR_micro_full_whittaker2021} dataset, with the addition of two
#'   columns giving the raw and processed age groups corresponding to each
#'   prevalence estimate.
#'
#' @keywords datasets
#'
#' @references Charles Whittaker, Hannah Slater, Rebecca Nash, Teun Bousema,
#'   Chris Drakeley, Azra C. Ghani, and Lucy C. Okell. Global patterns of
#'   submicroscopic Plasmodium falciparum malaria infection: insights from a
#'   systematic review and meta-analysis of population surveys. Lancet Microbe
#'   (2021) https://doi.org/10.1016/S2666-5247(21)00055-0
#'   (\href{https://www.sciencedirect.com/science/article/pii/S2666524721000550}{ScienceDirect})
#'
#' @source \href{https://github.com/cwhittaker1000/submicroscopic_malaria/blob/master/Data/SI_Systematic_Review_Results_R_Import.csv}{Systematic Review Results}
#' 
#' 
"PCR_micro_age_whittaker2021"

#------------------------------------------------
#' Data on prevalence and incidence broken down by age, used in Griffin et al. 2014 model fit
#' 
#' Data used to fit the model described in Griffin et al. (2010) and Griffin et
#' al. (2014). Data consists of multiple prevalence and incidence estimates,
#' broken down by age and study site. Alongside observation counts (i.e.
#' numerator and denominator values) data contains estimated distributions of
#' EIR and treatment rates that can be used to define priors in model fitting.
#' 
#'  @docType data
#'  
#'  @usage data(prev_inc_griffin2014)
#'
#' @format A dataframe of 262 rows and 15 columns. Columns are defined as follows:
#'   \itemize{
#'     \item \code{country}: country in which study was conducted.
#'     \item \code{reference}: the original reference for this data. All references can be found in Griffin et al. (2010) and Griffin et al. (2014).
#'     \item \code{site_name}: name of the study site.
#'     \item \code{site_index}: a unique numerical index given to each site (sites are contained within studies).
#'     \item \code{age0}, \code{age1}: lower and upper ends of this age group. Brackets are open on the right, for example if \code{age0 = 1} and \code{age1 = 2} then this corresponds to individuals of one year of age (rather than one and two year olds).
#'     \item \code{type}: whether incidence or prevalence estimate.
#'     \item \code{numer}, \code{denom}: numerator and denominator values. If prevalence data then denominator is total population size, if incidence data then denominator is total time at risk.
#'     \item \code{case_detection}: method of case detection, for example active case detection (ACD) vs. passive case detection (PCD).
#'     \item \code{meanEIR}, \code{sd_hi}, \code{sd_low}: estimate of the EIR in the region. Estimates are given in terms of a mean estimate and a lower and upper 95% interval. These values case be used to define a prior distribution on the EIR.
#'     \item \code{alpha(prop treated)}, \code{beta(prop treated)}: shape parameters of a Beta prior on treatment rates in the population. These can be used as a prior in model fitting, or the value \code{alpha / (alpha + beta)} can be used as a point estimate of the proportion of cases receiving treatment.
#'   }
#'
#'  @keywords datasets
#' 
#' @references Jamie T. Griffin, Deirdra Hollingsworth, Lucy C. Okell, Thomas S.
#'   Churcher, Michael White, Wes Hinsley, Teun Bousema, Chris J. Drakeley, Neil
#'   M. Ferguson, Maria-Gloria Basanes and Azra C. Ghani. Reducing Plasmodium
#'   falciparum Malaria Transmission in Africa: A Model-Based Evaluation of
#'   Intervention Strategies. PLoS Medicine (2010)
#'   doi:10.1371/journal.pmed.1000324
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/20711482/}{PubMed})
#'   
#'   Jamie T. Griffin, Neil M. Ferguson and Azra C. Ghani. Estimates of the
#'   changing age-burden of Plasmodium falciparum malaria disease in sub-Saharan
#'   Africa. Nature Communications (2014) DOI: 10.1038/ncomms4136
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/24518518/}{PubMed})
#' 
#' 
"prev_inc_griffin2014"
