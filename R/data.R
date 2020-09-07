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
#' @usage data(EIRprev_beier1999)
#'
#' @format A data frame of x rows and two columns. Column one, "annual_EIR" is
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
#' data(EIRprev_beier1999)
#' plot( EIRprev_beier1999$annual_EIR, EIRprev_beier1999$prevalence,
#' xlab = "Annual EIR (Pf)",
#' ylab = "Parasite prevalence (Pf)",
#' main = "Annual EIR vs Prevalence, Beier et al., 1999")
#' 
"EIRprev_beier1999"


#------------------------------------------------
#' Measures of annual EIR and prevalence within the same site
#'
#' Data from Hay et al. (2005) of measures of the prevalence  and annual Entomological
#' Innoculation Rate (APfEIR) of P. falciparum in Africa between 1984 and 2004 . Data were read from Supplement S2.
#' APfEIR was defined as falciparum infected bites per adult per night indoors, using human biting rates that were averaged over 
#' one year and standardized to human bait catch equivalents on adults.
#'
#' @docType data
#'
#' @usage data(EIRprev_hay2005)
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
#' data(EIRprev_hay2005)
#' plot( EIRprev_hay2005$annual_EIR, EIRprev_hay2005$prevalence,
#' xlab = "Annual EIR (Pf)",
#' ylab = "Parasite prevalence (Pf)",
#' main = "Annual EIR vs Prevalence, Hay et al., 2005")
#' 
#' 
"EIRprev_hay2005"

#------------------------------------------------
#' Measures of prevalence by PCR and Microscopy 
#'
#' Data from preprint by Whittaker et al. (BioRxiv, 2020). 
#' Data were obtained from [https://github.com/cwhittaker1000/Sub_Patent_Malaria_Analysis/]
#' 
#' @docType data
#'
#' @usage data(PCRMicro)
#'
#' @keywords datasets
#'
#' @references Global & Temporal Patterns of Submicroscopic Plasmodium falciparum Malaria Infection
#' Charles Whittaker, Hannah Slater, Teun Bousema, Chris Drakeley, Azra Ghani, Lucy Okell
#' bioRxiv 554311; doi: https://doi.org/10.1101/554311
#' https://doi.org/10.1101/554311
#'   (\href{https://www.ncbi.nlm.nih.gov/pubmed/15608702}{PubMed}
#'
#' @source \href{https://github.com/cwhittaker1000/Sub_Patent_Malaria_Analysis/blob/master/Data/SI_Systematic_Review_Results_R_Import.csv}{Systematic Review Results}
#' 
#' 
"PCRMicro"

