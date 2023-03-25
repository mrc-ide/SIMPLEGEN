#------------------------------------------------
#' Measures of annual EIR and prevalence (Beier et al.,1999)
#' 
#' Data from Beier et al.(1999) consisting of measures of the prevalence of
#' Plasmodium falciparum detected by microscopy in children under five years old
#' and annual Entomological Innoculation Rate (EIR) for the same year, from 31
#' sites in Africa. Data were read from Figure 1 in the original article using
#' WebPlotDigitizer [https://github.com/ankitrohatgi/WebPlotDigitizer].
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
#' Data from Hay et al. (2005) of measures of the prevalence and annual
#' Entomological Innoculation Rate (APfEIR) of P. falciparum in Africa between
#' 1984 and 2004 . Data were read from Supplement S2. APfEIR was defined as
#' falciparum infected bites per adult per night indoors, using human biting
#' rates that were averaged over one year and standardized to human bait catch
#' equivalents on adults.
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
#' Data on Plasmodium falciparum prevalence by microscopy vs. by PCR. Taken from
#' the systematic review by Whittaker et al. (Lancet Microbe, 2021), downloaded
#' from the associated
#' \href{https://github.com/cwhittaker1000/Sub_Patent_Malaria_Analysis/}{GitHub
#' repository}. This dataset contains prevalence in \emph{all ages} - see
#' \code{?PCR_micro_age_whittaker2021} for values broken down by age. Note, the
#' original repository contains more information about where these values come
#' from, and also more columns that were not considered relevant for our
#' purposes.
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
#' Data on prevalence and incidence broken down by age, used in Griffin et al.
#' 2014 model fit
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
#' @format A dataframe of 262 rows and 16 columns. Columns are defined as follows:
#'   \itemize{
#'     \item \code{country}: country in which study was conducted.
#'     \item \code{reference}: the original reference for this data. All references can be found in Griffin et al. (2010) and Griffin et al. (2014).
#'     \item \code{site_name}: name of the study site.
#'     \item \code{site_index}: a unique numerical index given to each site (sites are contained within studies).
#'     \item \code{age0}, \code{age1}: lower and upper ends of this age group. Brackets are open on the right, for example if \code{age0 = 1} and \code{age1 = 2} then this corresponds to individuals of one year of age (rather than one and two year olds).
#'     \item \code{type}: whether incidence or prevalence estimate.
#'     \item \code{numer}, \code{denom}: numerator and denominator values. If prevalence data then denominator is total population size, if incidence data then denominator is total time at risk.
#'     \item \code{case_detection}: method of case detection, for example active case detection (ACD) vs. passive case detection (PCD).
#'     \item \code{meanEIR}, \code{sd_hi}, \code{sd_low}: prior used when estimating the EIR in the region. Priors are given in terms of a mean and a lower and upper 95% interval. These values are only priors, and do not represent the fitted EIRs from the original paper.
#'     \item \code{alpha(prop treated)}, \code{beta(prop treated)}: shape parameters of a Beta prior on treatment rates in the population.
#'     \item \code{plot_index}: the order in which sites are plotted in the original paper. This is also the order in terms of posterior EIR.
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

#------------------------------------------------
#' 24-SNP barcode data at two time-points in Senegal (Bei et al., 2018)
#' 
#' Data from Bei et al. (2018). Here we give a brief summary of the data - see
#' the original paper for full details.
#' \cr
#' \cr
#' Samples were randomly selected from longitudinal cohorts collected from
#' Dielmo/Ndiop, Senegal, and focus on 2 distant time-points (2001-2002, and
#' 2014). The first time-point corresponds to a period of high transmission, and
#' the second a period of extremely low transmission. Samples were genotyped
#' using a 24-SNP barcode, and complexity of infection was estimated using the
#' COIL algorithm. Samples were also matched against a large database of
#' previously published and unpublished barcodes. For samples from the low
#' transmission period (2014), one of the three repeated barcode clusters (n =
#' 6) corresponded to a parasite type (haplotype 3), observed in Thiès in both
#' 2007 and 2010. The other 2 clusters (IP1, n = 2; IP2, n = 2) had not
#' previously been observed. Data were extracted using Tabula v1.2.1.
#'
#' @docType data
#'
#' @usage data(Bei_2018)
#'
#' @format A list of multiple data objects:
#'   \itemize{
#'     \item \code{EIR}: estimates of the EIR at each location and year
#'     \item \code{barcodes}: genetic data and associated sample characteristics
#'     \item \code{SNP_locations}: genomic locations of SNPs (key to be used with \code{barcodes})
#'   }
#'   
#'   \code{EIR}: A dataframe with 3 columns, giving the time, the sampling
#'   location, and the estimated EIR (see original paper for details). The EIR
#'   in 2014 in Ndiop was recorded as "0.0" in the paper, and so has been coded
#'   as "<0.05" here to indicate the precision of this estimate.
#'   
#'   \code{barcodes}: A dataframe with 34 columns. Gives sample characteristics
#'   (columns 1:7), the estimated COI and whether this indicates a
#'   monogenomic/polygenomic infection (columns 8:9), the individual SNPs
#'   (columns 10:33) and the corresponding haplotype, if known (column 34).
#'   
#'   \code{SNP_locations}: A dataframe that acts as a key relating the SNP codes
#'   present in \code{barcodes} to the corresponding genomic location.
#' 
#' @keywords datasets
#' 
#' @references
#' \insertRef{bei_dramatic_2018}{SIMPLEGEN}
#' 
#' @source
#'   \href{https://academic.oup.com/jid/article/217/4/622/4793403}{Data in main paper}
#'   \href{https://tabula.technology/}{Tabula data extractor}
#' 
#' @importFrom Rdpack reprompt
#' 
"Bei_2018"

#------------------------------------------------
#' 93-SNP barcode data from Thai-Myanmar border (Taylor et al., 2017)
#' 
#' Data from Taylor et al. (2017). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Samples were obtained from Shoklo Malaria Research Unit (SMRU) clinics
#' spanning approximately 120km of the border between Thailand and Myanmar, and
#' were collected between 2001-2014. 1731 samples were genotyped using a 93-SNP
#' barcode. Samples with 6 or more heteroallelic genotyping outcomes were
#' designated multiple infections (n = 558), and all other samples were denoted
#' single infections (n = 1173).
#'
#' @docType data
#'
#' @usage data(Taylor_2017)
#'
#' @format A dataframe with 98 columns, giving sample IDs, dates and locations
#'   (columns 1:3), estimated clonality (columns 4:5), and genotyping calls at
#'   all 93 SNPs (columns 6:98). Missing genetic data is denoted "N".
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{taylor_quantifying_2017}{SIMPLEGEN}
#'
#' @source
#'   \href{https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007065}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Taylor_2017"

#------------------------------------------------
#' 24-SNP barcode data from two urban cities in Nigeria (Bankole et al., 2018)
#' 
#' Data from Bankole et al. (2018). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Samples were obtained from two major urban cities, Ibadan and Enugu in
#' Nigeria, both high transmission settings. DNA was extracted from dried blood
#' impregnated filter paper and samples were genotyped using a 24-SNP barcode.
#' Samples were considered polygenomic if they had more than 2 heterozygous
#' sites (Enugu = 5/28, Ibadan = 5/37), and monogenomic otherwise.  Data were
#' extracted using Tabula v1.2.1.
#'
#' @docType data
#'
#' @usage data(Bankole_2018)
#'
#' @format A dataframe with 27 columns, giving sample IDs and locations (columns
#'   1:2), genotyping calls at all 24 SNPs (columns 3:26), and (unknown) (column
#'   27). Heterozygous genotyping calls are denoted "N", and missing data is
#'   denoted "X".
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{bankole_characterization_2018}{SIMPLEGEN}
#'
#' @source
#'   \href{https://malariajournal.biomedcentral.com/articles/10.1186/s12936-018-2623-8}{Data in main paper}
#'   
#'   \href{https://tabula.technology/}{Tabula data extractor}
#' 
#' @importFrom Rdpack reprompt
#' 
"Bankole_2018"

#------------------------------------------------
#' 24-SNP barcode data from 11 sentinel sites in Haiti (Daniels et al., 2020)
#' 
#' Data from Daniels et al. (2020). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Sample sites in Haiti included four sentinel sites in three departments
#' (Grand’Anse, Sud, and Nippes) from among 11 sites nationwide established for
#' anti-malarial molecular resistance marker surveillance. Individuals of all
#' ages seeking treatment at clinics located at the sentinel sites with symptoms
#' of malaria who also tested positive for malaria by either microscopy or rapid
#' diagnostic test (RDT) from March 2016 to December 2017 were considered
#' eligible to participate. Interviews were also conducted to collect patient
#' information (including age, gender, recent travel history). Genomic DNA was
#' isolated from dried blood spots and genotyped using a 24-SNP barcode. Samples
#' were designated polygenomic if multiple alleles were observed at two or more
#' positions. After quality filtering stages, what remained was 42 polygenomic
#' samples and 462 monogenomic samples.
#'
#' @docType data
#'
#' @usage data(Daniels_2020)
#'
#' @format A dataframe with 14 columns, giving sample characteristics (columns
#'   1:7), barcode information (columns 8:9), whether the sample was designated
#'   mono/polygenomic (column 10) and information on travel history (columns
#'   11:14).
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{daniels_genetic_2020-1}{SIMPLEGEN}
#'
#' @source
#'   \href{https://malariajournal.biomedcentral.com/articles/10.1186/s12936-020-03439-7}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Daniels_2020"

#------------------------------------------------
#' 24-SNP barcode data from two trials in Zambia (Daniels et al., 2020)
#' 
#' Data from Daniels et al. (2020). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Rapid diagnostic test-positive samples were obtained from household-level
#' data collection during two community-randomized trials conducted in the same
#' geographical region of Southern Province of Zambia between 2012 and 2016. the
#' first sample set (n = 836 children younger than 6 years) was collected during
#' the peak malaria transmission season (April–May) in both 2012 and 2013
#' (baseline) as part of a community-randomized controlled trial designed to
#' assess the impact of three rounds of an MTAT intervention that used RDTs for
#' testing and artemether–lumefantrine for treatment. The second sample set (n =
#' 784 individuals 3 months or older) was obtained from an 18-month longitudinal
#' cohort study, with monthly follow-up visits, conducted from December 2014 to
#' May 2016 (cohort). The cohort was designed to evaluate a cluster randomized
#' controlled trial for assessing the impact of four rounds of community-wide
#' MDA and household-level (focal) MDA (fMDA) with DHAp compared with that of no
#' mass treatment (control).
#' \cr
#' Samples were genotyped using a 24-SNP barcode. Samples were designated
#' polygenomic if multiple alleles were observed at two or more positions,
#' otherwise they were designated monogenomic. Samples with missing data at 5 or
#' more loci were deemed to have "failed" for the purposes of subsequent
#' analyses, but are included in the data anyway.
#'
#' @docType data
#'
#' @usage data(Daniels_2020B)
#'
#' @format A dataframe with 34 columns, giving the sample ID, collection date,
#'   study type, (TODO), study batch and arm (columns 1:6), the genomic data
#'   over all 24 SNPs (columns 7:30), and details of missingness and designated
#'   mono/polyclonality (columns 31:34). Heterozygous genotyping calls are
#'   identified by "N", and missing alleles are identified by "X". Note, the
#'   original downloaded data contained 18 samples with missing genotypic data,
#'   which have been removed to leave the 1620 samples referred to in the main
#'   paper.
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{daniels_evidence_2020}{SIMPLEGEN}
#'
#' @source
#'   \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7416975/}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Daniels_2020B"

#------------------------------------------------
#' 24-SNP barcode data from Richard Toll, Senegal (Daniels et al., 2020)
#' 
#' Data from Daniels et al. (2020). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Samples were obtained from routine case investigation carried out in Richard
#' Toll, Senegal between September 2012 and December 2015. Rapid diagnostic
#' tests (RDTs) were used to diagnose malaria cases either through
#' facility-based passive case detection (PCD) or through reactive case
#' detection (RACD). A standardized questionnaire was also filled out for all
#' participants to collect information on basic demographic information
#' including travel history. RDTs were used to genotype malaria infections using
#' a 24-SNP barcode. Samples were designated polygenomic if multiple alleles
#' were observed at two or more positions, otherwise they were designated
#' monogenomic. Samples with missing data at 5 or more loci were deemed to have
#' "failed" for the purposes of subsequent analyses, but are included in the
#' data anyway.
#'
#' @docType data
#'
#' @usage data(Daniels_2020C)
#'
#' @format A dataframe with 30 columns, giving sample ID and year (columns 1:2),
#'   genomic data at 24 SNPs (columns 3:26), and details of missingness and
#'   designated mono/polyclonality (columns 27:30). Heterozygous genotyping
#'   calls are identified by "N", and missing alleles are identified by "X".
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{daniels_genetic_2020}{SIMPLEGEN}
#'
#' @source
#'   \href{https://malariajournal.biomedcentral.com/articles/10.1186/s12936-020-03346-x}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Daniels_2020C"

#------------------------------------------------
#' 24-SNP barcode data from clonal outbreak in Panama (Obaldia et al., 2015)
#' 
#' Data from Obaldia et al. (2015). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' This study focused on 37 P. falciparum isolates collected during 2003–2008
#' from individuals in malaria-endemic provinces in eastern Panama, along with
#' 20 isolates collected during 2011–2012 from healthcare facilities in 3
#' malaria endemic sites in Colombia. Samples were sequenced using a 24-SNP
#' barcode. Samples were considered polyclonal if they contained multiple
#' alleles at two or more loci, otherwise they were considered monoclonal.
#' Samples were clustered into major groups using STRUCTURE software.
#'
#' @docType data
#'
#' @usage data(Obaldia_2015)
#'
#' @format A dataframe with 30 columns, giving sample characteristics (columns
#'   1:5), genomic data at 24 SNPs (columns 6:29), and major group membership
#'   (column 30). Heterozygous genotyping calls are identified by "N", and
#'   missing alleles are identified by "NA".
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{obaldia_clonal_2015}{SIMPLEGEN}
#'
#' @source
#'   \href{https://academic.oup.com/jid/article-lookup/doi/10.1093/infdis/jiu575}{Data extracted manually from paper}
#' 
#' @importFrom Rdpack reprompt
#' 
"Obaldia_2015"

#------------------------------------------------
#' 24-SNP barcode data showing population structure in Haiti (Charles et al.,
#' 2016)
#' 
#' Data from Charles et al. (2016). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Samples were collected through routine surveillance spanning 2006–2009 by the
#' Haitian Group for the Study of Kaposi’s Sarcoma and Opportunistic Infections
#' (GHESKIO) at 9 healthcare centers in various municipalities. Samples were
#' genotyped using a 24-SNP barcode. The paper states that samples that showed
#' >1 mixed-base SNP call or had >5 missing calls in the 24SNP molecular barcode
#' were removed from analysis, however, there is 1 sample (monogenomic sample
#' 24) that has 6 missing calls and has not been removed, and there is also 1
#' sample (polygenomic sample 8) that has exactly 1 mixed-base SNP call and 0
#' missing calls and yet has been removed. Therefore, it is possible that
#' samples were removed if they showed >0 mixed-base calls, or had >6 missing
#' calls (all data presented are consistent with this filtering). Both
#' monogenomic (n = 52) and polygenomic (n = 8) samples are available, although
#' samples excluded based on missingness are not available. Samples were
#' compared in terms of the proportion of shared bases, and were identified as
#' either identical, related, or unique based on this number.
#'
#' @docType data
#'
#' @usage data(Charles_2016)
#'
#' @format A list of two data objects:
#'   \itemize{
#'     \item \code{monoclonal}: the 42 samples identified as monoclonal
#'     \item \code{polyclonal}: the remaining 8 samples identified as polyclonal
#'   }
#'   
#'   \code{monoclonal}: A dataframe with 6 columns, giving sample
#'   characteristics (columns 1:3), barcode data (columns 4:5), and similarity
#'   category (column 6). location, and the estimated EIR (see original paper
#'   for details). Heterozygous genotyping calls are identified by "N", and
#'   missing alleles are identified by "X".
#'   
#'   \code{polyclonal}: A dataframe with 4 columns, giving sample
#'   characteristics (columns 1:3), and barcode data (column 4). As above,
#'   heterozygous genotyping calls are identified by "N", and missing alleles
#'   are identified by "X".
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{charles_plasmodium_2016}{SIMPLEGEN}
#'
#' @source
#'   \href{http://wwwnc.cdc.gov/eid/article/22/5/15-0359_article.htm}{Data extracted from tables}
#' 
#' @importFrom Rdpack reprompt
#' 
"Charles_2016"

#------------------------------------------------
#' 24-SNP barcode data showing clonal and epidemic transmission in Senegal
#' (Daniels et al., 2013)
#' 
#' Data from Daniels et al. (2013). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' Note - see \code{?Galinsky_2015} for an expanded version of this same set of
#' barcodes. All 529 Sample IDs in this dataset are also present in
#' Galinsky_2015, along with an additional 445 samples.
#' \cr
#' \cr
#' Samples were obtained annually from 2006–2011 from the Service de Lutte
#' Anti-Parasitaire (SLAP) clinic in Thies, Senegal. Samples were collected
#' passively; with patients over the age of 12 months admitted to this study
#' with self-reported acute fevers within 24 hours of visiting the clinic and no
#' recent anti-malarial use. Whole blood spots from 2006–2011 were preserved on
#' filter paper, and DNA was extracted and characterised using a 24-SNP
#' molecular barcode. Samples were excluded if they had missing data on more
#' than four SNP positions. Samples were identified as polygenomic if they had
#' multiallelic calls at more than one SNP position, otherwise monogenomic. It
#' is unclear whether the data presented is for monogenomic samples only, as one
#' sample (SenT081) has mixed calls at 2 sites and therefore is polygenomic, but
#' it seems unlikely that there was only one polygenomic sample observed in the
#' entire dataset.
#' 
#'
#' @docType data
#'
#' @usage data(Daniels_2013)
#'
#' @format A dataframe with 27 columns, giving sample characteristics (columns
#'   1:3), and genomic data at 24 SNPs (columns 4:27). Heterozygous genotyping
#'   calls are identified by "N", and missing alleles are identified by "NA".
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{daniels_genetic_2013}{SIMPLEGEN}
#'
#' @source
#'   \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0060780}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Daniels_2013"

#------------------------------------------------
#'  Multiplicity of infection in Mali using 24-SNP barcode data (Adomako-Ankomah
#'  et al., 2017)
#' 
#' Data from Adomako-Ankomah et al. (2017). Here we give a brief summary of the
#' data - see the original paper for full details.
#' \cr
#' \cr
#' A longitudinal cohort study was conducted in a seasonal and high-transmission
#' area of Mali, in which 500 subjects aged 1–65 years were followed for 1 year.
#' Blood samples were collected every 2 weeks. Multiclonality of Pf infection was
#' measured using a 24-SNP DNA barcoding assay at 4 time-points (two in wet season,
#' and two in dry season).
#' 
#'
#' @docType data
#'
#' @usage data(Adomako_Ankomah_2017)
#'
#' @format A dataframe with 26 columns, giving time point in the year (column
#'   1), multiplicity category, including number of failed assays (column 2),
#'   and frequencies for each of these classes for all 24 SNPs (columns 3:26).
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{adomako-ankomah_high_2017}{SIMPLEGEN}
#'
#' @source
#'   \href{https://dx.plos.org/10.1371/journal.pone.0170948}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Adomako_Ankomah_2017"

#------------------------------------------------
#' 24-SNP barcode data showing monogenomic and polygenomic infections in Senegal
#' (Galinsky et al., 2015)
#' 
#' Data from Galinsky et al. (2015). Here we give a brief summary of the
#' data - see the original paper for full details.
#' \cr
#' \cr
#' Samples are an expanded set of those present in another dataset (see
#' \code{?Daniels_2013} for further details of sample collection). The
#' difference is that samples here also include those that are likely
#' polygenomic infections.
#'
#' @docType data
#'
#' @usage data(Galinsky_2015)
#'
#' @format A dataframe with 26 columns, giving the sample ID and collection year
#'   (columns 1:2), and genotype calls at all 24 SNP loci (columns 3:26).
#'   Heterozygous genotyping calls are identified by "N", and missing alleles
#'   are identified by "X".
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{galinsky_coil_2015}{SIMPLEGEN}
#'
#' @source
#'   \href{https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-14-4}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Galinsky_2015"

#------------------------------------------------
#' SNP data collected from Uganda (Chang et al., 2017)
#' 
#' Data from Chang et al. (2017). Here we give a brief summary of the data - see
#' the original paper for full details.
#' \cr
#' \cr
#' Dried blood spots were taken from 2012-13 cross sectional surveys in several
#' provinces in Uganda. Households (n = 200) were randomly selected from each
#' province. All samples that had detectable asexual parasitemia were selected
#' for Sequenom SNP genotyping. Medium to high frequency SNPs (n = 128) with
#' high frequency in malaria populations (pf-community-project) were chosen,
#' leaving 105 SNPs after filtering variants with lower or missing frequencies.
#' Genotyping was based on the intensity of the SNPs. In addition, merozoite
#' surface protein 2 (msp2) genotyping was conducted on an age stratified subset
#' of the samples. Capillary electrophoresis was used to distinguish msp2 allele
#' sizes.
#' \cr
#' \cr
#' The data contains a list (size = 3) of data frames. The first data frame
#' contains SNP calls for 105 filtered SNPs. The values \{-1, 0, 0.5, 1\} denote
#' missing value / no call, heterozygous, or homozygous alleles. The second and
#' third data frame contain 95% credible intervals for COI and allele
#' frequencies using their described categorical and proportional methods for
#' modeling homozygous/heterozygous calls and with-in host allele frequency,
#' respectively.
#' 
#' @docType data
#' 
#' @usage data(Chang_2017)
#' 
#' @format A list of multiple data objects
#'   \itemize{
#'     \item \code{S1_Table_SNP_data}: SNP calls for each sample \{-1, 0, 0.5, 1\}.
#'     \item \code{95_cred_interval_coi}: Calculated COI of Uganda Samples
#'     using categorical and proportional methods described in paper.
#'     \item \code{95_cred_interval_allele_freq}: Calculated allele frequency
#'     of Uganda Samples using categorical and proportional methods described
#'     in paper.
#'   }
#'   
#'  \code{S1_Table_SNP_data}: A data frame of 107 columns. The sample id and 
#'  location are included in the first 2 columns. SNP calls (columns 3:107) are 
#'  either no call / missing data (-1), heterozygous (0.5) or homozygous (0 or 
#'  1).
#'  
#'  \code{S2_Table_95_cred_int_coi_uganda_samples}: A data frame of 15 columns 
#'  containing calculated complexity of infection summary statistics for
#'  categorical, proportional, and COIL methods.
#'  
#'  \code{S3_Table_95_cred_int_allele_freq}: A data frame of 25 columns
#'  containing summary statistics of allele frequency from categorical, 
#'  proportional, and COIL methods.
#'   
#' @keywords datasets
#'
#' @references 
#'   \insertRef{chang_real_2017}{SIMPLEGEN}
#' 
#' @source 
#'   \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5300274/}{Supplementary materials}
#'   
#' @importFrom Rdpack reprompt
#' 
"Chang_2017"

#------------------------------------------------
#' Microsatellite data across 26 loci, Eswatini (Roh et al., 2019)
#' 
#' Data from Roh et al. (2019). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' The study investigates local versus imported transmissions and how that
#' effects the genetic diversity in Eswatini. The parameters collected in the
#' epidemiological model include district of residence, season, case detection,
#' age, sex, and patient occupation. To compare genetic diversity, an 
#' overlapping subset of microsatellite markers were used to compare against
#' other countries with lower transmission rates (it was found that Eswatini had
#' a greater diversity of parasites (mean $H_E$ was consistently higher)). When
#' comparing the local and imported cases, population level diversity did not
#' differentiate greatly; however, a greater portion of local cases were
#' monoclonal and less complex. Differences in MOI and $F_ws$ between cases were
#' statistically significant but had marginal effect sizes and could not be used
#' to discriminate cases. 
#' \cr
#' \cr
#' Eligible cases were identified through Eswatini's national malaria surveillance
#' program. Samples were genotyped using microsatellites from 26 different loci,
#' which showed no selection based on $F_ws$ calculations. In their study, they
#' excluded data points where the probability of there being an allele was less
#' than 95%, but these are included in the data set.
#'
#' @docType data
#'
#' @usage data(Roh_2019)
#'
#' @format A data frame with 10 columns, giving sample ID and barcode
#' (columns 1:2), locus, locus indexes, and a combined id (columns 3:5),
#' peaks and fluorescence (columns 6:9), and probability of being a true allele
#' (column 10).
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{roh_high_2019}{SIMPLEGEN}
#'
#' @source
#'   \href{https://academic.oup.com/jid/article/220/8/1346/5514501}{Supplementary materials}
#' 
#' @importFrom Rdpack reprompt
#' 
"Roh_2019"

#------------------------------------------------
#' Microsatellite data across 26 loci, Northeast Namibia (Tessema et al., 2019)
#' 
#' Data from Tessema et al. (2019). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' The study collected samples from 4643 symptomatic outpatients in northeast
#' Namibia (Kavango East and Zambezi) with confirmed cases via rapid tests. 3871
#' cases were from Kavango East (March-June 2016) from 23 different health
#' facilities, and 772 from Zambezi (Feb 2015 and June 2016) from six
#' facilities.
#' \cr
#' \cr
#' DNA was extracted from dried blood spots using punches and strips from rapid 
#' diagnostic tests. If a facility had less than 100 rapid tests, all tests were
#' included for genotyping. If over 100 tests, then all cases with travel history
#' were included for genotyping, and up to 100 cases without history were also
#' included. In total, 2990 samples were genotyped using microsatellite markers,
#' in which 2585 were used for data analysis. 
#' 
#' @docType data
#'
#' @usage data(Tessema_2019)
#'
#' @format A data frame with 31 columns, giving case ID, Health Facility,
#' District, Region, Country (Columns 1:5), and microsatellite lengths (6:31).
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{tessema_using_2019}{SIMPLEGEN}
#'
#' @source
#'   \href{https://doi.org/10.7554/eLife.43510}{Supplementary file 1}
#' 
#' @importFrom Rdpack reprompt
#' 
"Tessema_2019"

#------------------------------------------------
#' Inferring local and crossborder transmission with human and mobility data, 
#' Kanunga District in Southwest Uganda (Briggs et al., 2019)
#' 
#' Data from Briggs et al. (2021). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' The study uses 80 randomly sampled households around a single health facility
#' and includes all children (6 mo - 10 years old) and at least one adult
#' caretaker from the households. DNA was acquired though dried blood spots, and
#' after two rounds of PCR amplification, the PCR products were sized using
#' capillary electrophoresis. Finally, allele length was calculated.
#' 
#' @docType data
#'
#' @usage data(Briggs_2021)
#'
#' @format A list of multiple data objects:
#'   \itemize{
#'     \item \code{final_samples_Dec}: contains allele lengths with barcodes. 
#'     \item \code{Kanungu_pairwise_comparison_df}: epidemiological data set with
#'     \code{final_samples_Dec} keys.
#'   }
#'  \code{Kanungu_pairwise_comparison_df}: Metadata for the allele sizes
#'  located in \code{final_samples_Dec}. 
#'
#' @keywords datasets
#'
#' @references
#' \insertRef{briggs_withinhousehold_2021}{SIMPLEGEN}
#'
#' @source
#'   \href{https://malariajournal.biomedcentral.com/articles/10.1186/s12936-021-03603-7}{Supplementary file 1}
#' 
#' @importFrom Rdpack reprompt
#' 
"Briggs_2021"

#------------------------------------------------
#' 250-SNP barcode reveals connectivity on Colombian-Pacific coast (Taylor et al., 2020)
#' 
#' Data from Taylor et al. (2020). Here we give a brief summary of the data -
#' see the original paper for full details.
#' \cr
#' \cr
#' The study by Taylor et al. (2020) uses data from a previously published study
#' by Echeverry et al. (2013). Samples were obtained from patients with
#' symptomatic uncomplicated malaria, and were collected between 1993 and 2007
#' from five cities in four provinces of Colombia. Samples were genotyped using
#' a 250-SNP barcode. These samples were all considered to be monoclonal
#' infections, therefore any heterozygous genotype calls were recoded as missing
#' data. Note that markers were re-ordered post-publication, but this did not
#' qualitatively alter any of the major conclusions
#' (\href{#'https://github.com/aimeertaylor/ColombianBarcode/blob/master/README.md}{see
#' Github README}).
#'
#' @docType data
#'
#' @usage data(Taylor_2020)
#'
#' @format A dataframe with 257 columns, giving the sample ID (column 1),
#'   multi-locus genotype ID (column 2), collection place and time (columns
#'   3:5), the genotype call at all 250 SNPs (columns 6:255), the number of
#'   heterozygous loci in the original data (now recoded as missing) (column
#'   256) and the collection year (column 257, matches info in column 5).
#'   Genotype values give 0 for minor allele, 1 for major allele, or NA for
#'   missing data. Samples were considered monoclonal and therefore any
#'   heterozygous calls were recoded as missing data.
#'   
#' @keywords datasets
#'
#' @references
#' \insertRef{taylor_identity-by-descent_2020}{SIMPLEGEN}
#' 
#' \insertRef{echeverry_long_2013}{SIMPLEGEN}
#'
#' @source
#'   \href{https://github.com/aimeertaylor/ColombianBarcode}{Associated Github repository}
#' 
#' @importFrom Rdpack reprompt
#' 
"Taylor_2020"
