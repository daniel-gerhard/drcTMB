# Data sets, copied from the package drcData (Christian Ritz, Jens C. Streibig and Daniel Gerhard, GPL-2);
# see data-raw/drcData.R

#' Germination of common chickweed
#'
#' Germination of 200 chickweed (\emph{Stellaria media}) seeds in one germination test, from a study of
#' chlorsulfuron resistant and sensitive biotypes, recorded in time intervals: an example of time-to-event data (\code{family = "event"} of
#' \code{\link{drmTMB}}).
#'
#' 200 seeds were placed in Petri dishes in a dark growth cabinet at 5 degrees Celsius, and the
#' germinated seeds were counted and removed at 34 inspection times. Chickweed has dormant seeds, so
#' that the proportion of germinated seeds stays well below 1. The last row, with \code{end = Inf},
#' holds the seeds that had not germinated by the end of the experiment.
#'
#' @format A data frame with 35 rows and 3 variables:
#' \describe{
#'   \item{start}{left endpoint of the monitoring interval (hours)}
#'   \item{end}{right endpoint of the monitoring interval (hours)}
#'   \item{count}{number of seeds germinated in the interval}
#' }
#' @source Lisa Borggaard (formerly Faculty of Life Sciences, University of Copenhagen); from the
#'   package drcData.
#' @references Ritz, C., Pipper, C. B. and Streibig, J. C. (2013) Analysis of germination data from
#'   agricultural experiments. \emph{European Journal of Agronomy}, \bold{45}, 1--6.
#' @examples
#' gmod <- drmTMB(count ~ start + end, data = chickweed, model = "loglogistic", family = "event",
#'                start = c(b5 = 0), fix = c(b5 = TRUE))
#' # times by which 10%, 50% and 90% of the germinating seeds have germinated
#' ed(gmod, c(10, 50, 90))
#' # time by which 10% of all seeds have germinated
#' ed(gmod, 0.1, type = "absolute")
"chickweed"

#' Immobility of daphnids
#'
#' The number of immobile daphnids out of 20, for several concentrations of a toxic substance; the
#' same daphnids were counted after 24 and after 48 hours.
#'
#' @format A data frame with 16 rows and 4 variables:
#' \describe{
#'   \item{dose}{concentration of the toxic substance}
#'   \item{no}{number of immobile daphnids}
#'   \item{total}{total number of daphnids}
#'   \item{time}{time of the count, a factor with levels \code{24h} and \code{48h}}
#' }
#' @source Nina Cedergreen, Faculty of Life Sciences, University of Copenhagen; from the package
#'   drcData.
#' @examples
#' dmod <- drmTMB(cbind(no, total - no) ~ dose, data = daphnids, model = "loglogistic",
#'                family = "binomial", fform = b3 + b4 ~ 0 + time,
#'                start = c(b1 = 0, b2 = 1, b4 = 1000, b5 = 0),
#'                fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
#' ed(dmod, 50, newdata = data.frame(time = c("24h", "48h")))
"daphnids"

#' Earthworm toxicity test
#'
#' The number of earthworms remaining in a container contaminated with a toxic substance, i.e. not
#' migrating to a neighbouring uncontaminated container. At dose 0, about half of the earthworms are
#' expected in each container, so the upper asymptote is about 0.5 rather than 1.
#'
#' @format A data frame with 35 rows and 3 variables:
#' \describe{
#'   \item{dose}{dose of the toxic substance}
#'   \item{number}{number of earthworms remaining in the contaminated container}
#'   \item{total}{total number of earthworms in the containers}
#' }
#' @source Nina Cedergreen, Faculty of Life Sciences, University of Copenhagen; from the package
#'   drcData.
#' @examples
#' emod <- drmTMB(cbind(number, total - number) ~ dose, data = earthworms, model = "loglogistic",
#'                family = "binomial", start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
#' emod
"earthworms"

#' Germination of three crops
#'
#' Germination of mungbean, rice and wheat seeds at temperatures between 10 and 40 degrees Celsius,
#' recorded daily for up to 18 days: time-to-event data with several curves.
#'
#' For each species and temperature, 20 seeds were placed in a Petri dish in a dark climate cabinet
#' (10, 16, 22, 28, 34 or 40 degrees Celsius; not all temperatures for all species), and the germinated
#' seeds were counted and removed daily.
#'
#' @format A data frame with 192 rows and 5 variables:
#' \describe{
#'   \item{temp}{temperature (degrees Celsius)}
#'   \item{species}{a factor with levels \code{mungbean}, \code{rice} and \code{wheat}}
#'   \item{start}{left endpoint of the monitoring interval (days)}
#'   \item{end}{right endpoint of the monitoring interval (days)}
#'   \item{germinated}{number of seeds germinated in the interval}
#' }
#' @source From the package drcData.
#' @references Ritz, C., Pipper, C. B. and Streibig, J. C. (2013) Analysis of germination data from
#'   agricultural experiments. \emph{European Journal of Agronomy}, \bold{45}, 1--6.
#' @examples
#' wheat <- subset(germination, species == "wheat" & temp == 22)
#' wmod <- drmTMB(germinated ~ start + end, data = wheat, model = "loglogistic", family = "event",
#'                start = c(b5 = 0), fix = c(b5 = TRUE))
#' ed(wmod, 50)
"germination"

#' Hormesis in lettuce plants
#'
#' The shoot biomass of lettuce (\emph{Lactuca sativa}) plants after 21 days in a nutrient solution
#' with isobutyl alcohol. Low concentrations stimulate the growth (hormesis).
#'
#' @format A data frame with 14 rows and 2 variables:
#' \describe{
#'   \item{conc}{concentration of isobutyl alcohol (mg/l)}
#'   \item{weight}{shoot biomass (g)}
#' }
#' @source van Ewijk, P. H. and Hoekstra, J. A. (1993) Calculation of the EC50 and its confidence
#'   interval when subtoxic stimulus is present. \emph{Ecotoxicology and Environmental Safety},
#'   \bold{25}, 25--32; from the package drcData.
#' @references van Ewijk, P. H. and Hoekstra, J. A. (1994) Curvature measures and confidence intervals
#'   for the linear logistic model. \emph{Applied Statistics}, \bold{43}, 477--487.
#' @examples
#' hmod <- drmTMB(weight ~ conc, data = lettuce, model = "braincousens")
#' peak(hmod)
"lettuce"

#' Liver tumour incidence after exposure to TCDD (NTP 2006)
#'
#' Liver tumour incidence in female Harlan Sprague-Dawley rats exposed to
#' 2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD); the current study of the EFSA example on informative
#' priors, see \code{\link{TCDD}} for the historical study, and \code{\link{bmdefsa}} and
#' \code{\link{fitpert}} for the analysis with an informative prior from the historical study.
#'
#' @format A data frame with 6 rows and 3 variables:
#' \describe{
#'   \item{conc}{concentration of TCDD (ng/kg)}
#'   \item{total}{number of rats}
#'   \item{incidence}{number of rats with liver tumours}
#' }
#' @source National Toxicology Program (2006) NTP technical report on the toxicology and
#'   carcinogenesis studies of 2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD) (CAS No. 1746-01-6) in
#'   female Harlan Sprague-Dawley rats (gavage studies). \emph{National Toxicology Program Technical
#'   Report Series}, \bold{521}, 4--232; from the package drcData.
#' @examples
#' # no tumours in the controls: background risk b1 fixed at 0
#' lmod <- drmTMB(cbind(incidence, total - incidence) ~ conc, data = liver.tumor,
#'                model = "loglogistic", family = "binomial", start = c(b1 = 0, b2 = 1, b5 = 0),
#'                fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
#' bmd(lmod, bmr = 0.1, def = "excess")
"liver.tumor"

#' Effect of ferulic acid on the root growth of ryegrass
#'
#' The root length of perennial ryegrass (\emph{Lolium perenne} L.) for several concentrations of
#' ferulic acid; part of a study of the joint action of phenolic acids on root growth inhibition.
#'
#' @format A data frame with 24 rows and 2 variables:
#' \describe{
#'   \item{rootl}{root length (cm)}
#'   \item{conc}{concentration of ferulic acid (mM)}
#' }
#' @source Inderjit, Streibig, J. C. and Olofsdotter, M. (2002) Joint action of phenolic acid mixtures
#'   and its significance in allelopathy research. \emph{Physiologia Plantarum}, \bold{114}, 422--428;
#'   from the package drcData.
#' @examples
#' mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
#'               start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
#' ed(mod, c(10, 50, 90))
"ryegrass"

#' Potency of two herbicides in white mustard
#'
#' The dry matter of white mustard (\emph{Sinapis alba}) for doses of the herbicides glyphosate and
#' bentazone. The asymptotes of the two herbicides can be assumed equal, while their slopes and ED50s
#' differ.
#'
#' @format A data frame with 68 rows and 3 variables:
#' \describe{
#'   \item{Dose}{dose (g/ha)}
#'   \item{Herbicide}{a factor with levels \code{Bentazone} and \code{Glyphosate}}
#'   \item{DryMatter}{dry matter (g/pot)}
#' }
#' @source Christensen, M. G., Teicher, H. B. and Streibig, J. C. (2003) Linking fluorescence induction
#'   curve and biomass in herbicide screening. \emph{Pest Management Science}, \bold{59}, 1303--1310;
#'   from the package drcData.
#' @examples
#' smod <- drmTMB(DryMatter ~ Dose, data = S.alba, model = "loglogistic",
#'                fform = b3 + b4 ~ Herbicide, parlink = c(b4 = "log"),
#'                start = c(b5 = 0), fix = c(b5 = TRUE))
#' ed(smod, 50, newdata = data.frame(Herbicide = c("Bentazone", "Glyphosate")))
"S.alba"

#' Toxicity of four types of selenium
#'
#' The number of dead flies for concentrations of four types of selenium.
#'
#' @format A data frame with 25 rows and 4 variables:
#' \describe{
#'   \item{type}{type of selenium (1 to 4)}
#'   \item{conc}{total selenium concentration}
#'   \item{total}{number of flies}
#'   \item{dead}{number of dead flies}
#' }
#' @source Jeske, D. R., Xu, H. K., Blessinger, T., Jensen, P. and Trumble, J. (2009) Testing for the
#'   equality of EC50 values in the presence of unequal slopes with application to toxicity of selenium
#'   types. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{14},
#'   469--483; from the package drcData.
#' @examples
#' semod <- drmTMB(cbind(dead, total - dead) ~ conc, data = selenium, model = "loglogistic",
#'                 family = "binomial", fform = b3 + b4 ~ factor(type),
#'                 start = c(b2 = 1, b5 = 0), fix = c(b2 = TRUE, b5 = TRUE))
#' ed(semod, 50, newdata = data.frame(type = 1:4))
"selenium"

#' Liver tumour incidence after exposure to TCDD (Kociba et al. 1978)
#'
#' Liver tumour incidence in Sprague-Dawley rats exposed to 2,3,7,8-tetrachlorodibenzo-p-dioxin
#' (TCDD); the historical study of the EFSA example on informative priors, see
#' \code{\link{liver.tumor}} for the current study, and \code{\link{bmdefsa}} and \code{\link{fitpert}}
#' for the analysis.
#'
#' @format A data frame with 4 rows and 3 variables:
#' \describe{
#'   \item{conc}{concentration of TCDD (ng/kg)}
#'   \item{total}{number of rats}
#'   \item{incidence}{number of rats with liver tumours}
#' }
#' @source Kociba, R., Keyes, D., Beyer, J., Carreon, R., Wade, C., Dittenber, D., Kalnins, R.,
#'   Frauson, L., Park, C., Barnard, S., Hummel, R. and Humiston, C. (1978) Results of a two-year
#'   chronic toxicity and oncogenicity study of 2,3,7,8-tetrachlorodibenzo-p-dioxin in rats.
#'   \emph{Toxicology and Applied Pharmacology}, \bold{46}, 279--303; from the package drcData.
#' @examples
#' tmod <- drmTMB(cbind(incidence, total - incidence) ~ conc, data = TCDD, model = "loglogistic",
#'                family = "binomial", start = c(b2 = 1, b5 = 0), fix = c(b2 = TRUE, b5 = TRUE))
#' bmd(tmod, bmr = 0.1, def = "excess")
"TCDD"

#' Vinclozolin in an androgen receptor reporter gene assay
#'
#' Luminescence of Chinese hamster ovary cells exposed to concentrations of vinclozolin for 22 hours;
#' the experiment was repeated six times on different days.
#'
#' @format A data frame with 53 rows and 3 variables:
#' \describe{
#'   \item{exper}{experiment (day), a factor with 6 levels}
#'   \item{conc}{concentration of vinclozolin}
#'   \item{effect}{luminescence}
#' }
#' @source Nellemann, C., Dalgaard, M., Lam, H. R. and Vinggaard, A. M. (2003) The combined effects of
#'   vinclozolin and procymidone do not deviate from expected additivity in vitro and in vivo.
#'   \emph{Toxicological Sciences}, \bold{71}, 251--262; from the package drcData.
#' @examples
#' vmod <- drmTMB(effect ~ conc, data = vinclozolin, model = "loglogistic",
#'                rform = b1 + b2 ~ 1 | exper, start = c(b5 = 0), fix = c(b5 = TRUE),
#'                lower = c(0, 0, -Inf, 0, -Inf))  # lower bounds of b1-b5
#' # ED50 of the population curve, with profile likelihood interval
#' confint(ed(vmod, 50, random = FALSE))
"vinclozolin"
