% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetDNA.R
\name{metdnaR}
\alias{metdnaR}
\title{MetDNA}
\usage{
metdnaR(ms1.data.pos = "data.csv", ms1.data.neg = "data.csv",
  sample.info.pos = "sample.info.csv", sample.info.neg = "sample.info.csv",
  mz.tol = 25, rt.tol.for.ms1.ms2.match = 10, pos.path = ".",
  neg.path = ".", instrument = c("SciexTripleTOF", "AgilentQTOF",
  "OtherQTOF", "ThermoOrbitrap"), ms2.type = c("mgf", "msp", "mzXML"),
  polarity = c("positive", "negative", "both"), column = c("hilic", "rp"),
  ce = c("30", "10", "20", "35,15", "40", "50"), ms2.match.plot = FALSE,
  prefer.adduct = ifelse(polarity == "positive", "M+H", "M-H"),
  use.default.md = TRUE, threads = 3, max.isotope = 4, rt.tol1 = 3,
  rt.tol2 = 30, cor.tol = 0, int.tol = 500, dp.tol = 0.5,
  max.step = 3, score.cutoff = 0, remain = FALSE, remain.per = 0.5,
  seed.neighbor.match.plot = FALSE, candidate.num = 5, group,
  uni.test = c("t", "wilcox"), correct = TRUE, p.cutoff = 0.01,
  use.old.null = FALSE, species = c("hsa", "dme", "mmu", "rat", "bta",
  "gga", "dre", "cel", "sce", "ath", "smm", "pfa", "tbr", "eco", "ppu", "syf"),
  use.all.kegg.id = FALSE, only.mrn.annotation = FALSE, check.data = TRUE,
  ms2.annotation = TRUE, mrn.annotation = TRUE, dn.analysis = FALSE,
  pathway.enrichment = TRUE, parameter = NULL)
}
\arguments{
\item{ms1.data.pos}{The name of positive ms1 data.}

\item{ms1.data.neg}{The name of negative ms1 data.}

\item{sample.info.pos}{The name of positive sample.info.}

\item{sample.info.neg}{The name of negative sample.info.}

\item{mz.tol}{mz tol for ms1 and ms2 data matching.}

\item{rt.tol.for.ms1.ms2.match}{RT tol for ms1 and ms2 data matching.}

\item{pos.path}{If polarity is both, pos.path must be provided.}

\item{neg.path}{If polarity is both, neg.path must be provided.}

\item{instrument}{The instrument you used to acquire data. "AgilentQTOF",
"SciexTripleTOF", "OtherQTOF" and "ThermoOrbitrap" are supported.}

\item{ms2.type}{"mgf" or "msp".}

\item{polarity}{The polarity of mode.}

\item{column}{The column.}

\item{ce}{The collision energy.}

\item{ms2.match.plot}{Output MS/MS match plot or nor.}

\item{prefer.adduct}{Which adduct you want to use for RT prediction.
tparam use.default.md Use default molecular descriptors for RT prediction or not.}

\item{use.default.md}{Use default molecular descriptors to predict RTs.}

\item{threads}{The number of threads.}

\item{max.isotope}{The number of isotope peaks}

\item{rt.tol1}{The RT tolerance for isotope and adduct annotation. (second)
Default is 3.}

\item{rt.tol2}{The RT tolerance for metabolite annotation. (\%) Default
is 30\%.}

\item{cor.tol}{The correlation tolerance.}

\item{int.tol}{The intensity ratio tolerance of isotope annotation.}

\item{dp.tol}{The tolerance of dot product.}

\item{max.step}{The max number of reaction step.}

\item{score.cutoff}{Score cutoff of annotations.}

\item{remain}{Remain some seeds as validation or not.}

\item{remain.per}{The percentage of remained seeds. Default is 50\%.}

\item{seed.neighbor.match.plot}{Output seed neighbor match plot or not.}

\item{candidate.num}{How many candidates for peaks are outputted. Default is 3.}

\item{group}{The group you want to use.}

\item{uni.test}{The method of univariate test.}

\item{correct}{Correct p value or not.}

\item{p.cutoff}{The cutoff of p value. Default is 0.01.}

\item{use.old.null}{Use old null distribution of not in module analysis. Default is FALSE}

\item{species}{The species. "hsa" is Homo sapiens (human),
"dme" is Drosophlia melanogaster (fruit fly),
"mmu" is Mus musculus (mouse), "rat" is Rattus norvegicus (rat),
"bta" is Bos taurus (cow), "gga" is Gallus domesticus (chicken),
"dre" is Danio rerio (zebrafish), "cel" is Caenorharomyces elegans (nematode),
 "sce" is Saccharomyces cerevisaiae (yeast), "ath" is Arabidopsis thaliana (thale cress),
 "smm" is Schistosoma mansoni, "pfa" is Plasmodum falciparum 3D7,
 "tbr" is Trypanosoma brucei, "eco" is Escherichia coli K-12 MG1655,
"ppu" is Pseudomonas putida KT2440, "syf" is Synechococcus elongatus.}

\item{use.all.kegg.id}{Use all annotations from KEGG database.}

\item{only.mrn.annotation}{Only use MRN annotations for dysregulated network analysis or not.}

\item{check.data}{Check data or not.}

\item{ms2.annotation}{MS/MS annotation or not.}

\item{mrn.annotation}{MRN based annotation or not.}

\item{dn.analysis}{Dysregulated network analysis or not.}

\item{pathway.enrichment}{Pathway enrichment analysis or not.}

\item{parameter}{The parameter path. Default is NULL.}
}
\description{
Metabolite annotation and dysregulated network analysis.
}
\author{
Xiaotao Shen
\email{shenxt@sioc.ac.cn}
}
