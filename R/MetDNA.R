#' @title MetDNA
#' @description Metabolite annotation and dysregulated network analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param ms1.data.pos The name of positive ms1 data.
#' @param ms1.data.neg The name of negative ms1 data.
#' @param sample.info.pos The name of positive sample.info.
#' @param sample.info.neg The name of negative sample.info.
#' @param mz.tol mz tol for ms1 and ms2 data matching.
#' @param rt.tol.for.ms1.ms2.match RT tol for ms1 and ms2 data matching.
#' @param pos.path If polarity is both, pos.path must be provided.
#' @param neg.path If polarity is both, neg.path must be provided.
#' @param instrument The instrument you used to acquire data. "AgilentQTOF",
#' "SciexTripleTOF", "OtherQTOF" and "ThermoOrbitrap" are supported.
#' @param ms2.type "mgf" or "msp".
#' @param polarity The polarity of mode.
#' @param column The column.
#' @param ce The collision energy.
#' @param ms2.match.plot Output MS/MS match plot or nor.
#' @param prefer.adduct Which adduct you want to use for RT prediction.
#' tparam use.default.md Use default molecular descriptors for RT prediction or not.
#' @param use.default.md Use default molecular descriptors to predict RTs.
#' @param threads The number of threads.
#' @param max.isotope The number of isotope peaks
#' @param rt.tol1 The RT tolerance for isotope and adduct annotation. (second)
#' Default is 3.
#' @param rt.tol2 The RT tolerance for metabolite annotation. (\%) Default
#' is 30\%.
#' @param cor.tol The correlation tolerance.
#' @param int.tol The intensity ratio tolerance of isotope annotation.
#' @param dp.tol The tolerance of dot product.
#' @param max.step The max number of reaction step.
#' @param score.cutoff Score cutoff of annotations.
#' @param remain Remain some seeds as validation or not.
#' @param remain.per The percentage of remained seeds. Default is 50\%.
#' @param seed.neighbor.match.plot Output seed neighbor match plot or not.
#' @param candidate.num How many candidates for peaks are outputted. Default is 3.
#' @param group The group you want to use.
#' @param uni.test The method of univariate test.
#' @param correct Correct p value or not.
#' @param p.cutoff The cutoff of p value. Default is 0.01.
#' @param use.old.null Use old null distribution of not in module analysis. Default is FALSE
#' @param species The species. "hsa" is Homo sapiens (human),
#' "dme" is Drosophlia melanogaster (fruit fly),
#' "mmu" is Mus musculus (mouse), "rat" is Rattus norvegicus (rat),
#' "bta" is Bos taurus (cow), "gga" is Gallus domesticus (chicken),
#' "dre" is Danio rerio (zebrafish), "cel" is Caenorharomyces elegans (nematode),
#'  "sce" is Saccharomyces cerevisaiae (yeast), "ath" is Arabidopsis thaliana (thale cress),
#'  "smm" is Schistosoma mansoni, "pfa" is Plasmodum falciparum 3D7,
#'  "tbr" is Trypanosoma brucei, "eco" is Escherichia coli K-12 MG1655,
#' "ppu" is Pseudomonas putida KT2440, "syf" is Synechococcus elongatus.
#' @param use.all.kegg.id Use all annotations from KEGG database.
#' @param only.mrn.annotation Only use MRN annotations for dysregulated network analysis or not.
#' @param ms2.annotation MS/MS annotation or not.
#' @param mrn.annotation MRN based annotation or not.
#' @param dn.analysis Dysregulated network analysis or not.
#' @param check.data Check data or not.
#' @param pathway.enrichment Pathway enrichment analysis or not.
#' @param parameter The parameter path. Default is NULL.
#' @export

#
# metdnaR(ms1.data.pos = "data.csv",
#           ms1.data.neg = "data.csv",
#           sample.info.pos = "sample.info.csv",
#           sample.info.neg = "sample.info.csv",
#           pos.path = "V:/workreport/Shen Xiaotao/MetDNA demo data/metDNA/fly20190921/POS",
#           neg.path = "V:/workreport/Shen Xiaotao/MetDNA demo data/metDNA/fly20190921/NEG",
#           polarity = "both",
#           column = "rp",
#           ce = "30",
#           use.default.md = TRUE,
#           group = c("W03", "W30"),
#           uni.test = "t",
#           correct = FALSE,
#           p.cutoff = 0.05,
#           species = "mmu",
#           dn.analysis = FALSE,
#           instrument = "OtherQTOF",
#           pathway.enrichment = TRUE,
#        parameter = NULL)


setGeneric(name = "metdnaR",
           def = function(ms1.data.pos = "data.csv",
                          ms1.data.neg = "data.csv",
                          sample.info.pos = "sample.info.csv",
                          sample.info.neg = "sample.info.csv",
                          mz.tol = 25,
                          rt.tol.for.ms1.ms2.match = 10,#second
                          pos.path = ".",
                          neg.path = ".",
                          instrument = c("SciexTripleTOF", "AgilentQTOF",
                                         "OtherQTOF", "ThermoOrbitrap"),
                          ms2.type = c("mgf", "msp", "mzXML"),
                          polarity = c("positive", "negative", "both"),
                          column = c("hilic", "rp"),
                          ce = c("30", "10", "20", "35,15", "40", "50"),
                          ms2.match.plot = FALSE,
                          ###
                          prefer.adduct = ifelse(polarity == "positive", "M+H", "M-H"),
                          use.default.md = TRUE,
                          threads = 3,
                          max.isotope = 4,
                          rt.tol1 = 3,#second
                          rt.tol2 = 30,#%
                          cor.tol = 0,
                          int.tol = 500,
                          dp.tol = 0.5,
                          max.step = 3,
                          score.cutoff = 0,
                          remain = FALSE,
                          remain.per = 0.5,
                          seed.neighbor.match.plot = FALSE,
                          candidate.num = 5,
                          ###
                          group,
                          uni.test = c("t", "wilcox"),
                          correct = TRUE,
                          p.cutoff = 0.01,
                          use.old.null = FALSE,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE,
                          check.data = TRUE,
                          ms2.annotation = TRUE,
                          mrn.annotation = TRUE,
                          dn.analysis = FALSE,
                          pathway.enrichment = TRUE,
                          parameter = NULL
           ){


             cat("\n")
             packageStartupMessage("MetDNA version 1.1.2.\n")
             #--------------------------------------------------------------------
             #parameter decision
             ms2.type <- match.arg(ms2.type)
             polarity <- match.arg(polarity)
             column <- match.arg(column)
             instrument <- match.arg(instrument)
             ce <- match.arg(ce)
             uni.test <- match.arg(uni.test)
             species <- match.arg(species)

             if(!is.null(parameter)){
               cat("Use the old parameter from the MetDNA.parameters.csv\n")
               parameter <- try({read.csv(parameter, stringsAsFactors = FALSE)}, silent = TRUE)

               if(class(parameter) == "try-error"){
                 cat("The parameter you provide is wrong, please check it.\n")
                 stop(parameter[[1]])
               }

               mz.tol <- as.numeric(parameter$value[parameter$parameter=="mz.tol"])
               rt.tol.for.ms1.ms2.match <- as.numeric(parameter$value[parameter$parameter=="rt.tol.for.ms1.ms2.match"])
               ms2.type <- parameter$value[parameter$parameter=="ms2.type"]
               polarity <- parameter$value[parameter$parameter=="polarity"]
               column <- parameter$value[parameter$parameter=="column"]
               ce <- parameter$value[parameter$parameter=="ce"]
               prefer.adduct <- parameter$value[parameter$parameter=="prefer.adduct"]
               threads <- as.numeric(parameter$value[parameter$parameter=="threads"])
               rt.tol1 <- as.numeric(parameter$value[parameter$parameter=="rt.tol1"])
               rt.tol2 <- as.numeric(parameter$value[parameter$parameter=="rt.tol2"])
               dp.tol <- as.numeric(parameter$value[parameter$parameter=="dp.tol"])
               group <- parameter$value[parameter$parameter=="group"]
               group <- strsplit(x = group, split = ",")[[1]]
               uni.test <- parameter$value[parameter$parameter=="uni.test"]
               correct <- as.logical(parameter$value[parameter$parameter=="correct"])
               p.cutoff <- as.numeric(parameter$value[parameter$parameter=="p.cutoff"])
               species <- parameter$value[parameter$parameter=="species"]
               instrument <- parameter$value[parameter$parameter=="instrument"]
             }

             ##path decision
             if(polarity == "both"){
               path <- paste(strsplit(pos.path, split = "/")[[1]][-length(strsplit(pos.path, split = "/")[[1]])], collapse = "/")
               path <- file.path(path, "POS and NEG")
             }else{
               path <- ifelse(polarity == "positive", pos.path, neg.path)
             }
             dir.create(path)

             ##creat run.log
             cat("MetDNA run log\n", file = file.path(path, "run.log.txt"))
             cat("MetDNA version: 1.1.2\n", file = file.path(path, "run.log.txt"), append = TRUE)
             cat(as.character(Sys.time()), file = file.path(path, "run.log.txt"), append = TRUE)
             cat("\n", file = file.path(path, "run.log.txt"), append = TRUE)

             ###check the sample.info.pos and sample.info.neg is same or not.
             old.sample.info.pos.name = sample.info.pos
             old.sample.info.neg.name = sample.info.neg
             if(polarity == "both"){
               temp.error <- errorDisplay(
                 sample.info.pos <- read.csv(file.path(pos.path, old.sample.info.pos.name), stringsAsFactors = FALSE),
                 error.info = paste("Error: There is no", old.sample.info.pos.name, "in your pos.path. Please check it.\n"))

               if(class(temp.error) == "try-error") stop("Error")

               temp.error <- errorDisplay(
                 sample.info.neg <- read.csv(file.path(pos.path, old.sample.info.neg.name), stringsAsFactors = FALSE),
                 error.info = paste("Error: There is no", old.sample.info.pos.name, "in your neg.path. Please check it.\n"))
               if(class(temp.error) == "try-error") stop("Error")

               if(!all(sample.info.pos == sample.info.neg)){
                 cat("Error: The sample.info in POS and NEG should be same.\n", file = file.path(path, "run.log.txt"), append = TRUE)
                 stop("The sample.info in POS and NEG should be same.\n")
               }
             }

             ##save parameters
             MetDNA.parameters <- c(mz.tol,
                                    rt.tol.for.ms1.ms2.match,
                                    ms2.type,
                                    polarity,
                                    column,
                                    ce,
                                    prefer.adduct,
                                    threads,
                                    rt.tol1,
                                    rt.tol2,
                                    dp.tol,
                                    paste(group, collapse = ","),
                                    uni.test,
                                    correct,
                                    p.cutoff,
                                    species,
                                    instrument)
             MetDNA.parameters <- data.frame(c("mz.tol",
                                               "rt.tol.for.ms1.ms2.match",
                                               "ms2.type",
                                               "polarity",
                                               "column",
                                               "ce",
                                               "prefer.adduct",
                                               "threads",
                                               "rt.tol1",
                                               "rt.tol2",
                                               "dp.tol",
                                               "group",
                                               "uni.test",
                                               "correct",
                                               "p.cutoff",
                                               "species",
                                               "instrument"),
                                             MetDNA.parameters,
                                             stringsAsFactors = FALSE)

             MetDNA.parameters <- rbind(c("Version", "1.0.1"), MetDNA.parameters)
             colnames(MetDNA.parameters) <- c('parameter', 'value')
             write.csv(MetDNA.parameters, file.path(path, "MetDNA.parameters.csv"))

             if(polarity == "both"){
             MetDNA.parameters.pos <- MetDNA.parameters
             MetDNA.parameters.pos[4,2] <- "positive"
             write.csv(MetDNA.parameters.pos, file.path(pos.path, "MetDNA.parameters.csv"))

             MetDNA.parameters.neg <- MetDNA.parameters
             MetDNA.parameters.neg[4,2] <- "negative"
             write.csv(MetDNA.parameters.neg, file.path(neg.path, "MetDNA.parameters.csv"))
             }


             ###check data
             if(polarity == "positive" | polarity == "both"){
               cat("========================================================\n")
               cat("MetDNA analysis for positive mode data.\n")
               cat("========================================================\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA analysis for positive mode data.\n", file = file.path(path, "run.log.txt"), append = TRUE)

               if(polarity == "positive") {
                 analysis.report <- TRUE
                 }else{
                   analysis.report <- TRUE
                 }


             #   metdnaForOneMode(ms1.data = ms1.data.pos,
             #                    sample.info = old.sample.info.pos.name,
             #                    mz.tol = mz.tol,
             #                    rt.tol.for.ms1.ms2.match = rt.tol.for.ms1.ms2.match,#second
             #                    path = pos.path,
             #                    log.path = path,
             #                    instrument = instrument,
             #                    ms2.type = ms2.type,
             #                    polarity = "positive",
             #                    column = column,
             #                    ce = ce,
             #                    ms2.match.plot = ms2.match.plot,
             #                    candidate.num = candidate.num,
             #                    ###
             #                    prefer.adduct = prefer.adduct,
             #                    use.default.md = use.default.md,
             #                    threads = threads,
             #                    max.isotope = max.isotope,
             #                    rt.tol1 = rt.tol1,#second
             #                    rt.tol2 = rt.tol2,#%
             #                    cor.tol = cor.tol,
             #                    int.tol = int.tol,
             #                    dp.tol = dp.tol,
             #                    max.step = max.step,
             #                    score.cutoff = score.cutoff,
             #                    remain = remain,
             #                    remain.per = remain.per,
             #                    seed.neighbor.match.plot = seed.neighbor.match.plot,
             #                    ###
             #                    group = group,
             #                    uni.test = uni.test,
             #                    correct = correct,
             #                    p.cutoff = p.cutoff,
             #                    use.old.null = use.old.null,
             #                    species = species,
             #                    use.all.kegg.id = use.all.kegg.id,
             #                    only.mrn.annotation = only.mrn.annotation,
             #                    check.data = check.data,
             #                    ms2.annotation = ms2.annotation,
             #                    mrn.annotation = mrn.annotation,
             #                    dn.analysis = ifelse(polarity == "both", FALSE, dn.analysis),
             #                    pathway.enrichment = ifelse(polarity == "both", FALSE, pathway.enrichment),
             #                    analysis.report = analysis.report)
             }


             if(polarity == "negative" | polarity == "both"){
               cat("========================================================\n")
               cat("========================================================\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA analysis for negative mode data.\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA analysis for negative mode data.\n")

               if(polarity == "positive") {
                 analysis.report <- TRUE
               }else{
                 analysis.report <- TRUE
               }

               # metdnaForOneMode(ms1.data = ms1.data.neg,
               #                  sample.info = old.sample.info.neg.name,
               #                  mz.tol = mz.tol,
               #                  rt.tol.for.ms1.ms2.match = rt.tol.for.ms1.ms2.match,#second
               #                  path = neg.path,
               #                  log.path = path,
               #                  ms2.type = ms2.type,
               #                  polarity = "negative",
               #                  column = column,
               #                  ce = ce,
               #                  ms2.match.plot = ms2.match.plot,
               #                  candidate.num = candidate.num,
               #                  ###
               #                  prefer.adduct = prefer.adduct,
               #                  use.default.md = use.default.md,
               #                  threads = threads,
               #                  max.isotope = max.isotope,
               #                  rt.tol1 = rt.tol1,#second
               #                  rt.tol2 = rt.tol2,#%
               #                  cor.tol = cor.tol,
               #                  int.tol = int.tol,
               #                  dp.tol = dp.tol,
               #                  max.step = max.step,
               #                  score.cutoff = score.cutoff,
               #                  remain = remain,
               #                  remain.per = remain.per,
               #                  seed.neighbor.match.plot = seed.neighbor.match.plot,
               #                  ###
               #                  group = group,
               #                  uni.test = uni.test,
               #                  correct = correct,
               #                  p.cutoff = p.cutoff,
               #                  use.old.null = use.old.null,
               #                  species = species,
               #                  use.all.kegg.id = use.all.kegg.id,
               #                  only.mrn.annotation = only.mrn.annotation,
               #                  check.data = check.data,
               #                  ms2.annotation = ms2.annotation,
               #                  mrn.annotation = mrn.annotation,
               #                  dn.analysis = ifelse(polarity == "both", FALSE, dn.analysis),
               #                  pathway.enrichment = ifelse(polarity == "both", FALSE, pathway.enrichment),
               #                  analysis.report = analysis.report
               # )
             }



             if(polarity == "both"){
               cat("\n")
               cat("========================================================\n")
               cat("========================================================\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA analysis for positive and negative mode data.\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA analysis for positive and negative mode data.\n")

               # metdnaFor2Mode(ms1.data.pos = ms1.data.pos,
               #                ms1.data.neg = ms1.data.neg,
               #                sample.info.pos = old.sample.info.pos.name,
               #                sample.info.neg = old.sample.info.neg.name,
               #                mz.tol = mz.tol,
               #                rt.tol.for.ms1.ms2.match = rt.tol.for.ms1.ms2.match,#second
               #                path = path,
               #                log.path = path,
               #                pos.path = pos.path,
               #                neg.path = neg.path,
               #                ms2.type = ms2.type,
               #                column = column,
               #                ce = ce,
               #                candidate.num = candidate.num,
               #                ###
               #                use.default.md = use.default.md,
               #                threads = threads,
               #                max.isotope = max.isotope,
               #                rt.tol1 = rt.tol1,#second
               #                rt.tol2 = rt.tol2,#%
               #                cor.tol = cor.tol,
               #                int.tol = int.tol,
               #                dp.tol = dp.tol,
               #                max.step = max.step,
               #                score.cutoff = score.cutoff,
               #                # candidate.num = candidate.num,
               #                ###
               #                group = group,
               #                uni.test = uni.test,
               #                correct = correct,
               #                p.cutoff = p.cutoff,
               #                use.old.null = use.old.null,
               #                species = species,
               #                use.all.kegg.id = use.all.kegg.id,
               #                only.mrn.annotation = only.mrn.annotation,
               #                dn.analysis = dn.analysis,
               #                pathway.enrichment = pathway.enrichment
               # )


               #
               ###delete the Dysregulated_network_analysis_result and Pathway_enrichment_analysis_result
               ##in the POS and NEG folder
               # unlink(x = file.path(pos.path, "Dysregulated_network_analysis_result"), recursive = TRUE)
               # unlink(x = file.path(pos.path, "Analysis_report"), recursive = TRUE)
               # unlink(x = file.path(pos.path, "Pathway_enrichment_analysis_result"), recursive = TRUE)
               # unlink(x = file.path(neg.path, "Analysis_report"), recursive = TRUE)
               # unlink(x = file.path(neg.path, "Dysregulated_network_analysis_result"), recursive = TRUE)
               # unlink(x = file.path(neg.path, "Pathway_enrichment_analysis_result"), recursive = TRUE)
               # unlink(x = file.path(path, "Dysregulated_network_analysis_result"), recursive = TRUE)
               # changeFile(path = file.path(path, "Pathway_enrichment_analysis_result"))
               # unlink(x = file.path(pos.path, "MetDNA.parameters.csv"), recursive = TRUE)
               # unlink(x = file.path(neg.path, "MetDNA.parameters.csv"), recursive = TRUE)
               #
               # unlink(x = file.path(pos.path, "MS2_match_result"), recursive = TRUE)
               # unlink(x = file.path(neg.path, "MS2_match_result"), recursive = TRUE)


               ##add P-value to Quantitative.pathway.metabolite.result
               # temp.error <- try(temp.data <- readr::read_csv(file.path(path, "Pathway_enrichment_analysis_result",
               #                                        "Quantitative.pathway.metabolite.result.csv")), silent = TRUE)
               #
               # if(class(temp.error) != "try-error"){
               #   temp.data <- as.data.frame(temp.data)
               # load(file.path(path, "Pathway_enrichment_analysis_result", "intermediate_data", "p.value.pos"))
               #   load(file.path(path, "Pathway_enrichment_analysis_result", "intermediate_data", "p.value.neg"))
               #
               #   p.value <- c(p.value.pos, p.value.neg)
               #   p <- NULL
               #   p <- p.value[match(temp.data$peak.name, names(p.value))]
               #
               #   temp.data <- data.frame(temp.data[,1:4], p,
               #                           temp.data[,-c(1:4)], stringsAsFactors = FALSE)
               #
               #   write.csv(temp.data, file.path(path, "Pathway_enrichment_analysis_result",
               #                                  "Quantitative.pathway.metabolite.result.csv"), row.names = FALSE)
               #
               # }

             }




             #remove some files
             # if(polarity == "positive"){
             #    unlink(x = file.path(pos.path, "Dysregulated_network_analysis_result"), recursive = TRUE)
             #    unlink(x = file.path(pos.path, "MS2_match_result"), recursive = TRUE)
             #    changeFile(path = file.path(pos.path, "Pathway_enrichment_analysis_result"))
             #
             #    temp.error <- try(temp.data <- readr::read_csv(file.path(pos.path, "Pathway_enrichment_analysis_result",
             #                                                             "Quantitative.pathway.metabolite.result.csv")), silent = TRUE)
             #
             #    if(class(temp.error) != "try-error"){
             #      temp.data <- as.data.frame(temp.data)
             #      load(file.path(pos.path, "Pathway_enrichment_analysis_result", "intermediate_data", "p.value"))
             #      p <- NULL
             #      p <- p.value[match(temp.data$peak.name, names(p.value))]
             #
             #      temp.data <- data.frame(temp.data[,1:4], p, temp.data[,-c(1:4)],
             #                              stringsAsFactors = FALSE)
             #
             #      write.csv(temp.data, file.path(pos.path, "Pathway_enrichment_analysis_result",
             #                                     "Quantitative.pathway.metabolite.result.csv"), row.names = FALSE)
             #
             #    }
             #  }

             # if(polarity == "negative"){
             #   unlink(x = file.path(neg.path, "Dysregulated_network_analysis_result"), recursive = TRUE)
             #   unlink(x = file.path(neg.path, "MS2_match_result"), recursive = TRUE)
             #   changeFile(path = file.path(neg.path, "Pathway_enrichment_analysis_result"))
             #
             #   temp.error <- try(temp.data <- readr::read_csv(file.path(neg.path, "Pathway_enrichment_analysis_result",
             #                                                            "Quantitative.pathway.metabolite.result.csv")), silent = TRUE)
             #
             #   if(class(temp.error) != "try-error"){
             #     temp.data <- as.data.frame(temp.data)
             #     load(file.path(neg.path, "Pathway_enrichment_analysis_result", "intermediate_data", "p.value"))
             #     p <- NULL
             #     p <- p.value[match(temp.data$peak.name, names(p.value))]
             #
             #     temp.data <- data.frame(temp.data[,1:4], p, temp.data[,-c(1:4)],
             #                             stringsAsFactors = FALSE)
             #
             #     write.csv(temp.data, file.path(neg.path, "Pathway_enrichment_analysis_result",
             #                                    "Quantitative.pathway.metabolite.result.csv"), row.names = FALSE)
             #
             #   }
             # }


             cat("\n")
             cat("========================================================\n", file = file.path(path, "run.log.txt"), append = TRUE)
             cat("##MetDNA is done.\n", file = file.path(path, "run.log.txt"), append = TRUE)
             cat("MetDNA is done.\n")
           })






.onAttach <- function(libname, pkgname){
  packageStartupMessage("
More information can be found in http://metdna.zhulab.cn/index.
If you have any questions, please send email to shenxt1990@163.com.
Authors: Xiaotao Shen and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).
Maintainer: Xiaotao Shen.\n2018-11-22
News:
Version 1.1.2
--------------
o Fixed a small bug, when there is no marker with ID for pathway enrichment.")
}


packageStartupMessage("
More information can be found in http://metdna.zhulab.cn/index.
If you have any questions, please send email to shenxt1990@163.com.
Authors: Xiaotao Shen and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).
Maintainer: Xiaotao Shen.\n2018-11-22
Version 1.1.2
--------------
o Fixed a small bug, when there is no marker with ID for pathway enrichment.")










