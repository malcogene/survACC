options(stringsAsFactors = FALSE) 
options(mySeedV=76543217)

loadUrl <-
  function(url,
           downloadPath = NA,
           sep = c("RData", " ", "," , "\t", ";", "xls", "gsheet"),
           onedrive = T,
           ...) {
    cat(
      'onedrive: copy link\n googlesheet: share-> Anyone with the link\n sep: "RData", ..."xls", "gsheet"\n'
    )
    if (!is.na(downloadPath))  {
      tmpFile <- downloadPath
      
    } else {
      tmpFile <-
        file.path(getwd(), paste0(substr(Sys.time(), 1, 10), '.rda'))
    }
    if (onedrive) {
      url2 <- gsub("e=.*", "download=1", url)
    } else {
      url2 <- url
    }
    
    download.file(url2, tmpFile, mode = "wb") 
    sep <- match.arg(sep)
    if (sep == "RData") {
      print(tmpFile)
      tmpFile <-  gsub("\\\\", "/", tmpFile)
      justLoaded <- try(load(tmpFile), silent = T)
      ;
      try(assign(justLoaded, eval(as.symbol(justLoaded)), .GlobalEnv), silent = T)
      ;
      if (class(justLoaded) == "try-error") {
        justLoaded <-
          try(read.delim(tmpFile, ...), silent = T)
        ; message("Need 'sep' argument, is it txt file?")
      }
    } else if (sep == "xls") {
      is.installed('readxl')
      justLoaded <- try(readxl::read_excel(tmpFile, ...), silent = T)
      
    } else if (sep == "gsheet") {
      is.installed("gsheet")
      cat('gsheet should be public, click share-> Anyone with the link\n')
      justLoaded <- gsheet::gsheet2tbl(url, ...)
    } else {
      justLoaded <- try(read.delim(tmpFile, sep = sep, ...), silent = T)
    }
    justLoaded
  }




mapg <-
  function(keys,
           probeSeq = F,
           probePlatform = c("hgu133a", "hgu133plus2", "GPL10558")) {
    require(org.Hs.eg.db) 
    cat.box(
      "mapG pipeline" ,
      " 1. EntrezID \n 2. HGNC official symbols were mapped to entrez ID (cf. in case of asNA => show alias to HGNC)"
    )
    
    if (length(keys) == 0)
      return(NA)
    keys <-  gsub("///.*", "", keys)
    temp <-
      tryCatch(
        as.numeric(keys),
        warning = function(x)
          print("Assumes that terms are all HGNC approved symbols")
      )
    probePlatform = match.arg(probePlatform)
    
    if (probeSeq) {
      if (any(probePlatform == "hgu133a")) {
        if (!is.numeric(temp)) {
          require(hgu133a.db)
          columns(hgu133a.db)
          a = mapIds(hgu133a.db, keys, c("PROBEID"), "SYMBOL", multiVals =
                       "list")
          print(a)
        } else {
          a = mapIds(hgu133a.db, keys, c("PROBEID"), "ENTREZID", multiVals = "list")
          print(a)
        }
        require(hgu133aprobe) 
        data(hgu133aprobe)
        return(hgu133aprobe[which(hgu133aprobe[, "Probe.Set.Name"] == a[[1]][1]), ])
        
      }
      if (any (probePlatform == "hgu133aplus2")) {
        if (!is.numeric(temp)) {
          require(hgu133plus2.db)
          a = mapIds(hgu133plus2.db,
                     keys,
                     c("PROBEID"),
                     "SYMBOL",
                     multiVals = "list")
          print(a)
        } else {
          a = mapIds(hgu133a.db, keys, c("PROBEID"), "ENTREZID", multiVals = "list")
          print(a)
        }
        require(hgu133aprobe)
        data(hgu133aprobe)
        return(hgu133aprobe[which(hgu133aprobe[, "Probe.Set.Name"] == a[[1]][1]), ])
      }
      
      if (any (probePlatform == "GPL10558")) {
        load("~/rawTemp/GPL10558-50081.RData")
        colnames(illu.seq)
        print(keys)
        return(illu.seq[which(illu.seq[, "Symbol"] == keys),])
      }
    }
    if (is.numeric(temp)) {
      hgnc <-
        mapIds(
          org.Hs.eg.db,
          keys = keys,
          column = "SYMBOL",
          keytype = "ENTREZID",
          multiVals = "asNA"
        )
      duplicated.inx <<- which(duplicated(hgnc))
      cat.box(
        "Handling multiVals = asNA",
        paste0(
          paste(keys[duplicated.inx], collapse = ","),
          "\n\nIndex of above ID wewe now assigned as global variable 'duplicated.inx'\n\n"
        )
      )
      return(hgnc) 
    } else {
      entrez <-
        mapIds(
          org.Hs.eg.db,
          keys = keys,
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "asNA"
        )
      
      multi.hgnc <-
        mapIds(
          org.Hs.eg.db,
          keys = keys,
          column = "SYMBOL",
          keytype = "ALIAS",
          multiVals = "list"
        )
      a = as.list(keys)
      b = multi.hgnc
      c = c()
      for (i in 1:length(a)) {
        c[i] <- isTRUE(a[[i]] == b[[i]])
      }
      diff.aliasToHGNC <- multi.hgnc[which(c == F)]
      attributes(entrez)$keytype <- "SYMBOL"
      attributes(entrez)$multiVals <- "asNA"
      attributes(entrez)$aliasToMulti.hgnc <- multi.hgnc
      attributes(entrez)$diff.aliasToHGNC <- diff.aliasToHGNC
      
      if (any(is.na(entrez))) {
        first.hgnc <-
          mapIds(
            org.Hs.eg.db,
            keys = keys,
            column = "SYMBOL",
            keytype = "ALIAS"
          )
        aliasToFirstHGNCToEntrez <-
          mapIds(
            org.Hs.eg.db,
            keys = first.hgnc,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "asNA"
          )
        entrez[which(is.na(entrez))] <-
          aliasToFirstHGNCToEntrez[which(is.na(entrez))]
        return(entrez)
      }   else  {
        return(entrez)
      }
    }
  }




factorizeDf <- function(df, lengthOfUniqueCut = 10) {
  if (!is.null(df)) {
    df.tmp <-
      data.frame(lapply(df, function(x) {
        if (length(unique(x)) < lengthOfUniqueCut) {
          as.factor(x)
        } else {
          x
        }
      }))
    colnames(df.tmp) <-
      colnames(df)
    rownames(df.tmp) <- rownames(df)
    df.tmp
  }
}



.w.MetaDE.merge <-
  function (x, MVperc = 0)
  {
    if (length(grep("///", x[[1]]$symbol)) != 0) {
      inxMultipleEntrezMapped <- grep("///", x[[1]]$symbol)
      x[[1]]$x <- x[[1]]$x[-inxMultipleEntrezMapped,]
      x[[1]]$symbol <-
        x[[1]]$symbol[-inxMultipleEntrezMapped]
    } 
    
    x.symbol <- lapply(x, function(y)
      toupper(rownames(y[[1]])))
    id.count <- table(unlist(x.symbol))
    n <- length(x.symbol)
    common.id <- names(which(id.count >= (1 - MVperc) * n))
    match.id <- function(study) {
      exprs <- study
      rownames(exprs) <- toupper(rownames(exprs))  
      id <- rownames(exprs)
      diff.id <- setdiff(common.id, id)
      n.sample <- ncol(exprs)
      margin.na <-
        matrix(NA, ncol = n.sample, nrow = length(diff.id))
      colnames(margin.na) <- colnames(exprs)
      rownames(margin.na) <- diff.id
      exprs <- rbind(exprs, margin.na)
      index <- match(common.id, rownames(exprs))
      exprs2 <- exprs[index,]
      return(exprs2)
    }
    K <- length(x)
    for (i in 1:K) {
      x[[i]][[1]] <- match.id(x[[i]][[1]])
    }
    return(x)
  }




sd.filter <-
  function (x,
            min.std = 0.25,
            zero.filter = 10,
            delQuantile = c(0, 0),
            plot = T) {
    if (delQuantile[1] == 0 && delQuantile[2] == 0) {
      predim = dim(x)[1]
      if (zero.filter != 0) {
        x = x[rowSums(x == 0) < zero.filter,]
        
        print(
          zero <-
            sprintf(
              "Zero filter: %s x[rowSums(==0)< %s, ] were filtered out",
              (predim - dim(x)[1]),
              zero.filter
            )
        )
      }
      
      tmp <- logical(dim(x)[1])
      if (is.numeric(min.std)) {
        for (i in 1:dim(x)[1]) {
          tmp[i] <- sd(x[i,], na.rm = T)
        }
        inx <- tmp > min.std
        print(tmp)
        print(inx)
        inx[is.na(inx)] <- TRUE
        cat(sprintf("\n sd(gene), %s) were filtered\n", min.std))
        cat(paste(sum(!inx), "genes excluded.\n"))
      }
      if (plot) {
        plot(sort(tmp), xlab = "Ordered genes", ylab = "Sd")
        abline(h = min.std, col = "red")
      }
      filtered <- x[inx,]
    } else {
      rank = order(apply(x, 1, mean, na.rm = T), decreasing = T)
      mean_mv <- x[which(rank > quantile(rank, delQuantile[1])), ]
      sd.rank = order(tmp2 <-
                        apply(x[rownames(mean_mv), ], 1, sd, na.rm = T), decreasing = T)
      sd_mv <-
        x[rownames(mean_mv), ][which(sd.rank > quantile(sd.rank, delQuantile[2])), ]
      cat(
        sprintf(
          "\ %s quantile of mean(gene),  %s quantile of sd(gene)) were filtered\n",
          delQuantile[1],
          delQuantile[2]
        )
      )
      filtered <- sd_mv
      if (plot) {
        plot(sort(tmp2), xlab = "Ordered genes", ylab = "Sd")
        abline(h = delQuantile[2], col = "red")
      }
    }
    print(zero)
    filtered
  }




summ.gene <-
  function(x,
           filter = c(0.25, 0.25),
           pl = c("both", "summ.probe", "gmatch"),
           o = c("IQR", "average"),
           forCentromics = F) {
    o <- match.arg(o)
    pl <- match.arg(pl)
    
    if (forCentromics) {
      filter = NULL
      pl = "both"
    }
    
    if (pl == "summ.probe") {
      ii <- MetaDE::MetaDE.match(x, o)
      for (i in 1:length(x)) {
        colnames(ii[[i]][[1]]) <- colnames(x[[i]][[1]])
      }
      
      ii = check.tar.list(ii)
      return(ii)
    }
    
    if (pl == "gmatch") {
      iii <- .w.MetaDE.merge(x)
      if (!is.null(filter))
        iii <- wrapper.MetaDE.filter(iii, c(filter[1], filter[2]))
      
      cat(" MetaDE.merge were performed")
      cat("\n")
      cat(
        sprintf(
          "=<quantile(mean(gene), %s), =<quantile(sd(gene), %s) were filtered",
          filter[1],
          filter[2]
        )
      )
      
      iii = check.tar.list(iii)
      return(iii)
    }
    if (pl == "both") {
      ii <-
        tryCatch(
          MetaDE::MetaDE.match(x, o),
          error = function(e) {
            print("Are you sure you need IQR summ.probe? ")
            x
          }
        )
      
      iii <- .w.MetaDE.merge(ii)
      print(names(iii))
      for (i in 1:length(ii)) {
        colnames(iii[[i]][[1]]) <- colnames(x[[i]][[1]])
      }
      if (!is.null(filter))
        iii <-  wrapper.MetaDE.filter(iii, c(filter[1], filter[2]))
      cat(" MetaDE.match(IQR probe summarize)+MetaDE.merge were performed")
      cat("\n")
      cat(
        sprintf(
          "=<quantile(mean(gene), %s), =<quantile(sd(gene), %s) were filtered",
          filter[1],
          filter[2]
        )
      )
      iii = check.tar.list(iii)
      if (forCentromics) {
        if (!is.null(x$pheno))
          iii$pheno <- x$pheno
        if (!is.null(x$meta))
          iii$meta <- x$meta
      }
      return(iii)
    }
    
  }



factorizeQuantile <- function(vec, quantileSeq = c(0, 1 / 2, 1)) {
  quantile <- quantile(vec, quantileSeq)
  factor(.bincode(vec, quantile, include.lowest = T))
}




require(sva)
.w.mergeStudyData <-
  function(ematList,
           sampleMetadata,
           batchColname = 'study',
           covariateName = NA,
           batchCorrection = TRUE,
           parPrior = TRUE,
           merge.mathod = c("combat", "reference.combat", "tdm", "c.bind", "qn"),
           ref.batch = 1,
           ...) {
    merge.mathod <- match.arg(merge.mathod) 
    geneIds = Reduce(intersect, lapply(ematList, function(x)
      rownames(x)))
    ematList2 = foreach(studyName = names(ematList)) %do% {
      ematNow = ematList[[studyName]][geneIds, ]
    }
    if (batchCorrection) {
      ematListScaled = lapply(ematList2, function(emat)
        (emat - mean(emat)) / sd(emat))
      ematMerged = do.call(cbind, ematListScaled)
      ematListScaled.pre2 <<- ematListScaled
      ematMerged.ref =  do.call(cbind, ematListScaled[1:(length(ematListScaled) -
                                                           1)])
      ematMerged.tar =  ematListScaled[[length(ematListScaled)]]
      
      merge.mathod.pre <<- merge.mathod
      
      if (is.na(covariateName)) {
        covariateInfo = model.matrix( ~ rep_len(1, ncol(ematMerged)))
      } else {
        covariateInfo = model.matrix( ~ sampleMetadata[colnames(ematMerged), covariateName])
      }
      
      if (length(unique(sampleMetadata[colnames(ematMerged), batchColname])) >
          1) {
        if (merge.mathod == "combat") {
          ematMergedNorm = sva::ComBat(
            ematMerged,
            batch = sampleMetadata[colnames(ematMerged), batchColname],
            mod = covariateInfo,
            par.prior = parPrior
          )
        }
        if (merge.mathod == "reference.combat") {
          ematMergedNorm = .w.ComBat(
            ematMerged,
            batch = sampleMetadata[colnames(ematMerged), batchColname],
            mod = covariateInfo,
            par.prior = parPrior,
            ref.batch =  ref.batch
          )
        }
        if (merge.mathod == "tdm") {
          ematMergedNorm = .w.tdm(ematMerged.ref, ematMerged.tar)
        }
        if (merge.mathod == "c.bind") {
          ematMergedNorm = cbind(ematMerged.ref, ematMerged.tar)
        }
        if (merge.mathod == "qn") {
          ematMergedNorm.tmp = batchCorrect(ematMerged.ref, ematMerged.tar, methods = "qn")
          ematMergedNorm = cbind(ematMergedNorm.tmp$qn$ref,
                                 ematMergedNorm.tmp$qn$tar)
          ematMergedNorm.pre3 <<- ematMergedNorm
          
        }
      } else {
        ematMergedNorm = ematMerged
      }
      return(ematMergedNorm)
    } else {
      return(do.call(cbind, ematList2))
    }
  }









cv.merged <-
  function(ematMerged,
           sampleMetadata,
           weights,
           alphas,
           nFolds = 10,
           foldid = NA,
           nRepeats = 3,
           yName = 'class',
           clinVarColnames = NA,
           GlobalOp = NA,
           seed = 1234,
           popSize = 30,
           maxiter = 5,
           type.min = "lambda.min",
           ...) {
    args = list(...)
    dim.temp <<- dim(ematMerged)
    
    if (!is.null(args[['family']]) & args[['family']] == 'cox') {
      y = as.matrix(sampleMetadata[colnames(ematMerged), yName])
      colnames(y) = c('time', 'status')
    } else {
      y = sampleMetadata[colnames(ematMerged), yName]
    }
 
    if (is.na(clinVarColnames[1])) {
      x = scale(t(ematMerged), center = TRUE, scale = FALSE)
    } else {
      clinVarTmp = data.frame(lapply(sampleMetadata[colnames(ematMerged), clinVarColnames], factor))
      clinDummy = model.matrix(~ 0 + ., data = clinVarTmp)
      x = cbind(scale(t(ematMerged), center = TRUE, scale = FALSE), clinDummy)
    }
    
    if (is.na(GlobalOp)) {
      if (is.na(foldid[1])) {
        cvFitList = list()
        for (ii in 1:nRepeats) {
          foldid = sample(rep(seq(nFolds), length = ncol(ematMerged)))
          cvFitList[[ii]] = foreach(alpha = alphas) %dopar% {
            set.seed(seed + alpha)
            print("alpha.seq:")
            
            cv.glmnet(
              x,
              y,
              weights = weights[colnames(ematMerged)],
              foldid = foldid,
              alpha = alpha,
              ...
            )
          }
        }
        args.temp1 <<- list(...)
        return(cvFitList)
      } else {
        cvFitList = foreach(alpha = alphas) %dopar% {
          set.seed(seed + alpha)
          cv.glmnet(
            x,
            y,
            weights = weights[colnames(ematMerged)],
            foldid = foldid[colnames(ematMerged)],
            alpha = alpha,
            
            ...
          )
        }
        return(cvFitList)
      }
      
    } else if (GlobalOp == "epsgo") {
      if (is.na(foldid[1]))
        foldid = sample(rep(seq(nFolds), length = ncol(ematMerged)))
      
      y.classes <- y
      bounds <- t(data.frame(alpha = c(0, 1)))
      colnames(bounds) <- c("lower", "upper")
      
      epsgo.cvFitList <-
        c060::epsgo(
          Q.func = ".w.tune.glmnet.interval",
          parms.coding = "none",
          seed = seed,
          show = "none",
          bounds = bounds,
          x = x,
          y = y.classes,
          foldid = foldid,
          weights = weights,
          
          type.min = type.min,
          ...
        )
      
      return(epsgo.cvFitList)
      
    } else if (GlobalOp == "GA") {
      if (is.na(foldid[1]))
        foldid = sample(rep(seq(nFolds), length = ncol(ematMerged)))
      
      GA.results <- GA.solution (
        x = x,
        y = y,
        weights = weights,
        foldid = foldid,
        seed = seed,
        popSize = popSize,
        maxiter = maxiter,
        
        ...
      )
      return(GA.results)
      
    } else {
      
    }
  }





.w.tune.glmnet.interval <-
  function (parms,
            x,
            y,
            weights,
            offset = NULL,
            lambda = NULL,
            type.measure = c("mse", "deviance", "class", "auc", "mae"),
            seed = 12345,
            nfolds = 10,
            foldid = foldid,
            grouped = TRUE,
            type.min = c("lambda.min", "lambda.1se"),
            family,
            verbose = FALSE,
            ...)
  {
    alpha <- parms[1]
    names(alpha) <- NULL
    if (verbose)
      print(paste("alpha=", alpha))
    set.seed(seed)
    cv <- cv.glmnet(
      x = x,
      y = y,
      family = family,
      alpha = alpha,
      offset = NULL,
      lambda = NULL,
      type.measure = type.measure,
      weights = weights,
      nfolds = nfolds,
      foldid = foldid,
      grouped = grouped,
      keep = T,
      intercept = T,
      standardize = F
    )
    opt.lambda <- ifelse(type.min == "lambda.min", cv$lambda.min,
                         cv$lambda.1se)
    q.val <- cv$cvm[which(cv$lambda == opt.lambda)]
    fit <- glmnet(
      x = x,
      y = y,
      family = family,
      alpha = alpha,
      lambda = opt.lambda
    )
    ret <-
      list(
        q.val = q.val,
        model = list(
          alpha = alpha,
          lambda = opt.lambda,
          nfolds = nfolds,
          cvreg = cv,
          fit = fit
        )
      )
    return(ret)
  }






require(foreach)
require(sva)
require(glmnet)
require(c060)
require(pROC)
require(Epi)
require(plot3D)
require(ggplot2)
require(precrec)
require(reshape2)
require(RColorBrewer)
require(dplyr)
require(beeswarm)
require(ROCR)
require(pamr)
require(pheatmap)
require(gridExtra)
require(GSAR)
require(org.Hs.eg.db)
require(AnnotationDbi)

options(mySeedV=76543217)

biPDSmlv2 <- function (d.merged,
                       devCohorts,
                       pw.in = NULL,
                       path = "",
                       PDS.min_exp = -10,
                       PDS.min_std = 0.2,
                       className = 'class',
                       classesTrain = c("C", "E"),
                       familyName = "binomial",
                       screen.name = "",
                       CVal = c("LOOCV", "LOSOCV", "nfold"),
                       nfold = 10,
                       GlobalOp = c("epsgo", "GA",  "none"),
                       alphas.seq = seq(0.03, 1, length = 10),
                       type.measure = "deviance",
                       type.min = c("lambda.min", "lambda.1se"),
                       geneSymbol = TRUE,
                       populations = 30,
                       generation = 5,
                       penalty.factor = rep(1, dim(d.merged[[1]]$x)[1]),
                       authors = paste(names(d.merged), "authors"),
                       seed = 29365917,
                       prior.pds.results = NULL,
                       prior.predsList = NULL,
                       browseTable = F,
                       selectedStat = T,
                       preCalpredsList = NULL,
                       valBatchCorrectMethod = "combat",
                       ...) {

  
  try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
  
  op <- par(no.readonly = T)
  set.seed(seed)
  
  GlobalOp = match.arg(GlobalOp)
  GlobalOp.pre <- GlobalOp
  type.min = match.arg(type.min)
  CVal = match.arg(CVal)
  
  valCohorts = 1 - devCohorts
  pdf.options(family = 'Helvetica', useDingbats = F)
  
  ifelse(!dir.exists("./results"),
         dir.create("./results"),
         "Folder exists already")
  
  if (!is.null(pw.in)) {
    pw.in.pre <- NULL
    pw.in.pre$pathway <- as.list(names(pw.in))
    pw.in.pre$entrez_gene_ids <- pw.in
    pw.in <- pw.in.pre
    PDS = T
  } else {
    PDS = F
  }
  
  for (i in 1:length(d.merged))
    d.merged[[i]][[2]] <-  d.merged[[i]][[2]][, 1]
  study <-
    foreach(a = 1:length(d.merged), .combine = "c") %do% {
      rep(names(d.merged)[a], length(d.merged[[a]][[2]]))
    }
  sample <-
    foreach(a = 1:length(d.merged), .combine = "c") %do% {
      colnames(d.merged[[a]][[1]])
    }
  class <-
    foreach(a = 1:length(d.merged), .combine = "c") %do% {
      d.merged[[a]][[2]]
    }
  ematList <-
    foreach(a = 1:length(d.merged)) %do% {
      d.merged[[a]][[1]]
    }
  names(ematList) <- names(d.merged)
  if (any(rownames(ematList[[a]]) == "")) {
    for (a in 1:length(ematList)) {
      ematList[[a]] = ematList[[a]][-which(rownames(ematList[[a]]) == ""), ]
    }
  }
  meta2Tmp <-
    data.frame(
      study = as.character(study),
      sample = as.character(sample),
      class = as.factor(class),
      stringsAsFactors = F
    )
  str(ematList)
  
  
  meta <-
    data.frame(
      study = as.character(names(ematList)),
      devCohorts = devCohorts,
      valCohorts = valCohorts,
      stringsAsFactors = F
    )
  intercept = TRUE
  rownames(meta) = meta[, 'study']
  meta[, 'devCohorts'] = meta[, 'devCohorts'] == 1
  meta[, 'valCohorts'] = meta[, 'valCohorts'] == 1
  rownames(meta2Tmp) = meta2Tmp[, 'sample']
  meta2 = meta2Tmp[meta2Tmp[, 'study'] %in% meta[, 'study'], ]
  mergedClass = .w.mergeStudyData(ematList[meta[meta[, 'devCohorts'], 'study']], meta2, covariateName =
                                 NA)
  
  devNames = meta[meta[, 'devCohorts'], 'study']
  idxDev = (meta2[, 'study'] %in% devNames) &
    (meta2[, className] %in% classesTrain)
  devSampleNames = rownames(meta2)[idxDev]
  meta2Dev = meta2[idxDev, ]
  meta2Dev[, className] = factor(meta2Dev[, className], levels = classesTrain)
  devMerged = mergedClass[, meta2Dev[, 'sample']]
  s.Title <- sprintf("%sAND%s", devNames[1], length(devNames) - 1)
  glmnetArgs = makeGlmnetArgs(meta2Dev)
  
  if (!is.null(prior.pds.results) || PDS) {
    if (!is.null(prior.pds.results)) {
      exp <- devMerged
      PDS.res = prior.pds.results
      insel <- meta2[which(meta2[, 'study'] %in% devNames), ]
      PDSmatrix <- do.call(rbind.data.frame, PDS.res$scores)
      colnames(PDSmatrix) <- colnames(exp)
      PDSmatrix <- as.matrix(PDSmatrix)
      PDS = T
      gene_ids <- PDS.res$gene_ids
      
    } else {
      exp <- devMerged
      pw.in <- pw.in
      gene_ids <- unique(do.call("c", pw.in$entrez_gene_ids))
      nor <- rep(F, dim(exp)[2])
      nor[which(colnames(exp) %in% meta2$sample[which(meta2$class == "C")])] <-
        T
      normal.cv.merged <- nor
      exp.cv.merged <- exp
      
      if (file.exists(file.path(
        "./results",
        sprintf("%s%s_PDS.res.RData", s.Title, path)
      ))) {
        load(file.path(
          "./results",
          sprintf("%s%s_PDS.res.RData", s.Title, path)
        ))
      } else {
        PDS.res <-
          wrap.PDS(
            exp.cv.merged,
            rownames(exp.cv.merged),
            pw.in$entrez_gene_ids,
            pw.in$pathway,
            normal.cv.merged,
            attempts = 2,
            min_exp = PDS.min_exp ,
            min_std = PDS.min_std
          )
      }
      PDSmatrix <- do.call(rbind.data.frame, PDS.res$scores)
      colnames(PDSmatrix) <- colnames(exp)
      PDSmatrix <- as.matrix(PDSmatrix)
      insel <- meta2[which(meta2[, 'study'] %in% devNames), ]
      PDS.res$PDSmatrix <- PDSmatrix
      PDS.res$insel <- insel
      PDS.res$gene_ids <- gene_ids
      if (is.null(prior.pds.results) &&
          !file.exists(file.path(
            "./results",
            sprintf("%s%s_PDS.res.RData", s.Title, path)
          )))
        save(PDS.res, file = file.path(
          "./results",
          sprintf("%s%s_PDS.res.RData", s.Title, path)
        ))
      cat(sprintf("\n\nPDSmatrix dim: %s\n\n", dim(PDSmatrix)))
    }
  }
  
  if (PDS) {
    Input.D <- "PDS"
  } else {
    Input.D <- "expr"
  }
  
  if (Input.D == "PDS") {
    Input.D.matrix <- PDSmatrix
    geneSymbol = F
  } else {
    Input.D.matrix <- devMerged
  }
  
  if (CVal == "LOOCV") {
    nRepeats = 1
    nFolds = ncol(Input.D.matrix)
    foldid = NA
  } else if (CVal == "LOSOCV") {
    nRepeats = 2
    nFolds = NA
    foldid = glmnetArgs$foldid
  } else {
    nRepeats = 1
    nFolds = nfold
    foldid = NA
  }
  
  if (GlobalOp == "none") {
    alphas = alphas.seq
    GlobalOp = NA
  }
  
  CVconfusion <-
    function(cvFit,
             lambda,
             mergedMatrix,
             meta2,
             className = 'class',
             classLevels = NA) {
      if (is.na(classLevels[1])) {
        classLevels = names(cvFit$glmnet.fit$beta)
      }
      if (class(cvFit$glmnet.fit)[1] == "lognet") {
        cvProbs = cvFit$fit.preval[, which.min(abs(cvFit$lambda - lambda)), drop =
                                     F]
        cvProbs <- cbind(1 - cvProbs, cvProbs)
      } else {
        cvProbs <-
          cvFit$fit.preval[, , which.min(abs(cvFit$lambda - lambda))]
      }
      rownames(cvProbs) <- colnames(mergedMatrix)
      colnames(cvProbs) <- cvFit$glmnet.fit$classnames
      cvProb <- cvProbs
      preds = colnames(cvProbs)[apply(cvProbs, MARGIN = 1, function(x)
        which.max(x))]
      predsFactor = factor(preds, levels = classLevels)
      trueClasses =  factor(meta2[colnames(mergedMatrix), className], levels =
                              classLevels)
      cv.preds <- preds
      cv.predsFactor <- predsFactor
      cv.trueClasses <- trueClasses
      cv.trueClasses.pre <<-
        list(trueClasses,
             predsFactor,
             meta2,
             colnames(mergedMatrix),
             className,
             classLevels)
      confus.cv <- table(trueClasses, predsFactor)
      list(
        cvProb = cvProb,
        cv.preds = cv.preds,
        cv.predsFactor = cv.predsFactor,
        cv.trueClasses = cv.trueClasses,
        confus.cv = confus.cv
      )
    }
 
  
  res <- list()
  n.fold = sprintf("nfold%s", nFolds)
  Exp. <-
    sprintf("%s.%s.%s.%s.%s.%s.%s",
            Input.D,
            s.Title,
            path,
            CVal,
            n.fold,
            GlobalOp,
            screen.name)
  savePath = file.path(getwd(), "results", Exp.)
  if (!file.exists(sprintf("%s_mergedMeta.RData", savePath))) {
    mergedMeta = list(
      x = Input.D.matrix,
      geneMatrix = devMerged,
      y = meta2[colnames(Input.D.matrix), "class"],
      meta = meta2
    )
    save(mergedMeta, file = sprintf("%s_mergedMeta.RData", savePath))
  }
  
  set.seed(seed)
  
  sel.cvFitList <-
    cv.merged(
      Input.D.matrix,
      meta2,
      yName = className,
      weights = glmnetArgs$weights,
      nRepeats = nRepeats,
      nFolds = nFolds,
      alphas = alphas,
      family = familyName,
      intercept = intercept,
      keep = TRUE,
      GlobalOp = GlobalOp,
      type.measure = type.measure,
      type.min = type.min,
      standardize = F,
      popSize = populations,
      maxiter = generation,
      seed = seed,
      foldid = foldid,
      penalty.factor = penalty.factor,
      ...
    )
  save(sel.cvFitList, file = sprintf("%s_sel.cvFitList.RData", savePath))
  
  
  GlobalOp <- GlobalOp.pre
  require(plot3D)
  if (GlobalOp == "none") {
    Crude3Dmatrix = NULL
    if (CVal == "LOSOCV")
      sel.cvFitList = list(sel.cvFitList)
    for (i in 1:length(sel.cvFitList[[1]])) {
      Crude3Dmatrix = rbind(Crude3Dmatrix,
                            c(
                              alphas[i],
                              sel.cvFitList[[1]][[i]]$lambda.min,
                              min(sel.cvFitList[[1]][[i]]$cvm)
                            ))
    }
    D = as.data.frame(Crude3Dmatrix)
    scatter3D(
      D$V1,
      D$V2,
      D$V3,
      phi = 10,
      theta = 33,
      bty = "g",
      type = "p",
      cex.lab = 1.5,
      xlab = "alpha.grid",
      ylab = "lambda.min",
      zlab = "cvm.min",
      pch = 20,
      cex = c(1.2, 1.2, 1.2),
      xlim = c(0, 1.2),
      ylim = c(-2, 30),
      zlim = c(min(D$V3) * .99, max(D$V3) * 1.02)
    )
    p1 <- recordPlot()
    D.min = D[which(D[, 3] == min(D[, 3])), ]
    alpha.min = which(D[, 3] == min(D[, 3]))
    alpha.sol = D[alpha.min, 1]
    alpha = alpha.sol
    cvFit = sel.cvFitList[[1]][[alpha.min]]
    cvFit.pre <<-  cvFit
    cvFit.pre$fit.preval
    fitResult = cvFit$glmnet.fit
    lambda = cvFit$lambda.min
  }
  if (GlobalOp == "GA") {
    DD = as.data.frame(sel.cvFitList[[3]])
    DD.min = DD[which(DD[, 3] == min(DD[, 3])), ]
    DD.min
    scatter3D(
      DD$V1,
      DD$V2,
      DD$V3,
      phi = 10,
      theta = 35,
      bty = "g",
      type = "p",
      cex.lab = 1.5,
      xlab = "alpha",
      ylab = "lambda.min",
      zlab = "cvm.min",
      pch = 20,
      cex = c(1.2, 1.2, 1.2),
      xlim = c(0, 1.2),
      ylim = c(-2, 30),
      zlim = c(min(DD$V3) * .99, max(DD$V3) * 1.02),
      main = ""
    )
    
    p1 <- recordPlot()
    alpha = round(DD[, 1], 3)
    GA.sol = round(sel.cvFitList[[2]], 3)
    cvFit = sel.cvFitList[[1]][[which(alpha == GA.sol[[1]])[1]]]
    alpha = GA.sol
    fitResult = cvFit$glmnet.fit
    lambda = cvFit$lambda.min
  }
  if (GlobalOp == "epsgo") {
    sumfit <- summary(sel.cvFitList, verbose = TRUE)
    plot(sumfit)
    p0 <- recordPlot()
    sumfit.pre <<- sumfit
    alpha = sumfit$opt.alpha
    cvFit = sumfit$opt.models$model$cvreg
    fitResult = cvFit$glmnet.fit
    lambda = sumfit$opt.lambda
    epsgo.cvFitList <- sel.cvFitList
    Crude3Dmatrix = NULL
    for (i in 1:length(epsgo.cvFitList$model.list)) {
      Crude3Dmatrix = rbind(
        Crude3Dmatrix,
        c(
          epsgo.cvFitList$model.list[[i]]$model$alpha,
          epsgo.cvFitList$model.list[[i]]$model$cvreg$lambda.min,
          min(epsgo.cvFitList$model.list[[i]]$model$cvreg$cvm)
        )
      )
    }
    DDD = as.data.frame(Crude3Dmatrix)
    scatter3D(
      DDD$V1,
      DDD$V2,
      DDD$V3,
      phi = 7,
      theta = 35,
      bty = "g",
      type = "p",
      cex.lab = 1.5,
      xlab = "alpha.grid",
      ylab = "lambda.min",
      zlab = "cvm.min",
      pch = 20,
      cex = c(1.2, 1.2, 1.2),
      xlim = c(0, 1.2),
      ylim = c(-2, 30),
      zlim = c(min(DDD$V3) * .99, max(DDD$V3) * 1.02)
    )
    p1 <- recordPlot()
    DDD.min = DDD[which(DDD[, 3] == min(DDD[, 3])), ]
  }
  
  
  alpha.lambda.sol.temp <- round(c(alpha, lambda), 4)
  cat(sprintf(
    "\n\nalpha.lambda.sol: \n %s\n",
    paste0(alpha.lambda.sol.temp, collapse = ", ")
  ))
  
  coefDf <-
    data.frame(as.matrix(coef(fitResult, s = lambda)[coef(fitResult, s = lambda)[, 1] != 0, ]))
  
  save(coefDf, file = sprintf("%s_coefDf.rda", savePath))
  coefDf <-
    data.frame('Nonzero.Features' = rownames(coefDf)[-1], Coef = coefDf[-1 , 1])
  coefDf <- coefDf[order(coefDf[, "Coef"], decreasing = T),]
  print(data.frame(
    'Nonzero.Features' = substr(coefDf[, 1], 1, 45),
    Coef = coefDf[, 2]
  ))
  attributes(coefDf)$meta <-
    paste0(savePath,
           "_alphaAndlamda_",
           paste(alpha.lambda.sol.temp, collapse = "_"))
  nonzero <- coefDf
  CVconf <-
    CVconfusion(cvFit,
                lambda,
                Input.D.matrix,
                meta2,
                className = className,
                classesTrain)
  
  cvProb <- CVconf$cvProb
  cv.trueClasses <- CVconf$cv.trueClasses
  cv.predsFactor <- CVconf$cv.predsFactor
  cv.preds <- CVconf$cv.preds
  confus.cv <- CVconf$confus.cv
  
  cv <- new(
    "binary",
    trueclass  = as.numeric(cv.trueClasses) - 1,
    predclass =  as.numeric(cv.predsFactor) - 1,
    predprob = cvProb[, 2]
  )
  cv.statistics <-
    round(as.data.frame(c(
      AUROC = AUC(cv), BRIER = brier(cv), APRFMscore(cv)
    ), drop = F), 3)
  rownames(cv.statistics) <- "CV_statistics"
  cv.statistics
  cvPrList <- list()
  cvClList <- list()
  for (i in 1:length(devNames)) {
    inx10 <- meta2$sample[which(meta2$study %in% devNames[i])]
    inx11 <- which(devSampleNames %in% inx10)
    cvPrList[[i]] <- cvProb[, 2][inx10]
    cvClList[[i]] <- cv.trueClasses[inx11]
  }
  cvPrList <- c(cvPrList, list(cvProb[, 2]))
  cvClList <- c(cvClList, list(cv.trueClasses))
  names(cvPrList) <- c(devNames, "Overall")
  names(cvClList) <- c(devNames, "Overall")
  
  msmdat <-
    mmdata(
      cvPrList,
      cvClList,
      modnames =  gsub("^s.*_", "", names(cvClList)) ,
      dsids = 1:length(cvClList)
    )
  sscurves.cv <- precrec::evalmod(msmdat)
  print(autoplot(sscurves.cv))
  
  eval.md <- precrec::evalmod(msmdat)
  aucs <- precrec::auc(sscurves.cv)
  prauc.cv <- subset(aucs, curvetypes == "PRC")
  print(knitr::kable(aucs))
  AUPRC.cv = c(prauc.cv$aucs[1:length(prauc.cv$aucs)])
  
  cv.each.val <- list()
  cutoff = .5
  for (i in names(cvClList)) {
    cv.each.val[[i]] <-
      data.frame(
        predsFactor = as.factor(ifelse(cvPrList[[i]] > cutoff, "E", "C")),
        trueClasses = cvClList[[i]],
        predsProb = cvPrList[[i]]
      )
  }
  cv.stat <- extract.stat(cv.each.val)
  cv.stat.all <- cv.stat
  
  if (length(which(valCohorts == 1)) == 0) {
    if (Input.D == "PDS") {
      result.list <-
        list(
          Exp. = Exp.,
          savePath = savePath,
          pw.in = pw.in,
          PDS.res = PDS.res,
          insel.discovery = insel,
          coefDf = coefDf,
          alpha = alpha,
          lambda = lambda,
          meta2 = meta2,
          cv.stat.all = cv.stat.all,
          CVconf = CVconf,
          cv.plot.list = cv.plot.list,
          nonzero = nonzero,
          table.cv = out.st2
        )
    } else {
      result.list  <-
        list(
          Exp. = Exp.,
          savePath = savePath,
          coefDf = coefDf,
          alpha = alpha,
          lambda = lambda,
          meta2 = meta2,
          cv.stat.all = cv.stat.all,
          CVconf = CVconf,
          cv.plot.list = cv.plot.list,
          nonzero = nonzero,
          table.cv = out.st2
        )
    }
    
    res.n = sprintf('%s_alpha_%s_lambda_%s_res.RData',
                    Exp.,
                    round(alpha, 3),
                    round(lambda, 3))
    res.path = file.path(getwd(), "results", res.n)
    save(result.list, file = res.path)
    return(result.list)
  }
}


require(cluster)

resamplingOfPamSpearman <-
  function (d = NULL,
            maxK = 7,
            reps = 10,
            pItem = 0.8,
            title = "untitled_consensus_cluster",
            seed = 1234,
            verbose = T) {
    clusterAlg = "pam"
    distance = "spearman"
    innerLinkage = "complete"
    finalLinkage = "average"
    pFeature = 1
    corUse = "everything"
    weightsItem = NULL
    weightsFeature = NULL
    ml = NULL
    
    if (!class(d) %in% c("matrix"))
      stop("d must be a matrix")
    set.seed(seed)
    connectivityMatrix <-
      function (clusterAssignments, m, sampleKey) {
        names(clusterAssignments) <- sampleKey
        cls <-
          lapply(unique(clusterAssignments), function(i)
            as.numeric(names(clusterAssignments[clusterAssignments %in% i])))
        for (i in 1:length(cls)) {
          nelts <- 1:ncol(m)
          cl <- as.numeric(nelts %in% cls[[i]])
          updt <- outer(cl, cl)
          m <- m + updt
        }
        return(m)
      }
    sampleCols <-
      function (d,
                pSamp = NULL,
                pRow = NULL,
                weightsItem = NULL,
                weightsFeature = NULL)  {
        space <- ifelse(inherits(d, "dist"), ncol(as.matrix(d)),
                        ncol(d))
        sampleN <- floor(space * pSamp)
        sampCols <- sort(sample(space, sampleN, replace = FALSE,
                                prob = weightsItem))
        this_sample <- sampRows <- NA
        if (inherits(d, "matrix")) {
          if ((!is.null(pRow)) &&
              ((pRow < 1) || (!is.null(weightsFeature)))) {
            space = nrow(d)
            sampleN = floor(space * pRow)
            sampRows = sort(sample(space, sampleN, replace = FALSE,
                                   prob = weightsFeature))
            this_sample <- d[sampRows, sampCols]
            dimnames(this_sample) <- NULL
          } else {
            
          }
        }
        return(list(
          submat = this_sample,
          subrows = sampRows,
          subcols = sampCols
        ))
      }
    
    myPal <-  function (n = 10) {
      seq = rev(seq(0, 255, by = 255 / (n)))
      palRGB = cbind(seq, seq, 255)
      rgb(palRGB, maxColorValue = 255)
    }
    setClusterColors <- function (past_ct, ct, colorU, colorList) {
      newColors = c()
      if (length(colorList) == 0) {
        newColors = colorU[ct]
        colori = 2
      }
      else {
        newColors = rep(NULL, length(ct))
        colori = colorList[[2]]
        mo = table(past_ct, ct)
        m = mo / apply(mo, 1, sum)
        for (tci in 1:ncol(m)) {
          maxC = max(m[, tci])
          pci = which(m[, tci] == maxC)
          if (sum(m[, tci] == maxC) == 1 & max(m[pci,]) ==
              maxC & sum(m[pci,] == maxC) == 1) {
            newColors[which(ct == tci)] = unique(colorList[[1]][which(past_ct ==
                                                                        pci)])
          }
          else {
            colori = colori + 1
            newColors[which(ct == tci)] = colorU[colori]
          }
        }
      }
      return(list(newColors, colori, unique(newColors)))
    }
    
    CDF <- function (ml, breaks = 100) {
      areaK = c()
      cdfDf = list()
      for (i in 2:length(ml)) {
        v = triangle(ml[[i]], mode = 1)
        h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1 / breaks))
        h$counts = cumsum(h$counts) / sum(h$counts)
        thisArea = 0
        for (bi in 1:(length(h$breaks) - 1)) {
          thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
          bi = bi + 1
        }
        areaK = c(areaK, thisArea)
        cdfDf[[i]] = list(
          cdf = data.frame(
            K = i,
            x = h$mids,
            y = h$counts
          ),
          hist = list(
            v = v,
            h = h,
            breaks = breaks
          )
        )
      }
      deltaK = areaK[1]
      for (i in 2:(length(areaK))) {
        deltaK = c(deltaK, (areaK[i] - areaK[i - 1]) / areaK[i - 1])
      }
      for (i in 2:(length(deltaK) + 1)) {
        cdfDf[[i]] = c(cdfDf[[i]], list(deltaK = deltaK[i - 1]))
      }
      cdfDf
    }
    triangle <- function (m, mode = 1) {
      n = dim(m)[1]
      nm = matrix(0, ncol = n, nrow = n)
      fm = m
      nm[upper.tri(nm)] = m[upper.tri(m)]
      fm = t(nm) + nm
      diag(fm) = diag(m)
      nm = fm
      nm[upper.tri(nm)] = NA
      diag(nm) = NA
      vm = m[lower.tri(nm)]
      if (mode == 1) {
        return(vm)
      } else if (mode == 3) {
        return(fm)
      } else if (mode == 2) {
        return(nm)
      }
    }
    ccRun <-
      function (d = d,
                maxK = NULL,
                repCount = NULL,
                diss = inherits(d, "dist"),
                pItem = NULL,
                pFeature = NULL,
                innerLinkage = NULL,
                distance = NULL,
                clusterAlg = NULL,
                weightsItem = NULL,
                weightsFeature = NULL,
                verbose = NULL,
                corUse = NULL) {
        m = vector(mode = "list", repCount)
        ml = vector(mode = "list", maxK)
        n <- ncol(d)
        mCount = mConsist = matrix(c(0), ncol = n, nrow = n)
        ml[[1]] = c(0)
        if (is.null(distance))
          distance <- "euclidean"
        acceptable.distance <- c(
          "euclidean",
          "maximum",
          "manhattan",
          "canberra",
          "binary",
          "minkowski",
          "pearson",
          "spearman"
        )
        main.dist.obj <-
          as.dist(1 - cor(d, method = distance, use = corUse))
        attr(main.dist.obj, "method") <- distance
        
        for (i in 1:repCount) {
          if (verbose) {
            message(paste("random subsample", i))
          }
          sample_x = sampleCols(d, pItem, pFeature, weightsItem,
                                weightsFeature)
          
          this_dist = NA
          boot.cols <- sample_x$subcols
          this_dist <- as.matrix(main.dist.obj)[boot.cols,
                                                boot.cols]
          this_dist <- as.dist(this_dist)
          attr(this_dist, "method") <- attr(main.dist.obj,
                                            "method")
          
          this_cluster = NA
          mCount <- connectivityMatrix(rep(1, length(sample_x[[3]])),
                                       mCount, sample_x[[3]])
          for (k in 2:maxK) {
            if (verbose) {
              message(paste("  k =", k))
            }
            if (i == 1) {
              ml[[k]] = mConsist
            }
            this_assignment = NA
            if (clusterAlg == "pam") {
              this_assignment <- pam(
                x = this_dist,
                k,
                diss = TRUE,
                metric = distance,
                cluster.only = TRUE
              )
            } else {
              this_assignment <- get(clusterAlg)(this_dist, k)
            }
            ml[[k]] <-
              connectivityMatrix(this_assignment, ml[[k]], sample_x[[3]])
          }
        }
        res = vector(mode = "list", maxK)
        for (k in 2:maxK) {
          tmp = triangle(ml[[k]], mode = 3)
          tmpCount = triangle(mCount, mode = 3)
          res[[k]] = tmp / tmpCount
          res[[k]][which(tmpCount == 0)] = 0
        }
        message("end fraction")
        return(res)
      }
    
    if (is.null(ml) == TRUE) {
      ml <-
        ccRun(
          d = d,
          maxK = maxK,
          repCount = reps,
          diss = inherits(d, "dist"),
          pItem = pItem,
          pFeature = pFeature,
          innerLinkage = innerLinkage,
          clusterAlg = clusterAlg,
          weightsFeature = weightsFeature,
          weightsItem = weightsItem,
          distance = distance,
          verbose = verbose,
          corUse = corUse
        )
    }
    res = list()
    res[[1]] = 1:2
    for (tk in 2:maxK) {
      if (verbose) {
        message(paste("consensus ", tk))
      }
      fm = ml[[tk]]
      hc = hclust(as.dist(1 - fm), method = finalLinkage)
      message("clustered")
      ct = cutree(hc, tk)
      names(ct) = colnames(d)
      if (class(d) == "dist") {
        names(ct) = colnames(as.matrix(d))
      }
      res[[tk]] = list(
        consensusMatrix = fm,
        consensusTree = hc,
        consensusClass = ct,
        ml = ml[[tk]]
      )
    }
    cdf = CDF(ml)
    res[["cdf"]] = cdf
    return(res)
  }






consensusP <-
  function(d,
           maxK = 7,
           dendColFunc = NULL,
           seed = 1234,
           reps = 10,
           verbose = T,
           heatmap = T,
           fixedK = NULL,
           ...) {
    try(dev.off(), silent = T)
    consensus.pam = resamplingOfPamSpearman(
      d,
      maxK = maxK,
      reps = reps,
      verbose = verbose,
      seed = seed,
      ...
    )
    cdfDf2 = do.call("rbind", lapply(consensus.pam[["cdf"]], function(x)
      x$cdf))
    cdfDf2 = factorizeDf(cdfDf2)
    cdfFnLower <- cdfFnUpper <- PAC <- list()
    for (i in 1:length(levels(cdfDf2$K)) + 1) {
      cdfFnLower[[i]] <-
        cdfDf2[which(cdfDf2$K == i & cdfDf2$x == 0.095), , drop = F]$y
      cdfFnUpper[[i]] <-
        cdfDf2[which(cdfDf2$K == i & cdfDf2$x == 0.895), , drop = F]$y
      PAC[[i]] = round(cdfFnUpper[[i]] - cdfFnLower[[i]], 2)
    }
    orderInx <- order(do.call("c", PAC))
    for (k in 2:length(PAC))
      if (PAC[[k]] == min(do.call("c", PAC)))
        optK = k
    if (!is.null(fixedK))
      optK = fixedK
    list(res = consensus.pam, optK = optK)
  }





require(VennDiagram)
require(MetaDE)
require(RankProd) 
require(org.Hs.eg.db)

.metaGene <-
  function(s,
           FDRcutoff = 0.05,
           nperm = NULL,
           RP.ind = T,
           heatmap.sig.genes = F,
           NewRankRrodMethods = F,
           bedArgs = list(label.sig = c(0.01, 2)),
           mcircosArgs = list(
             switchtype = "updown",
             pianoPlot = T,
             pianoPlot.height = .07
           ),
           top = 20,
           seed = seed,
           fillcatcol = c(
             "#40A4D8",
             "#33BEB7",
             "#B2C224",
             "#FECC2F",
             "#FBA127",
             "#F66320",
             "#DB3937",
             "#A463D7",
             "#0C5BCE"
           ),
           path = "") {
    check.tar.list <- function(s) {
      names.tmp = gsub(".*_", "", names(s))
      if (length(inx <-
                 which(rownames(s[[1]]$x) %in% c("", "NA"))) != 0) {
        s = lapply(s, function(x) {
          x$x <- x$x[-inx,]
          x
        })
        names(s) <- names.tmp
        s
      } else {
        names(s) <- names.tmp
        s
      }
    }
    s = check.tar.list(s)
    set.seed(seed)
    sprintf("The total number of sample from %s studies was %s",
            length(s),
            length(unlist(do.call(
              "c", lapply(s, function(x)
                x$y)
            ))))
    if (is.null(nperm) &
        length(unlist(do.call("c", lapply(s, function(x)
          x$y)))) >= 50)
      nperm = 20
    sprintf("nperm : %s", nperm)
    
    for (i in 1:length(s))
      s[[i]]$y <-
      as.numeric(as.factor(s[[i]]$y[, 1])) - 1  
    mainTitle = sprintf("%s and %s others %s", names(s)[1], (length(s) - 1), path)
    d.filtered <- s
    
    if (RP.ind) {
      RP.adv.out.ind.metaDE <- list()
      for (i in 1:length(d.filtered)) {
        RP.adv.out.ind.metaDE[[i]] <-
          MetaDE.rawdata(
            d.filtered[i],
            meta.method = c("rankProd"),
            rth = length(d.filtered),
            nperm = nperm
          )
        RP.adv.out.ind.metaDE[[i]]$up <-
          as.data.frame(RP.adv.out.ind.metaDE[[i]]$AveFC[RP.adv.out.ind.metaDE[[i]]$AveFC < 0 &
                                                           RP.adv.out.ind.metaDE[[i]]$FDR.up < FDRcutoff]) 
        RP.adv.out.ind.metaDE[[i]]$down <-
          as.data.frame(RP.adv.out.ind.metaDE[[i]]$AveFC[RP.adv.out.ind.metaDE[[i]]$AveFC > 0 &
                                                           RP.adv.out.ind.metaDE[[i]]$FDR.down < FDRcutoff]) 
        RP.adv.out.ind.metaDE[[i]]$updown  <-
          unique(c(
            rownames(RP.adv.out.ind.metaDE[[i]]$up),
            rownames(RP.adv.out.ind.metaDE[[i]]$down)
          ))
      }   
    } else {
      RP.adv.out.ind.metaDE <- NULL
    }
    
    if (length(d.filtered) == 1) {
      cat("\nAssign variable for save single cohort RP results\n")
      RP.adv.out.ind.metaDE <- RP.adv.out.ind.metaDE[[1]]
      temp.FC <-
        data.frame(updown = RP.adv.out.ind.metaDE$updown)
      rownames(temp.FC) <- temp.FC[, 1]
      
      temp.FC2 <-
        merge.rowname(temp.FC,
                      data.frame(Log2FC = -RP.adv.out.ind.metaDE$AveFC))
      
      RP.adv.out.ind.metaDE$pfp = data.frame(FDR.up = RP.adv.out.ind.metaDE$FDR.up, FDR.down =
                                               RP.adv.out.ind.metaDE$FDR.down)
      temp.FC2 <- merge.rowname(temp.FC2, RP.adv.out.ind.metaDE$pfp)
      temp.FC2[, 2] <- format(round(temp.FC2[, 2], 4), nsmall = 4)
      temp.FC2[, 3] <- round(temp.FC2[, 3], 4)
      temp.FC2[, 3] <- format.pval(temp.FC2[, 3], eps = .000001)
      
      temp.FC2[, 4] <- round(temp.FC2[, 4], 4)
      temp.FC2[, 4] <- format.pval(temp.FC2[, 4], eps = .000001)
      
      if (length(which(rownames(temp.FC2) %in% c("NA", ""))) != 0) {
        temp.FC2 <-
          temp.FC2[-which(rownames(temp.FC2) %in% c("NA", "")), ]
      }
      
      temp.FC2[, 1] <- mapg(rownames(temp.FC2))
      
      colnames(temp.FC2) <-
        c("Gene symbol",
          "Log2 FC (Class 2/Class 1)",
          "FDR (RP-up) ",
          "FDR (RP-down)")
      temp.FC2[1:10, ]
      return(temp.FC2)
    }
    MetaDE.Res <-
      MetaDE.rawdata(
        d.filtered,
        ind.method = rep("modt", length(d.filtered)),
        meta.method = c("Fisher", "Stouffer") ,
        rth = length(d.filtered),
        paired = rep(F, length(d.filtered)),
        nperm = nperm,
        asymptotic = T
      )
    MetaDE.Res.FEM <-
      MetaDE.rawdata(
        d.filtered,
        ind.method = rep("modt", length(d.filtered)),
        meta.method = c("FEM"),
        rth = length(d.filtered),
        paired = rep(F, length(d.filtered)),
        nperm = nperm,
        asymptotic = T
      )
    MetaDE.Res.REM <-
      MetaDE.rawdata(
        d.filtered,
        ind.method = rep("modt", length(d.filtered)),
        meta.method = c("REM"),
        rth = length(d.filtered),
        paired = rep(F, length(d.filtered)),
        nperm = nperm,
        asymptotic = T
      )
    MetaDE.Res.RP <-
      MetaDE.rawdata(
        d.filtered,
        meta.method = c("rankProd"),
        rth = length(d.filtered),
        nperm = nperm,
        asymptotic = F
      ) 
    MetaDE.Res.RP$pfp <-
      data.frame(MetaDE.Res.RP$FDR.up, MetaDE.Res.RP$FDR.down)
    
    MetaDE.Res.RP$up <-
      as.matrix(sort(MetaDE.Res.RP$AveFC[MetaDE.Res.RP$AveFC < 0 &
                                           MetaDE.Res.RP$pfp[, 1] < FDRcutoff])) # class1 < class2
    MetaDE.Res.RP$down <-
      as.matrix(sort(MetaDE.Res.RP$AveFC[MetaDE.Res.RP$AveFC > 0 &
                                           MetaDE.Res.RP$pfp[, 2] < FDRcutoff], decreasing = T))   # class1 < class2
    MetaDE.Res.RP$updown  <-
      unique(c(rownames(MetaDE.Res.RP$up), rownames(MetaDE.Res.RP$down)))
    Meta.Results.Raw <-
      list(
        MetaDE.Res = MetaDE.Res,
        MetaDE.Res.FEM = MetaDE.Res.FEM,
        MetaDE.Res.REM = MetaDE.Res.REM ,
        MetaDE.Res.RP = MetaDE.Res.RP
      )
    
    metaFDR <- MetaDE.Res$meta.analysis$FDR
    metaFDR.FEM <- MetaDE.Res.FEM$meta.analysis$FDR
    metaFDR.REM <- MetaDE.Res.REM$meta.analysis$FDR
    metaFDR.RP.up <- MetaDE.Res.RP$up
    metaFDR.RP.down <- MetaDE.Res.RP$down
    
    Meta.Results.FDR <-
      list(
        metaFDR = metaFDR,
        metaFDR.FEM = metaFDR.FEM,
        metaFDR.REM = metaFDR.REM,
        metaFDR.RP.up = metaFDR.RP.up,
        metaFDR.RP.down = metaFDR.RP.down
      )
    
    
    if (heatmap.sig.genes)
      heatmap.sig.genes(
        MetaDE.Res,
        meta.method = "Fisher",
        fdr.cut = FDRcutoff,
        color = "GR"
      )
    try(dev.off(), silent = T)
    mainTitle = mainTitle
    
    Fisher = rownames(metaFDR[which(metaFDR[, "Fisher"] < FDRcutoff), ])
    Stouffer = rownames(metaFDR[which(metaFDR[, "Stouffer"] < FDRcutoff), ])
    FEM = rownames(metaFDR.FEM[which(metaFDR.FEM[, "FEM"] < FDRcutoff), , drop =
                                 F])
    REM = rownames(metaFDR.REM[which(metaFDR.REM[, "REM"] < FDRcutoff), , drop =
                                 F])
    RP = MetaDE.Res.RP$updown
    
    x = list (
      Fisher = Fisher,
      Stouffer = Stouffer,
      FEM = FEM,
      REM = REM,
      RP = RP
    )
    x <- lapply(x, function(tt) {
      na.omit(mapg(tt))
    })
    null.inx <- c()
    for (i in 1:length(x)) {
      if (is.null(x[[i]])) {
        print(paste(sprintf("%s", names(x[i])), ": NULL"))
        null.inx <- c(null.inx, i)
      }
    }
    if (!is.null(null.inx))
      x <- x[-null.inx]
    .venn.diagram(x, mainTitle = mainTitle, fillcatcol = fillcatcol)  
    overlap <- VennDiagram::calculate.overlap(x)
    RP.adv.out <- MetaDE.Res.RP
    RP.adv.out.ind  <- RP.adv.out.ind.metaDE
    
    if (NewRankRrodMethods) {
      for (i in 1:length(s)) {
        s[[i]]$y <- as.numeric(as.factor(s[[i]]$y[, 1])) - 1
        rownames(s[[i]]$x) <-  gsub("///.*", "", rownames(s[[i]]$x))
      } 
      mainTitle = sprintf("   %s and %s others", names(s)[1], (length(s) -
                                                                 1))
      if (length(do.call("c", lapply(s, function(x)
        x$y))) >= 100) {
        RandomPairs = 100
      } else {
        RandomPairs = NA
      }
      tt <- list()
      tt.origin <- c()
      tt.cl <- c()
      for (i in  1:length(s))
      {
        tt[[i]] <- s[[i]][[1]]
        tt.origin <- c(tt.origin , rep(i, dim(s[[i]][[1]])[2]))
        tt.cl <- c(tt.cl, s[[i]]$y)
      }
      ttt <- do.call(cbind, tt)
      RP.adv.out <-
        RP.advance(
          ttt,
          tt.cl,
          tt.origin,
          logged = T,
          rand = 123,
          RandomPairs = RandomPairs
        ) 
      RP.adv.out.ind = list()
      pfp.cut.off <-
        FDRcutoff 
      for (i in  1:length(s)) {
        RP.adv.out.ind[[i]] <-
          RP.advance(
            s[[i]][[1]],
            s[[i]]$y,
            rep(1, length(s[[i]]$y)),
            RandomPairs = RandomPairs,
            logged = T,
            rand = 123
          )
        RP.adv.out.ind[[i]]$up <-
          RP.adv.out.ind[[i]]$AveFC[RP.adv.out.ind[[i]]$pfp[, 1] < pfp.cut.off, , drop =
                                      F] #  class1 < class2
        RP.adv.out.ind[[i]]$down <-
          RP.adv.out.ind[[i]]$AveFC[RP.adv.out.ind[[i]]$pfp[, 2] < pfp.cut.off, , drop =
                                      F] #  class1 > class2
        RP.adv.out.ind[[i]]$updown <-
          rbind(RP.adv.out.ind[[i]]$up, RP.adv.out.ind[[i]]$down)
      }   
      RP.adv.out$up <-
        RP.adv.out$AveFC[RP.adv.out$pfp[, 1] < pfp.cut.off, , drop = F] #  class1 < class2
      RP.adv.out$down <-
        RP.adv.out$AveFC[RP.adv.out$pfp[, 2] < pfp.cut.off, , drop = F] #  class1 > class2
      RP.adv.out.ind[[i]]$down <-
        RP.adv.out.ind[[i]]$AveFC[RP.adv.out.ind[[i]]$pfp[, 2] < pfp.cut.off, , drop =
                                    F]
      RP.adv.out$updown <- rbind(RP.adv.out$up, RP.adv.out$down)
      
    }  
    temp.study <- paste0('Study_', 1:length(RP.adv.out.ind))
    if (NewRankRrodMethods) {
      for (i in 1:length(temp.study))
        do.call("<-", list(temp.study[i], rownames(RP.adv.out.ind[[i]]$updown)))
      Study_1
    } else {
      for (i in 1:length(temp.study))
        do.call("<-", list(temp.study[i], RP.adv.out.ind[[i]]$updown))
      Study_1
    }
    meta.top20.pre <-
      rbind(as.matrix(RP.adv.out$up[1:top, ]), as.matrix(RP.adv.out$down[top:1, ]))
    GENENAME = mapIds(
      org.Hs.eg.db,
      keys = rownames(meta.top20.pre),
      column = "GENENAME",
      keytype = "ENTREZID",
      multiVals = "asNA"
    )
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    meta.top20 <- merge.rowname(meta.top20.pre, data.frame(GENENAME))
    meta.top20$EntrezID <- rownames(meta.top20.pre)
    meta.top20$GeneSymbol <- .mapid(meta.top20$EntrezID)
    meta.top20 <- meta.top20[, c(3, 4, 1, 2)]
    meta.top20[, 3] <- -meta.top20[, 3]
    meta.top20[, 4] <- firstup(meta.top20[, 4])
    colnames(meta.top20) <-
      c("Entrez ID",
        "Gene symbol",
        "Log2 FC (Class 2/Class 1)",
        "Gene name")
    
    table.top.out <- .table(meta.top20, add.rownames = T)
    table.top.out[, 3] = textItalic()
    temp.FC <-
      data.frame(updown = RP.adv.out$updown)
    rownames(temp.FC) <- temp.FC[, 1]
    temp.FC2 <-
      merge.rowname(temp.FC, data.frame(Log2FC = -RP.adv.out$AveFC))
    temp.FC2 <- merge.rowname(temp.FC2, RP.adv.out$pfp)
    temp.FC2[, 2] <- format(round(temp.FC2[, 2], 4), nsmall = 4)
    temp.FC2[, 3] <- round(temp.FC2[, 3], 4)
    temp.FC2[, 3] <- format.pval(temp.FC2[, 3], eps = .000001)
    
    temp.FC2[, 4] <- round(temp.FC2[, 4], 4)
    temp.FC2[, 4] <- format.pval(temp.FC2[, 4], eps = .000001)
    
    colnames(temp.FC2) <-
      c("Gene symbol",
        "Log2 FC (Class 2/Class 1)",
        "FDR (RP-up) ",
        "FDR (RP-down)")
    
    meta <- rownames(temp.FC2)
    GainLossMetaList <- list()
    for (i in 1:length(temp.study))
      GainLossMetaList[[i]] <- eval(as.symbol(temp.study[i]))
    names(GainLossMetaList) <- names(s)  
    GainLossMetaList <- c(GainLossMetaList, list(meta))
    GainLossMetaList.GeneSymbol <-
      lapply(GainLossMetaList, function(x) {
        na.omit(mapg(x))
      })
    if (length(temp.study) > 5)  {
      GainLossMetaList.GeneSymbol.unique <-
        unique(do.call("c", GainLossMetaList.GeneSymbol[1:length(temp.study)]))
      GainLossMetaList.GeneSymbol.list <- NULL
    } else {
      GainLossMetaList.GeneSymbol.unique <- NULL
      GainLossMetaList.GeneSymbol.list <-
        GainLossMetaList.GeneSymbol[1:length(temp.study)]
    }
    temp.FC2$gene.symbol <- mapg(rownames(temp.FC2))
    temp.FC2 <- na.omit(temp.FC2)
    rownames(temp.FC2) <- temp.FC2$gene.symbol
    temp.FC2$Direction <-
      ifelse(as.numeric(temp.FC2[, 2]) > 0, "UP", "DOWN")
    metaInfo <- temp.FC2[, c(6, 2:4), drop = F]
    vennGainLossMeta <-
      .vennGainLossMeta(
        GainLossMetaList.GeneSymbol[[length(GainLossMetaList.GeneSymbol)]],
        GainLossMetaList.GeneSymbol.unique,
        GainLossMetaList.GeneSymbol.list,
        metaInfo = metaInfo,
        fillcatcol = fillcatcol
      )   
    vennresults <-
      list(
        RP = RP.adv.out,
        RP.up = RP.adv.out$up,
        RP.down = RP.adv.out$down,
        FiveMetaMethodsOverlap = overlap,
        meta.top20 = meta.top20,
        Meta.Results.Raw = Meta.Results.Raw,
        Meta.Results.FDR = Meta.Results.FDR,
        Meta.Results.FDR005 = x,
        GainLossMetaList = GainLossMetaList,
        GainLossMetaList.GeneSymbol = GainLossMetaList.GeneSymbol,
        metaInfo = metaInfo
      )
    vennresults$RP$AveFC <- as.matrix(vennresults$RP$AveFC)
    metaVennRes <- vennresults
    vennresults.temp <<- vennresults
    mainTitle = sys.time(mainTitle)
    save(metaVennRes, file = file.path(getwd(), paste0(mainTitle, "_metaVennRes.RData")))
    res = list(metaVennRes = metaVennRes, bedMETA = bedMETA)
    res
  }  



require(VennDiagram)
require(grid)

.venn.diagram <-
  function(x,
           mainTitle = "5_meta_methods",
           fillcatcol = NULL) {
    if (is.null(fillcatcol))
      fillcatcol = c(
        "dodgerblue",
        "darkorange1",
        "seagreen3",
        "orchid3",
        "goldenrod1",
        "darkblue",
        "darkred",
        "darkgreen",
        "gray50",
        "red",
        "blue"
      )
    try(dev.off(), silent = T)
    temp <- venn.diagram(
      x = x,
      filename = NULL,
      col = "black",
      fill = fillcatcol[1:length(x)],
      alpha = 0.50,
      cat.col = fillcatcol[1:length(x)],
      cat.cex = 1.5,
      ext.text = T,
      ext.length = 0.6,
      label.col = "gray10",
      lwd = rep(0.5, length(x)),
      lty = c(2:length(x), 1),
      cex = 1,
      fontface = "bold",
      fontfamily = "sans",
      cat.fontfamily = "sans",
      print.mode = c("raw"),
      cat.fontface = "bold",
      margin = 0.05,
      main = NULL,
      sub = NULL,
      main.cex = .1,
      sub.cex = .1
    )
    grid.draw(temp)
    pdf(file = file.path(getwd(), sprintf('%s.Venn.pdf', mainTitle)),
        6,
        5,
        useDingbats = F)
    grid.draw(temp)
    dev.off()
  }









SGSES <-
  function(df,
           pw.in,
           es.score = c("deltaD.max", "deltaD.sum"),
           scaleFreeConstant = 1) {
    es.score = match.arg(es.score)
    if (es.score == "deltaD.sum")
      scaleFreeConstant = 0.25
    es <- matrix(NA, length(pw.in), dim(df)[2])
    for (i in 1:dim(df)[2]) {
      df.tar = df[order(df[, i, drop = F], decreasing = T), i, drop = F]
      rankedNamedVec = df.tar[, 1]
      names(rankedNamedVec) = rownames(df.tar)
      for (j in 1:length(pw.in)) {
        subsetOfVec <- pw.in[[j]]
        subsetOfVec <- intersect(subsetOfVec, names(rankedNamedVec))
        PDF.Hit <- PDF.Miss <- rep(0, length(rankedNamedVec))
        Hits <- names(rankedNamedVec) %in% subsetOfVec
        PDF.Hit[Hits] <- abs(rankedNamedVec[Hits]) ^ scaleFreeConstant
        PDF.Hit <-  PDF.Hit / sum(PDF.Hit)
        CDF.Hit <- cumsum(PDF.Hit)
        PDF.Miss[!Hits] <-
          1 / (length(rankedNamedVec) - length(subsetOfVec))
        CDF.Miss <- cumsum(PDF.Miss)
        deltaD = CDF.Hit - CDF.Miss
        if (es.score == "deltaD.max") {
          deltaD.res = deltaD[which(abs(deltaD) == max(abs(deltaD)))]
        } else {
          deltaD.res = sum(deltaD)
        }
        es[j, i] <- deltaD.res
      }
    }
    rownames(es) <- names(pw.in)
    colnames(es) <- colnames(df)
    es
  }







batchCorrect <-
  function(ref,
           tar = NULL,
           methods = c("eb", "qn", "ebRef", "qnRef", "tdmRef"),
           personalized = F,
           breaks.n = 40) {
    require(sva)
    require(preprocessCore)
    require(foreach)
    methods <-  match.arg(methods, several.ok = T)
    print(sprintf("Methods : %s", methods))
    res <- list()
    op <- par(no.readonly = T)
    if (class(ref) == c("list")) {
      if (length(ref) == 1) {
        res[["ebCorrectedRef"]] <- ref[[1]]
        return(res)
      }
      geneIds = Reduce(intersect, lapply(ref, function(x)
        rownames(x)))
      ematList2 = foreach(a = names(ref)) %do% {
        ematNow = ref[[a]][geneIds, ]
      }
      ematListScaled = lapply(ematList2, function(emat)
        (emat - mean(emat)) / sd(emat)) 
      ematMerged = do.call("cbind", ematListScaled)
      colnames(ematMerged)
      mod <- data.frame("(Intercept)" = rep(1, ncol(ematMerged)))
      
      mod <-
        tryCatch({
          rownames(mod) <- colnames(ematMerged)
          mod
        }, error = function(e) {
          mod <- NULL
          mod
        })
      ii <- batch.vec <- list()
      for (i in 1:length(ref)) {
        ii[[i]] <- rep(i, dim(ref[[i]])[2])
        batch.vec[[i]] <- rep(names(ref)[i], dim(ref[[i]])[2])
      }
      batch <- do.call("c", ii)
      if(methods =="ebRef") { eb.ref.batch = 1 } else { eb.ref.batch = NULL }
      ebCorrectedRef = ComBat(ematMerged, batch = batch, mod = mod, ref.batch = eb.ref.batch)
      res[["ebCorrectedRef"]] <- ebCorrectedRef
      
      par(mfrow = c(2, 3))
      for (k0 in 1:length(ref)) {
        hist(
          ref[[k0]],
          prob = T,
          breaks = breaks.n,
          col = "skyblue",
          border = "#CCF0FF",
          xlab = "",
          main = sprintf("%s", names(ref)[k0])
        )
      }
      for (k in unique(batch)) {
        hist(
          ebCorrectedRef[, batch == k],
          prob = T,
          breaks = breaks.n,
          col = "skyblue",
          border = "#CCF0FF",
          xlab = "",
          main = sprintf("%s (corrected)", names(ref)[k])
        )
      }
      col = paste0(.col(wypark)(1:length(ref)), "80")
      names(col) <- names(res)
      hist(
        ebCorrectedRef,
        prob = T,
        breaks = breaks.n,
        col = "skyblue",
        border = "#CCF0FF",
        main = "",
        xlab = ""
      )
      lines(density(ebCorrectedRef), col = "#333333", xpd = T)
      for (k in unique(batch)) {
        lines(density(ebCorrectedRef[, batch == k]),
              col = col[k],
              xpd = T)
        
      }
      title(main = "EB-Corrected Distribution")
      legend(
        "topright",
        names(ref),
        col = col,
        lty = 1,
        bty = "n"
      ) # fill=c("#33333366",fill), pch=20
      par(op)
      
      methods <- "eb"
      ref <- ebCorrectedRef
      if (is.null(tar))
        return(res)
    }
    
    
    if (class(tar) == c("list")) {
      ref <- tar
      geneIds = Reduce(intersect, lapply(ref, function(x)
        rownames(x)))
      ematList2 = foreach(a = names(ref)) %do% {
        ematNow = ref[[a]][geneIds, ]
      }
      str(ematList2)
      ematListScaled = lapply(ematList2, function(emat)
        (emat - mean(emat)) / sd(emat))
      ematMerged = do.call("cbind", ematListScaled)
      colnames(ematMerged)
      mod <- data.frame("(Intercept)" = rep(1, ncol(ematMerged)))
      rownames(mod) <- colnames(ematMerged)
      ii <- list()
      for (i in 1:length(ref)) {
        ii[[i]] <- rep(i, dim(ref[[i]])[2])
      }
      batch <- do.call("c", ii)
      ebCorrectedTar = ComBat(ematMerged, batch = batch, mod = mod)
      res[["ebCorrectedTar"]] <- ebCorrectedTar
      methods <- "eb"
      tar <- ebCorrectedTar
    }
    
    inx <- intersect(rownames(ref), rownames(tar))
    ref <- ref[inx, , drop = F]
    tar.pre <- tar[inx, , drop = F]

    res.ind <- list()
    if (personalized &&
        !is.null(tar)) {
      n = ncol(tar.pre)
    } else {
      n = 1
    }
    for (j in 1:n) {
      if (!personalized && !is.null(tar))
        j = 1:ncol(tar.pre)
      if (!is.null(tar))  {
        tar <- tar.pre[, j, drop = F]
        dataMat <- cbind(ref, tar)
        batch <- c(rep(1, ncol(ref)), rep(2, ncol(tar)))
      } else {
        temp <- normalize.quantiles(ref)
        rownames(temp) <- rownames(ref)
        res[["qn"]] <- list(ref = temp)
        return(res)
      }
      
      dataMat.pre <<- dataMat
      batch.pre <<- batch
      
      
      if (any(methods %in% "eb")) {
        mod <- data.frame("(Intercept)" = rep(1, ncol(dataMat)))
        rownames(mod) <- colnames(dataMat)
        temp <- ComBat(dataMat, batch, mod = mod)
        res[["eb"]] <-
          list(ref = temp[, 1:ncol(ref)], tar = temp[, (ncol(ref) + 1):ncol(temp)])
      }
      if (any(methods %in% "ebRef")) {
        mod <- data.frame("(Intercept)" = rep(1, ncol(dataMat)))
        rownames(mod) <- colnames(dataMat)
        temp <- ComBat(dataMat, batch, mod = mod, ref.batch = 1)
        res[["ebRef"]] <-
          list(ref = temp[, 1:ncol(ref)], tar = temp[, (ncol(ref) + 1):ncol(temp)])
      }
      if (any(methods %in% "qn")) {
        temp <- normalize.quantiles(dataMat)
        rownames(temp) <-
          rownames(ref)
        colnames(temp) <- c(colnames(ref), colnames(tar))
        res[["qn"]] <-
          list(ref = temp[, 1:ncol(ref)], tar = temp[, (ncol(ref) + 1):ncol(temp)])
      }
      if (any(methods %in% "qnRef")) {
        ref.temp <- normalize.quantiles.determine.target(ref)
        tar.temp <- normalize.quantiles.use.target(tar, ref.temp)
        res[["qnRef"]] <- list(ref = ref, tar = tar.temp)
      }
      if (any(methods %in% "tdmRef")) {
        cat(
          sprintf(
            "%s ref, must be log2 transformed microarray data \n tar, must be NGS data (not log transformed) \n\n\n ", cat.caution
          )
        )
        tar.temp <- TDM(ref, tar, plot = F)
        res[["tdmRef"]] <- list(ref = ref, tar = tar.temp)
      }
      if (!personalized) {
        return(res)
      } else {
        res.ind[[j]] <- res
      }
    } 
    return(res.ind)
    ref.z <-
      lapply(res.ind, function(x)
        lapply(x, function(x)
          t(scale(t(
            x$ref
          )))))
    
    tar.z <-
      lapply(res.ind, function(x)
        lapply(x, function(x) {
          u = apply(x$ref, 1, mean)
          sd = apply(x$ref, 1, sd)
          z = (x$tar - u) / sd
        }))
    names(ref.z) <- names(tar.z) <- colnames(tar.pre)
    
    ind.ref.tar <- ind.ref.tar2 <- list()
    for (k in names(tar.z)) {
      for (l in names(tar.z[[1]])) {
        ind.ref.tar[[l]] = cbind(ref.z[[k]][[l]], l = tar.z[[k]][[l]])
      }
      ind.ref.tar2 [[k]] = ind.ref.tar
    }
    par(mfrow = c(1, 1))
    col = c(rainbow(length(res)))
    names(col) <- names(res)
    
    hist(
      ref,
      prob = T,
      col = "skyblue",
      border = "#CCF0FF",
      breaks = breaks.n,
      main = NULL,
      xpd = T
    )
    
    lines(density(ref), col = "skyblue",  xpd = T)
    for (i in names(res)) {
      lines(density(res[[i]]$tar), col = col[i],  xpd = T)
      
    }
    title(main = sprintf(
      "Comparison of %s \n to the Reference Distribution",
      paste(names(res), collapse = " ,")
    ))
    legend("topright",
           names(res),
           col = col,
           lty = 1,
           bty = "n") 
    par(op)
    
    res
    
  }











add.gs <-
  function(x,
           y,
           gs.list,
           func = c("eigene", "mean", "lassoCoxRiskScore", "RFRiskScore", "custom"),
           add.gs.quantile = list(median = c(0, 1 / 2, 1),
                                  tertiles = c(0, 1 / 3, 2 / 3, 1)),
           custom.func = NULL) {
    if (colnames(x) != rownames(y))
      stop("should be : colnames(x) == rownames(y)")
    require(survival)
    require(randomForestSRC)
    e <- new.env()
    eigene = F
    if (class(func) == "function") {
      func = func
    } else {
      func = match.arg(func)
      if (func == "eigene") {
        func = function(x, y) {
          svd(t(scale(t(x))))$v[, 1]
        }
        eigene = T
        cat("\nSame value as WGCNA::moduleEigengenes\n")
      } else if (func == "mean") {
        func = function(x, y) {
          colMeans(x)
        }
        cat("\nFunc: means\n")
      } else if (func == "lassoCoxRiskScore") {
        func = function(x, y) {
          re <- survG(x, y, 0, 1, glmnet.alpha = 1)
          
          assign("riskScoreFormula", re, envir = e)
          print(re)
          re$RiskScore[1, ]
        }
        gs.list = list(lassoCoxRiskScore = rownames(x))
        cat("covariates: X <- cbind(t(x), y[,-c(1:2), drop=F])")
      } else if (func == "RFRiskScore") {
        func = function(x, y) {
          RSF <-
            rfsrc(Surv(time, status) ~ .,
                  cbind(t(x), y),
                  ntree = 100,
                  importance = T)
          cat("sum of the risk output = Risk score\n")
          names(RSF$predicted.oob) <-
            rownames(y)
          RSF$predicted.oob
        }
        gs.list = list(RFRiskScore = rownames(x))
      } else {
        return(cat("\nfunc should be function\n"))
      }
    }
    
    signature.xy <- tar.x.out <- tar.y.out <- MissingGenes <- list()
    
    for (i in names(gs.list)) {
      inx11 <- gs.list[[i]]
      x.temp <- x[intersect(inx11, rownames(x)), , drop = F] 
      x.temp.tmp <<- x.temp
      
      MissingGenes[[i]] <- setdiff(inx11, rownames(x.temp))
      cat(sprintf(
        "\nMissing Genes: %s ===================\n\n",
        paste(MissingGenes[[i]], collapse = ", ")
      ))
      pc1 <- func(x.temp, y)
      print(pc1)
      
      if (eigene) {
        averExpr = rowMeans(scale(t(x.temp)), na.rm = T) # sample mean
        cor = cor(averExpr, pc1, use = "p")
        if (cor < 0)
          pc1 = -pc1
      }
      
      tar.x.out[[i]] <- pc1
      
      signature.xy[["x"]] <- x
      for (j in names(add.gs.quantile)) {
        quantile <- quantile(tar.x.out[[i]], add.gs.quantile[[j]])
        tar.y.out[[paste0(i, "-", j)]] <-
          factor(.bincode(tar.x.out[[i]], quantile, include.lowest = T))
      }
      names(tar.x.out[[i]]) <- names(gs.list)[i]
    }
    signature.xy[["y"]] <-
      Reduce(cbind, list(y, data.frame(tar.x.out), data.frame(tar.y.out)))
    signature.xy[["MissingGenes"]] <- MissingGenes
    
    if (exists("riskScoreFormula", envir = e))
      signature.xy[["riskScoreFormula"]] <- e$riskScoreFormula
    signature.xy
    
  }







survG <-
  function (x,
            y,
            outer.split = .3,
            n.outerCV = 2,
            glmnet.alpha = 1,
            glmnet.standardize = T,
            superpc = F,
            superpc.nPC = floor(dim(x)[2] / 2),
            plot = T,
            seed = "202008",
            showNullHR = F,
            showNullHR.n.perm = 10,
            ...) {
    require(glmnet)
    set.seed(seed)
    if (superpc) {
      tmp <-
        list(
          x = x,
          y = y[, "time"],
          censoring.status = y[, "status"],
          eaturenames = rownames(x)
        )
      tmp2 <-
        names(sort(
          abs(
            superpc::superpc.train(tmp, type = "survival")$feature.scores
          ),
          decreasing = T
        ))[1:superpc.nPC]
      x <- x[tmp2,]
    }
    n.genes <- nrow(x)
    n.samples <- ncol(x)
    
    splitSurvGroup <-
      function(riskScore,
               y,
               plot = F,
               cut = median(riskScore)) {
        gs <-
          cut(riskScore,
              breaks = c(-Inf, cut, Inf),
              c("Low Risk", "High Risk"))
        gs
        levels(gs)[1] == "Low Risk"
        y$gs <- gs
        if (plot) {
          par(
            mfrow = c(2, 2),
            mar = c(5, 4, 4, 1),
            bty = "n",
            xpd = F
          )
          plot(
            riskScore,
            riskScore,
            ylab = "riskScore",
            main = "Risk Group (median cut)",
            ask = F
          )
          points(riskScore[gs == "Low Risk"], riskScore[gs == "Low Risk"], col = "grey")
          points(riskScore[gs == "High Risk"], riskScore[gs == "High Risk"], col = "red")
          abline(cut, 0, lty = 2)
          legend(
            "bottomright",
            col = c("grey", "red"),
            c("Low Risk",  "High Risk"),
            pch = 1,
            bty = "n"
          )
          fit <-
            survfit(Surv(time, status == 1) ~ gs, data = y)
          fit
          summary(fit)
          plot(
            fit,
            ylab = "Survival",
            xlab = "Time",
            col = c("grey", "red"),
            main = "Kaplan-Meier curves "
          )
          boxplot(
            y[, "time"] ~ y[, "gs"],
            notch = T,
            varwidth = T,
            ylab = "Survival Time",
            main = "Low vs High Risk groups",
            col = c("grey", "red")
          )
          
        }
        if (dim(y)[2] > 3) {
          y <- cbind(y[, c("time", "status", "gs")], y[, c(3:(ncol(y) - 1))])
          SurFit <-
            coxph(Surv(time, status == 1) ~ ., data = y)
          SurFit
        } else {
          SurFit <- coxph(Surv(time, status == 1) ~ gs, data = y)
          SurFit
        }
        return(list(SurFit = SurFit, gs = gs))
      }
    
    
    if (outer.split == 0) {
      glmnet.surv <-
        function(x,
                 y,
                 glmnet.alpha = 1,
                 glmnet.standardize = T,
                 plot = T) {
          genes.names <- rownames(x)
          if (dim(y)[2] > 2) {
            X <- cbind(t(x), y[, -c(1:2), drop = F])
            penafac <-
              c(rep(0, ncol(y[, -c(1:2), drop = F])), rep(1, nrow(x)))
          } else {
            X <- t(x)
            penafac <- rep(1, ncol(X))
          }
          y[, "time"][y[, "time"] <= 0] <-
            quantile(y[, "time"], probs = 0.02)
          COXNET.cv <-
            cv.glmnet(
              x = as.matrix(X),
              y = Surv(as.vector(y[, "time"]), as.vector(y[, "status"]) == 1),
              family = "cox",
              alpha = glmnet.alpha,
              nlambda = 100,
              penalty.factor = penafac,
              standardize = glmnet.standardize
            )
          COXNET <-
            glmnet(
              x = as.matrix(X),
              y = Surv(y[, "time"], y[, "status"] ==  1),
              family = "cox",
              alpha = glmnet.alpha,
              nlambda = 100,
              penalty.factor = penafac,
              standardize = glmnet.standardize
            )
          lambda <- COXNET.cv$lambda.min
          seq.lambda <- COXNET.cv$lambda
          ind.lambda <- which(seq.lambda == COXNET.cv$lambda.min)
          betahat <- coef(COXNET.cv, s = "lambda.min")
          print(betahat)
          betaFP <- betahat[betahat[, 1] != 0,]
          if (!is.null(dim(betaFP)))
            stop("Too Many variables are selected !!! try to increase alpha ")
          
          if (dim(y)[2] > 2) {
            gm <- setdiff(names(betaFP), colnames(y[, -c(1:2), drop = F]))
          } else {
            gm <- names(betaFP)
          }
          print(gm)
          ng <- length(gm)
          lab <- COXNET.cv$lambda.min
          bb = 0
          
          while (ng < 1) {
            bb <- bb + 1
            betahat <- coef(COXNET.cv, s = seq.lambda[ind.lambda + bb])
            betaFP <- betahat[betahat[, 1] != 0,]
            if (dim(y)[2] > 2) {
              gm <- setdiff(names(betaFP), colnames(y[, -c(1:2), drop = F]))
            } else {
              gm <- names(betaFP)
            }
            lab <- seq.lambda[ind.lambda + bb]
            ng <- length(gm)
          }
          RiskScore <- betaFP[gm] %*% t(X[, gm])
          if (plot) {
            par(mfrow = c(2, 2), mar = c(5, 4, 4, 1))
            plot(
              log(COXNET.cv$lambda),
              COXNET.cv$cvm,
              main = paste("Alpha = ", glmnet.alpha, sep = ""),
              ylim = c(min(COXNET.cv$cvlo), max(COXNET.cv$cvup)),
              xlab = "log(lambda)",
              ylab = "Partial Likelihood Deviance",
              pch = 19,
              col = "red"
            )
            for (i in 1:length(COXNET.cv$cvm))
              lines(log(c(
                COXNET.cv$lambda[i], COXNET.cv$lambda[i]
              )), c(COXNET.cv$cvlo[i], COXNET.cv$cvup[i]))
            plot(COXNET, xvar = "lambda", label = T)
            abline(
              v = log(COXNET.cv$lambda.min),
              lwd = 2,
              lty = 2,
              col = "black"
            )
            abline(
              v = log(lab),
              lwd = 2,
              lty = 2,
              col = "red"
            )
          }
          Results <- splitSurvGroup(RiskScore, y, plot = plot)
          return(list(
            RiskScore = RiskScore,
            Results = Results,
            betaS = betaFP,
            gm = gm
          ))
        }
      y.tmp <<- y
      obs <-
        glmnet.surv(x,
                    y,
                    glmnet.alpha = glmnet.alpha,
                    glmnet.standardize = glmnet.standardize)
      
      
      if (showNullHR) {
        n.perm = showNullHR.n.perm
        HRlowPerm <- matrix(NA, nrow = n.perm, ncol = 1)
        ind.gene <- matrix(NA, nrow = n.perm, ncol = n.samples)
        for (i in 1:n.perm) {
          ind.gene[i, ] <- 1:n.samples
        }
        for (i in 1:n.perm) {
          ind.gene[i, ] <- sample(c(1:n.samples), replace = F)
        }
        for (i in 1:n.perm) {
          Temp <- NA
          Temp <-
            glmnet.surv(
              x[, ind.gene[i, ]],
              y,
              plot = F,
              glmnet.alpha = glmnet.alpha,
              glmnet.standardize = glmnet.standardize
            )
          Temp.tmp <<- Temp
          if ((!is.na(Temp))[1]) {
            HRlowPerm[i, ] <- summary(Temp$Results$SurFit)[[8]][1, 1]
          }
          if ((is.na(Temp))[1])
            HRlowPerm[i, ] <- NA
          message(i)
        }
        HRlowObs <- summary(obs$Results$SurFit)[[8]][1, 1]
        print(HRlowPerm)
        print(HRlowObs)
        
        HR <- HRlowPerm[, 1]
        HR <- na.exclude(HR)
        n <- n.perm
        vv = HRlowObs[1]
        pvalue <- sum(vv > HR) / n
        dotsCall <- substitute(list(...))
        ll <- eval(dotsCall)
        if (!hasArg("xlab"))
          ll$xlab <-
          paste("Estimated HR \n Emperical p-value: ",
                round(pvalue, 4),
                sep = "")
        if (!hasArg("ylab"))
          ll$ylab <- ""
        ll$main <-
          "Null Distribution of HR on Permuted Data \n for High risk group"
        if (!hasArg("cex.lab"))
          ll$cex.lab <- 0.8
        if (!hasArg("cex.main"))
          ll$cex.main <- 1
        if (!hasArg("col"))
          ll$col <- 1
        if (!hasArg("ylim"))
          ll$ylim <- c(0, 4)
        print(ll)
        ll$x <- density(HR, from = 0, to = (max(HR) + 0.25))
        do.call("plot", args = ll)
        abline(v = vv, col = 'red')
        #CI for permuated cases
        qq <- quantile(sort(HR), prob = c(0.05, 0.95))
        abline(v = qq[1],
               col = "#00000080",
               lty = 3)
        abline(v = qq[2],
               col = "#00000080",
               lty = 3)
        abline(v = median(HR), col = "#00000080")
      }
      
      return(obs)
      
    } else {
      glmnet.n.cv <- n.outerCV
      glmnet.fold <- floor(1 / outer.split)
      Run.Time <- rep(NA, glmnet.n.cv)
      n.train <- (n.samples - floor(n.samples / glmnet.fold))
      n.test <- floor(n.samples / glmnet.fold)
      cv.train <- matrix(0, glmnet.n.cv, n.train)
      cv.test <- matrix(0, glmnet.n.cv, n.test)
      n.g <- rep(0, glmnet.n.cv)
      lambda <- rep(NA, glmnet.n.cv)
      pld <- rep(NA, glmnet.n.cv)
      HRT <- matrix(NA, nrow = glmnet.n.cv, ncol = 3)
      HRTE <- matrix(NA, nrow = glmnet.n.cv, ncol = 3)
      message(sprintf("%s outer CV are being runing ...", n.outerCV))
      pIndex <- c(1:n.samples)
      y[, "time"][y[, "time"] <= 0] <-
        quantile(y[, "time"], probs = 0.02)
      if (dim(y)[2] > 2) {
        AllX <- cbind(t(x), y[, -c(1:2), drop = F])
        penafac <-
          c(rep(0, ncol(y[, -c(1:2), drop = F])), rep(1, nrow(x)))
      } else {
        AllX <- cbind(t(x))
        penafac <- rep(1, ncol(AllX))
      }
      
      perProgFact <- NULL
      coef.mat <-
        gene.mat <- matrix(0, nrow = glmnet.n.cv, ncol = nrow(x))
      for (i in 1:glmnet.n.cv) {
        message(i, appendLF = FALSE)
        cv.train[i,] <- sort(sample(pIndex, n.train, replace = F))
        cv.test[i,] <-
          c(1:n.samples)[-c(intersect(cv.train[i, ], c(1:n.samples)))]
        Stime = y[, "time"][cv.train[i,]]
        sen = y[, "status"][cv.train[i,]]
        X <- AllX[cv.train[i,],]
        Run.Time[i] <-
          system.time(
            COXNET.cv <-
              cv.glmnet(
                x = as.matrix(X),
                y = Surv(Stime, sen == 1),
                family = "cox",
                alpha = glmnet.alpha,
                nlambda = 100,
                penalty.factor = penafac
              )
          )[3]
        lab <- COXNET.cv$lambda.min
        seq.lambda <- COXNET.cv$lambda
        ind.lambda <- which(seq.lambda == COXNET.cv$lambda.min)
        pld[i] <-
          COXNET.cv$cvm[COXNET.cv$lambda == COXNET.cv$lambda.min]
        betahat <- coef(COXNET.cv, s = "lambda.min")
        betahat.tmp <<- betahat
        betaFP <- betahat[betahat[, 1] != 0,]
        BetaFP <- betaFP
        if (class(BetaFP)[1] == "dgCMatrix") {
          betaFP <- attributes(BetaFP)$x
          names(betaFP) <- attributes(BetaFP)$Dimnames[[1]]
        }
        if (dim(y)[2] > 2) {
          gm <- setdiff(names(betaFP), colnames(y[, -c(1:2)]))
        } else {
          gm <- names(betaFP)
        }
        print(gm)
        n.g[i] <- length(gm)
        bb = 0
        while (n.g[i] < 1) {
          bb <- bb + 1
          betahat <- coef(COXNET.cv, s = seq.lambda[ind.lambda +  bb])
          betaFP <- betahat[betahat[, 1] != 0,]
          BetaFP <- betaFP
          if (class(BetaFP)[1] == "dgCMatrix") {
            betaFP <- attributes(BetaFP)$x
            names(betaFP) <- attributes(BetaFP)$Dimnames[[1]]
          }
          if (dim(y)[2] > 2) {
            gm <- setdiff(names(betaFP), colnames(y[, -c(1:2), drop = F]))
          } else {
            gm <- names(betaFP)
          }
          n.g[i] <- length(gm)
          lab <- seq.lambda[ind.lambda + bb]
        }
        
        lambda[i] <- lab
        gene.mat[i, is.element(rownames(x), gm)] <- 1
        coef.mat[i, is.element(rownames(x), gm)] <- betaFP[gm]
        score2.lassoT <- betaFP[gm] %*% x[gm, cv.train[i, ]]
        score2.lassoTE <- betaFP[gm] %*% x[gm, cv.test[i,]]
        
        y.train <- y[cv.train[i,], , drop = F]
        Results1 <-
          splitSurvGroup(score2.lassoT, y.train, plot = plot)
        HRT[i,] <- summary(Results1$SurFit)[[8]][1, 1]
        y.test <- y[cv.test[i,], , drop = F]
        Results2 <-
          splitSurvGroup(score2.lassoTE, y.test, plot = plot)
        HRTE[i,] <- summary(Results2$SurFit)[[8]][1, 1]
      }
      
      DistHR <- data.frame(HRTrain = HRT[, 1], HRTest = HRTE[, 1])
      colnames(DistHR) <- c("Train", "Test")
      dotsCall <- substitute(list(...))
      ll <- eval(dotsCall)
      if (!hasArg("xlab"))
        ll$xlab <- ""
      if (!hasArg("ylab"))
        ll$ylab <- "HR estimate"
      ll$main <-
        "Distribution of HR on Test and Train Data \n for high risk group"
      if (!hasArg("cex.lab"))
        ll$cex.lab <- 0.8
      if (!hasArg("cex.main"))
        ll$cex.main <- 1
      if (!hasArg("col"))
        ll$col <- c("grey", "red")
      ll$x <- DistHR
      ll$notch = T
      ll$varwidth = T
      do.call(boxplot, args = ll)
      Freq = colSums(gene.mat)
      names(Freq) <- rownames(x)
      sFreq <- sort(Freq, decreasing = TRUE)
      sFreq <- sFreq[sFreq > 0]
      maxG <- length(sFreq)
      if (maxG > 30)
        maxG <- 30
      dotsCall <- substitute(list(...))
      lll <- eval(dotsCall)
      lll$height <- sFreq[1:maxG]
      if (!hasArg("xlab"))
        lll$xlab <- ""
      if (!hasArg("ylab"))
        lll$ylab <- "Frequency"
      if (!hasArg("main"))
        lll$main <- "Mostly Selected Genes"
      lll$col <- "darkblue"
      lll$names <- names(sFreq)[1:maxG]
      if (!hasArg("cex.lab"))
        lll$cex.lab <- 1
      if (!hasArg("las"))
        lll$las <- 2
      lll$cex.names <- 0.65
      do.call(barplot, args = lll)
      cat("Cross Valdiated Results for Lasso and Elastic Net based Predictive Gene signature\n")
      cat("Number of CV: ", length(lambda), "\n")
      cat("Estimated  quantiles of HR on test data\n")
      print(quantile(HRTE[, 1], probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
      cat("\n")
      cat("Estimated quantiles of HR on train data\n")
      print(quantile(HRT[, 1], probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
      cat("Mostly selected 30 features:\n")
      Freq = colSums(gene.mat)
      names(Freq) <- rownames(x)
      sFreq <- sort(Freq, decreasing = TRUE)
      sFreq <- sFreq[sFreq > 0]
      maxG <- length(sFreq)
      if (maxG > 30)
        maxG <- 30
      print(names(sFreq)[1:maxG])
      invisible(
        list(
          coef.mat = coef.mat,
          beta.train = betaFP[gm],
          Run.Time = Run.Time,
          lambda = lambda,
          n.g = n.g,
          gene.mat = gene.mat,
          HRT = HRT,
          HRTE = HRTE,
          pld = pld,
          xdata = x
        )
      )
    }
  }








