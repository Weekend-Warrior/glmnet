cv.coxnet <-
function (outlist, lambda, x, y, weights, offset, foldid, type.measure, 
    grouped, keep = FALSE) 
{
    typenames = c(deviance = "Partial Likelihood Deviance",
                  auc = 'AUC')
    if (type.measure == "default") 
        type.measure = "deviance"
    if (!match(type.measure, c("deviance", 'auc'), FALSE)) {
        warning("Only 'deviance' and 'auc' available for Cox models; changed to type.measure='deviance'")
        type.measure = "deviance"
    }
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    N = nrow(y)
    nfolds = max(foldid)
    
    if ((N/nfolds < 10) && type.measure == "auc") {
      warning("Too few (< 10) observations per fold for type.measure='auc' in cv.coxnet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",
              call. = FALSE)
      type.measure = "deviance"
    }
    
    if ((length(weights)/nfolds < 10) && !grouped) {
        warning("Option grouped=TRUE enforced for cv.coxnet, since < 3 observations per fold", 
            call. = FALSE)
        grouped = TRUE
    }
    cvraw = matrix(NA, nfolds, length(lambda))
###We dont want to extrapolate lambdas on the small side
  mlami=max(sapply(outlist,function(obj)min(obj$lambda)))
  which_lam=lambda >= mlami
  if (keep | type.measure == 'auc') predmat = matrix(NA, nrow(y), length(lambda))
  if (type.measure == 'auc') {
    good = matrix(0, nfolds, length(lambda))
  }
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        
        if (keep | type.measure == 'auc') {
          nlami = sum(which_lam)
          predmat[which, seq(nlami)] = predict(fitobj, x[which, , drop = FALSE], s=lambda[which_lam], newoffset = offset[which])
        }
        
        if(type.measure == 'auc'){
          good[i, seq(nlami)] = 1
          for (j in seq(nlami)) {
            cvraw[i, j] = auc(y[which, 'status'], predmat[which, j], weights[which])
          }
        }
        else if (grouped) {
            coefmat = predict(fitobj, type = "coeff",s=lambda[which_lam])
            plfull = coxnet.deviance(x = x, y = y, offset = offset, 
                weights = weights, beta = coefmat)
            plminusk = coxnet.deviance(x = x[!which, ], y = y[!which, 
                ], offset = offset[!which], weights = weights[!which], 
                beta = coefmat)
            cvraw[i, seq(along = plfull)] = plfull - plminusk
        }
        else {
            coefmat = predict(fitobj, type = "coeff",s=lambda[which_lam])
            plk = coxnet.deviance(x = x[which, ], y = y[which, 
                ], offset = offset[which], weights = weights[which], 
                beta = coefmat)
            cvraw[i, seq(along = plk)] = plk
        }
    }
    status = y[, "status"]
    if (type.measure == 'auc') {
      N = apply(good, 2, sum)
      weights = as.vector(tapply(weights, foldid, sum))
    }
    else {
      N = nfolds - apply(is.na(cvraw), 2, sum)
      weights = as.vector(tapply(weights * status, foldid, sum))
      cvraw = cvraw/weights
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
        w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
    if (keep) 
      out$fit.preval = predmat
    out
}
