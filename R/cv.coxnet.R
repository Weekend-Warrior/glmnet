cv.coxnet <-
function (outlist, lambda, x, y, weights, offset, foldid, type.measure, 
    grouped, keep = FALSE) 
{
    typenames = c(deviance = "Partial Likelihood Deviance")
    if (type.measure == "default") 
        type.measure = "deviance"
    if (!match(type.measure, c("deviance"), FALSE)) {
        warning("Only 'deviance'  available for Cox models; changed to type.measure='deviance'")
        type.measure = "deviance"
    }
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    nfolds = max(foldid)
    if ((length(weights)/nfolds < 10) && !grouped) {
        warning("Option grouped=TRUE enforced for cv.coxnet, since < 3 observations per fold", 
            call. = FALSE)
        grouped = TRUE
    }
    cvraw = matrix(NA, nfolds, length(lambda))
###We dont want to extrapolate lambdas on the small side
  mlami=max(sapply(outlist,function(obj)min(obj$lambda)))
  which_lam=lambda >= mlami
  if (keep) predmat = matrix(NA, nrow(y), length(lambda))
 
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        coefmat = predict(fitobj, type = "coeff",s=lambda[which_lam])
        if (grouped) {
            plfull = coxnet.deviance(x = x, y = y, offset = offset, 
                weights = weights, beta = coefmat)
            plminusk = coxnet.deviance(x = x[!which, ], y = y[!which, 
                ], offset = offset[!which], weights = weights[!which], 
                beta = coefmat)
            cvraw[i, seq(along = plfull)] = plfull - plminusk
        }
        else {
            plk = coxnet.deviance(x = x[which, ], y = y[which, 
                ], offset = offset[which], weights = weights[which], 
                beta = coefmat)
            cvraw[i, seq(along = plk)] = plk
        }
        if (keep) {
          preds = predict(fitobj, x[which, , drop = FALSE], s=lambda[which_lam], newoffset = offset[which])
          nlami = sum(which_lam)
          predmat[which, seq(nlami)] = preds
        }
    }
    status = y[, "status"]
    N = nfolds - apply(is.na(cvraw), 2, sum)
    weights = as.vector(tapply(weights * status, foldid, sum))
    cvraw = cvraw/weights
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
        w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
    if (keep) 
      out$fit.preval = predmat
    out
}
