#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2021
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "IBM SPSS, JKP"
# version__ = "1.0.1"

# History
# 06-may-2015 Original Version
# 31-oct-2021 fixes for strata
# 13-jul-2025 add CI's to output and implement confidence intervals

### MAIN ROUTINE ###
docomprisk = function(ftime, fstatus, cov1=NULL, cov2=NULL, tf="quad", 
    threshold=NULL, cengroup=NULL, failcode, cencode,
    maxiter=10, gtol=1e-06, missing="omit",
    dataset=NULL, plotit=TRUE, confint=0.95,
    quantiles = c(.25, .50, .75)) {
    # Estimate competing risk model

    setuplocalization("STATS_COMPRISK")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Competing Risks Regression")
    warningsprocname = gtxt("Competing Risks Regression: Warnings")
    omsid="STATSCOMPRISK"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(cmprsk), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "cmprsk"),dostop=TRUE)
        }
    )

    if (!is.null(dataset)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
                dostop=TRUE)
        }
        if (dataset %in% alldatasets) {
            warns$warn(gtxt("The output dataset name must not already be in use"),
                dostop=TRUE)
        }
    }

    if (is.null(cov1) && is.null(cov2)) {
        warns$warn(gtxt("At least one fixed or time-varying covariate is required"), dostop=TRUE)
    }
    
    alldata = c(ftime, fstatus, cov1, cov2, cengroup)
    nacode = ifelse(missing == "omit", na.omit, na.fail)
    allargs = as.list(environment())
    dta = tryCatch(spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE,
        factorMode="levels"),
        error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )
    ncases = nrow(dta)
    dta = dta[complete.cases(dta),]
    nmissing = ncases - nrow(dta)
    if (nmissing > 0 && identical(nacode, na.fail)) {
        warns$warn(gtxtf("Stopping: there are %d missing cases and Stop was specified", nmissing),
            dostop=TRUE)
    }
    allargs['nmissing'] = nmissing
    if (!is.null(cov1)) {
        cov1dta = data.frame(dta[,3:(2 + length(cov1))])
        names(cov1dta) = names(dta)[3:(2 + length(cov1))]
        # convert categorical covariates to dummies and build new data frame
        cov1dta = handlecatcovar(cov1dta)
    } else {
        cov1dta=NULL
    }
    if (!is.null(cov2)) {
        cov2dta = data.frame(dta[,(3 + length(cov1)):(2 + length(cov1) + length(cov2))])

        names(cov2dta) = names(dta)[(3 + length(cov1)):(2 + length(cov1) + length(cov2))]
        cov2dta = handlecatcovar(cov2dta)
        if (is.null(tf)) {
            warns$warn(gtxt("A time function is required if using time covariates"),
                dostop=TRUE)
        }
        assign("cov2dta", cov2dta, env=.GlobalEnv)
        if (tf == "quad") {
            tffunc = tfquad
        } else {
            tffunc = tfthresh
            if (is.null(threshold)) {
                warns$warn(gtxt("A threshold value is required if the threshold time function is used"), 
                                dostop=TRUE)
            }
            assign("threshold", threshold, envir=.GlobalEnv)
        }
    } else {
        cov2dta = NULL
    }

    nacode = ifelse(missing == "omit", na.omit, na.fail)
    args = list(ftime=dta[ftime], fstatus=dta[fstatus], failcode=failcode,
        cencode=cencode, na.action=nacode, gtol=gtol, maxiter=maxiter)
    ###print(args)
    # crr seems to have hardwired these names.  Renaming is required
    # or crr fails :-(
    names(args$ftime) = "ftime"    
    names(args$fstatus) = "fstatus"
    if (!is.null(cengroup)) {
        args$cengroup = dta[cengroup]
        names(args$cengroup) = "cengroup"
    }
    pargs = list()
    if (!is.null(cov1)) {
        args[["cov1"]] = cov1dta
        pargs[["cov1"]] = cov1dta
    }
    if (!is.null(cov2)) {
        args[["cov2"]] = cov2dta
        pargs[["cov2"]] = cov2dta
        args[["tf"]] = tffunc
    }
    ###print(args)
    res = tryCatch(do.call(crr, args),
        error = function(e) {
            warns$warn(e$message, dostop=TRUE)
        }
    )
    if (plotit || !is.null(dataset)) {
        doquantile = function(x, q=quantiles) {
            return(quantile(x, q, na.rm=TRUE))
        }

        pts1 = as.matrix(sapply(cov1dta, doquantile))
        pts2 = as.matrix(sapply(cov2dta, doquantile))
        if (nrow(pts1) > 0 && nrow(pts2) > 0) {
            pred = predict(res, cov1=pts1, cov2=pts2)
            evalpts = data.frame(pts1, pts2)
            names(evalpts) = names(res$coef)
        } else {
            if (nrow(pts1) > 0) {
                pred = predict(res, cov1=pts1)
                evalpts = data.frame(pts1)
            } else {
                pred = predict(res, cov1=pts2, cov2=pts2)
                evalpts = data.frame(pts2)
            }
        }
    } else {
        pred = NULL
        evalpts = NULL
    }
    # suspect that this will have the common can't find data problem :-(
    displayresults(allargs, res, pred, evalpts, confint, warns)

    if (!is.null(dataset)) {
        savepred(allargs, res, warns)
    }

}

tfquad = function(ft) {
    # quadratic function of time
    # It will be called only once for all the covariates
    # ft is the vector of times for value 1
    # in cov2 (TIMECOVAR), but the number of covariates is not passed
    # so we have to beg docmprsk for it :-()

    cols = ncol(get("cov2dta", envir=.GlobalEnv))
    ftx = rep(ft^2, cols)
    dim(ftx) = c(length(ft), cols)
    return(ftx)
}

tfthresh = function(ft) {
    # threshold function as above
    thresh = get("threshold", envir=.GlobalEnv)
    cols = ncol(get("cov2dta", envir=.GlobalEnv))
    ftx = rep(ifelse(ft >=thresh, 1, 0), cols)
    dim(ftx) = c(length(ft), cols)
    return(ftx)
}
    
handlecatcovar = function(df) {
    #convert any factors in df to dummies and return new data frame
    
    # get names of any factors
    catcovar = names(df)[sapply(df, is.factor)]
    if (length(catcovar) == 0) {
        return(df)
    }
    frml = as.formula(paste("~", paste(catcovar, collapse="+"), collapse=""))
    catdf = data.frame(model.matrix(frml, data=df))[-1]  # converted categorical df
    noncatcovar=names(df)[sapply(df, function(f) !is.factor(f))]
    if (length(noncatcovar) == 0) {
        return(catdf)
    }
    return(data.frame(df[noncatcovar], catdf))
}

displayresults = function(allargs, res, pred, evalpts, confint, warns) {
    # display results
    # allargs is the parameter set (estimation or prediction)
    ###print(pred)
    
    ressum = summary(res, conf.int=confint)
    tflabels = list(quad=gtxt("Quadratic"), threshold=gtxt("Threshold: %s"))
    if (!is.null(allargs$tf)) {
        tflabel = tflabels[[allargs$tf]]
        if (allargs$tf == "threshold") {
            tflabel = sprintf(tflabel, allargs$threshold)
        }
    }
    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    # summary results
    # input specifications
    # although groups can be specified (cengroup), separate results are not
    # produced.
    lbls = c(gtxt("Time Variable"),
             gtxt("Status Variable"),
             gtxt("Time Function"),
             gtxt("Group Variable"),
             gtxt("Failure Code"),
             gtxt("Censoring Code"),
             gtxt("Missing Value Treatment"),
             gtxt("Convergence"),
             gtxt("Maximum Number of Iterations"),
             gtxt("Number of Valid Cases"),
             gtxt("Number of Missing Cases"),
             gtxt("Pseudo Log Likelihood"),
             gtxt("Pseudo Likelihood Ratio Test"),
             gtxt("Log Likelihood D. F."),
             gtxt("Prediction Dataset")
    )

    vals = c(
            allargs$ftime,
            allargs$fstatus,
            ifelse(is.null(allargs$tf), gtxt("--NA--"), tflabel),
            ifelse(is.null(allargs$cengroup), gtxt("--NA--"), allargs$cengroup),
            allargs$failcode,
            allargs$cencode,
            ifelse(allargs$missing == "omit", gtxt("omit"), gtxt("Stop")),
            ifelse(res$converged, gtxt("Yes"), gtxt("NO")),
            allargs$maxiter,
            ressum$n - ressum$n.missing,
            ressum$n.missing,
            round(res$loglik, 5),
            round(ressum$logtest[[1]], 5),
            ressum$logtest[[2]],
            ifelse(is.null(allargs$dataset), gtxt("--NA--"), allargs$dataset)
    )

    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="COMPRESKSUMMARY", outline=gtxt("Summary"),
        caption = gtxt("Computations done by R package cmprsk by Bob Gray")
    )
    
    sumres = summary(res, conf.int=confint)
    coefs = data.frame(sumres$coef)
    sumnames = names(sumres$conf.int)
    cis = data.frame(sumres$conf.int)[3:4]
    cin = names(cis)
    cin1 = substr(cin[[1]], 2, nchar(cin[[1]]))
    cin2 = substr(cin[[2]], 2, nchar(cin[[2]]))
    ###save(sumres, coefs, sumnames, cis, cin, cin1, cin2, file="c:/comprisk/fromeh/summary.rdata")
    names(cis) = list(paste(gtxt("CI"), cin1),
        paste(gtxt("CI"), cin2)
    )
    names(coefs) = c(
        gtxt("Coefficient"), gtxt("Exp"), gtxt("Std. Error"), gtxt("Z"), gtxt("Sig."))
    coefs = cbind(coefs, cis)

    spsspivottable.Display(coefs,
        title=gtxt("Coefficients"),
        rowdim=gtxt("Variable"), hiderowdimtitle=FALSE,
        templateName="COMPRESCOEF",
        outline=gtxt("Coefficients"),
        caption=gtxtf("Failure code: %s", allargs$cencode))

    if (allargs$plotit) {
        plot(pred, main=gtxt("Subdistribution Functions"), xlab=gtxt("Time"),ylab="")
        grid()  # a light background grid
        
        # label quantiles with plot line style
        lty=c(gtxt("solid"), gtxt("dashed"), gtxt("dotted"), gtxt("dot dash"), 
            gtxt("long dash"), gtxt("two dash"))
        rr = row.names(evalpts)
        for (i in 1:length(rr)) {
            rr[i] = paste(rr[i], lty[(i-1) %% 6 + 1], sep=" - ")
        }
        row.names(evalpts) = rr
        spsspivottable.Display(t(evalpts), 
            title=gtxt("Covariate Evaluation Points"),
            rowdim=gtxt("Variable"), hiderowdimtitle=FALSE,
            coldim=gtxt("Quantile"), hidecoldimtitle=FALSE,
            templateName="COMPRESEVAL",
            outline=gtxt("Covariate Evaluation Points"))
    }
    
    spsspkg.EndProcedure()
}

savepred = function(allargs, res, warns) {
    # save residuals
    
    dict = list()
    resdf = data.frame(res$res)
    nam = names(resdf)
    for (n in 1:ncol(resdf)) {
        dict[[n]] = c(nam[n], "", 0, "F8.2", "scale")
    }
    dict = spssdictionary.CreateSPSSDictionary(dict)
    spssdictionary.SetDictionaryToSPSS(allargs$dataset, dict)
    tryCatch(spssdata.SetDataToSPSS(allargs$dataset, resdf),
        error=function(e) {warns$warn(e$message, dostop=TRUE)}
    )
    spssdictionary.EndDataStep()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_COMPRISK"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_COMPRISK"))
}


Run = function(args) {
    #Execute the STATS COMPRISK command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("FAILTIME", subc="", ktype="existingvarlist", var="ftime"),
        spsspkg.Template("FAILSTATUS", subc="", ktype="existingvarlist", var="fstatus"),
        spsspkg.Template("FIXEDCOVAR", subc="", ktype="existingvarlist", var="cov1", islist=TRUE),
        spsspkg.Template("TIMECOVAR", subc="", ktype="existingvarlist", var="cov2", islist=TRUE),
        spsspkg.Template("TIMEFUNC", subc="", ktype="str", var="tf",
            vallist=list("quad", "threshold")),
        spsspkg.Template("THRESHOLD", subc="", ktype="float", var="threshold"),
        spsspkg.Template("STRATA", subc="", ktype="existingvarlist", var="cengroup"),
        spsspkg.Template("FAILCODE", subc="", ktype="str", var="failcode"),
        spsspkg.Template("CENSORCODE", subc="", ktype="str", var="cencode"),
        
        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", var="maxiter"),
        spsspkg.Template("TOL", subc="OPTIONS", ktype="float", var="gtol",
            vallist=list(1e-10)),
        spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missing",
            vallist=list("omit", "fail")),
        spsspkg.Template("PLOT", subc="OPTIONS", ktype="bool", var="plotit"),
        spsspkg.Template("QUANTILES", subc="OPTIONS", ktype="float", var="quantiles", 
            islist=TRUE, vallist=list(0, 1)),
        spsspkg.Template("CONFINT", subc="OPTIONS", ktype="float", var="confint",
            vallist=list(0.0001, .9999)),
        
        spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset")      
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "docomprisk")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
