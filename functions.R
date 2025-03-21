#' Draw parameters for a specific age group from the parameter inputs table -----------
#' 
#' @param nsamples integer value specifying the number of samples to draw.
#' @param pars.table data.frame of input values for each parameter distribution.
#'
#' @return
#' @export
#'
#' @examples
get_parameters <- function(nsamples, pars.table) {
    
    lapply(1:nrow(pars.table), function(x){  
        data.table(
            # fixed
            sim  = 1:nsamples,
            risk = pars.table$risk[x],
            age  = pars.table$age[x],
            population = pars.table$population[x],
            proportion = pars.table$proportion[x],
            
            # CAR, hospitalization risk
            clinical.attack.rate = runif(n = nsamples, min = pars.table$car.lo[x],       max = pars.table$car.hi[x]),
            hosp.risk            = runif(n = nsamples, min = pars.table$hosp.risk.lo[x], max = pars.table$hosp.risk.hi[x]), 
            
            # care-seeking
            seek.care       = runif(n = nsamples, min = pars.table$seek.care.lo[x],  max = pars.table$seek.care.hi[x]), 
            seek.care.early = runif(n = nsamples, min = pars.table$seek.early.lo[x], max = pars.table$seek.early.hi[x]),
            
            # test sensitivity
            test.sens       = runif(n = nsamples, min = pars.table$test.sens.lo[x],      max = pars.table$test.sens.hi[x]),
            
            # antiviral compliance
            nai.compliance  = runif(n = nsamples, min = pars.table$nai.compliance.lo[x], max = pars.table$nai.compliance.hi[x]) 
        )
    } ) %>% bind_rows()
}

#' Get model parameters stratified by low/high risk ------------------------
#' 
#' If the population CHR is a weighted average of the high/low-risk CHRS, then
#' CHR_pop = prop_HR x RR_HR:LR x CHR_LR + prop_LR x CHR_LR
#' and CHR_LR = CHR_pop / (prop_HR x RR + prop_LR)
#'
#' @param pop numeric vector specifying population averages (length = number of age groups)
#' @param rr relative risk or risk ratio for high-risk vs low-risk parameters
#' @param prop.high numeric vector with prevalence of high-risk conditions among each age group
#' 
#' @return
#' @export
#'
#' @examples
get_risk_stratifications <- function(pop, rr, low, high, prop.high) {
    
    low.vals  <- pop / (prop.high * rr + (1 - prop.high))
    high.vals <- low.vals * rr
    
    return(list(low = low.vals, high = high.vals))
}

    
#' Simulate outcomes from the model  ------------------------
#' 
#' @param parameters data.frame or data.table of samples parameter values
#' @param only.highrisk logical TRUE/FALSE indicating whether to only include high risk individuals. Default = FALSE
#'
#' @return
#' @export
#'
#' @examples
get_outcomes <- function(parameters) {
    parameters %>% 
        mutate(
            ## Partition into branches ------------------
            # those with flu symptoms to include in model (assign to risk-category according to assumed proportions)
            flu.sympts = population * clinical.attack.rate * proportion,
            
            # those who do/don't seek care
            outpatient     = flu.sympts * seek.care,
            non.outpatient = flu.sympts * (1 - seek.care),
            
            # those who seek care early/late
            outpatient.early = seek.care.early * outpatient,
            outpatient.late  = (1 - seek.care.early) * outpatient,
            
            # outpatients who do/don't receive a positive test
            pos.test.early = perc.tested.early * test.sens * outpatient.early,
            pos.test.late  = perc.tested.late  * test.sens * outpatient.late,
            
            neg.test.early = perc.tested.early * (1 - test.sens) * outpatient.early,
            neg.test.late  = perc.tested.late  * (1 - test.sens) * outpatient.late,
            
            no.test.early = (1 - perc.tested.early) * outpatient.early,
            no.test.late  = (1 - perc.tested.late) * outpatient.late,
            
            # positive test outpatients prescribed / not prescribed NAIs
            nai.pos.early = prescribed.nai.pos.early * pos.test.early,
            nai.pos.late  = prescribed.nai.pos.late  * pos.test.late,
            
            non.nai.pos.early = (1 - prescribed.nai.pos.early) * pos.test.early,
            non.nai.pos.late  = (1 - prescribed.nai.pos.late)  * pos.test.late,
            
            # negative test outpatients prescribed / not prescribed NAIs
            nai.neg.early = prescribed.nai.neg.early * neg.test.early,
            nai.neg.late  = prescribed.nai.neg.late  * neg.test.late,
            
            non.nai.neg.early = (1 - prescribed.nai.neg.early) * neg.test.early,
            non.nai.neg.late  = (1 - prescribed.nai.neg.late)  * neg.test.late,
            
            # no test outpatients prescribed / not prescribed NAIs
            nai.no.early = prescribed.nai.no.early * no.test.early,
            nai.no.late  = prescribed.nai.no.late  * no.test.late,
            
            non.nai.no.early = (1 - prescribed.nai.no.early) * no.test.early,
            non.nai.no.late  = (1 - prescribed.nai.no.late)  * no.test.late,
            
            # positive test outpatients prescribed NAIs who do/don't comply
            comply.nai.pos.early = nai.compliance * nai.pos.early,
            comply.nai.pos.late  = nai.compliance * nai.pos.late,
            
            non.comply.nai.pos.early = (1 - nai.compliance) * nai.pos.early,
            non.comply.nai.pos.late  = (1 - nai.compliance) * nai.pos.late,
            
            # negative test outpatients prescribed NAIs who do/don't comply
            comply.nai.neg.early = nai.compliance * nai.neg.early,
            comply.nai.neg.late  = nai.compliance * nai.neg.late,
            
            non.comply.nai.neg.early = (1 - nai.compliance) * nai.neg.early,
            non.comply.nai.neg.late  = (1 - nai.compliance) * nai.neg.late,
            
            # no test outpatients prescribed NAIs who do/don't comply
            comply.nai.no.early = nai.compliance * nai.no.early,
            comply.nai.no.late  = nai.compliance * nai.no.late,
            
            non.comply.nai.no.early = (1 - nai.compliance) * nai.no.early,
            non.comply.nai.no.late  = (1 - nai.compliance) * nai.no.late,
            
            # Hospitalizations from different compartments -------------------
            # Non-outpatients (those who didn't seek care)
            hosp.non.outpatient = hosp.risk * non.outpatient,
            
            # Outpatients who weren't prescribed NAIs (positive, negative, no test / early, late care-seeking)
            hosp.non.prescribed.nai = hosp.risk * (non.nai.pos.early + non.nai.pos.late + non.nai.neg.early + non.nai.neg.late + non.nai.no.early + non.nai.no.late),
            
            # Prescribed NAIs who weren't compliant
            hosp.non.comply.nai = hosp.risk * (non.comply.nai.pos.early + non.comply.nai.pos.late + non.comply.nai.neg.early + non.comply.nai.neg.late + non.comply.nai.no.early + non.comply.nai.no.late),
            
            # Prescribed NAIs who did comply
            hosp.comply.nai = hosp.risk * ( (1 - nai.effect.early) * (comply.nai.pos.early + comply.nai.neg.early + comply.nai.no.early) + (1 - nai.effect.late) * (comply.nai.pos.late + comply.nai.neg.late + comply.nai.no.late) ),
            
            # All hospitalizations
            hosp.total = hosp.non.outpatient + hosp.non.prescribed.nai + hosp.non.comply.nai + hosp.comply.nai,
            
            # Number of antiviral prescriptions filled -------------------
            nai.total = nai.pos.early + nai.pos.late + nai.neg.early + nai.neg.late + nai.no.early + nai.no.late
        ) 
}



#' Summarize (mean and 95th percentile) function -------------------------------
#' 
#' @param dat data.frame
#' @param var Character specifying variable to be summarized
#' @param groups Character string specifying columns to group by
#'
#' @return
#' @export
#'
#' @examples
get_summary <- function(dat, var, groups, level = 0.95) {
    dat %>% group_by(!!!syms(groups)) %>% 
        summarize(mean = mean(get(var)), 
                  lower = quantile(get(var), (1 - 0.95)/2), 
                  upper = quantile(get(var), 1 - (1 - 0.95)/2), .groups = "drop")
}


#' Plotting function -----------------------------------------
#' 
#' @param txt Numerical value specifying the plot text size. Default value is 12.
#' @param ... Other optional arguments passed to theme()
#'
#' @return
#' @export
#'
#' @examples
get_theme <- function(txt = 10, ...) {
    theme_light() + theme(axis.text = element_text(size = txt),
                          axis.title = element_text(size = txt),
                          legend.text = element_text(size = txt),
                          legend.title = element_text(size = txt),
                          strip.text = element_text(size = txt - 1, color = "black"),
                          title = element_text(size = txt - 1),
                          strip.background = element_rect(fill = "white", color = "black"),
                          panel.grid.minor = element_blank(),
                          panel.grid.major.x = element_blank(),
                          ...)
}



#' Renaming parameters for plotting ---------------------------
#' 
#' @param dat data.frame or data.table of results to plot
#' @param reorder logical TRUE/FALSE indicating whether the new names should appear in a specific order. Default value is TRUE.
#'
#' @return
#' @export
#'
#' @examples
rename_pars <- function(dat, reorder = TRUE) {
    dat <- dat %>% mutate(parname = recode(param, 
                                    high.risk = "Fraction high-risk",
                                    clinical.attack.rate = "Clinical attack rate",
                                    hosp.risk = "Risk of hospitalization",
                                    seek.care = "Seek care", 
                                    seek.care.early = "Seek care within 48 hours", 
                                    
                                    perc.tested.early = "Fraction tested (< 48hrs)", 
                                    perc.tested.late = "Fraction tested (> 48hrs)", 
                                    test.sens = "Test sensitivity",
                                    
                                    diagnosis.early = "Positive test (< 48hrs)",
                                    diagnosis.late  = "Positive test (> 48hrs)",
                                    
                                    prescribed.nai.pos.early = "Prescribed NAIs (< 48hrs, pos test)", 
                                    prescribed.nai.pos.late  = "Prescribed NAIs (> 48hrs, pos test)", 
                                    prescribed.nai.neg.early = "Prescribed NAIs (< 48hrs, no pos test)", 
                                    prescribed.nai.neg.late  = "Prescribed NAIs (> 48hrs, no pos test)", 
                                    
                                    nai.compliance = "Full compliance", 
                                    nai.effect.early = "NAI effect (< 48hrs)", 
                                    nai.effect.late = "NAI effect (> 48hrs)"
                                    ) )
    
    if(reorder) dat$parname <- factor(dat$parname, levels = c("Fraction high-risk", "Clinical attack rate",  "Risk of hospitalization",
                                                              "Seek care", "Seek care within 48 hours", 
                                                              "Fraction tested (< 48hrs)", "Fraction tested (> 48hrs)", "Test sensitivity", 
                                                              "Positive test (< 48hrs)", "Positive test (> 48hrs)",
                                                              "Prescribed NAIs (< 48hrs, pos test)", "Prescribed NAIs (> 48hrs, pos test)",
                                                              "Prescribed NAIs (< 48hrs, no pos test)", "Prescribed NAIs (> 48hrs, no pos test)", 
                                                              "Full compliance",  "NAI effect (< 48hrs)", "NAI effect (> 48hrs)"))
    return(dat)
}




