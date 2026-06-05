PanelEstimate <- function (sets, panel.data, number.iterations = 1000, df.adjustment = FALSE, 
    confidence.level = 0.95, moderator = NULL, se.method = "bootstrap", 
    pooled = FALSE, include.placebo.test = FALSE, parallel = FALSE, 
    num.cores = 1) 
{
    if (!inherits(panel.data, "PanelData")) 
        stop("Please provide a PanelData object.")
    if (pooled && !identical(se.method, "bootstrap")) {
        se.method = "bootstrap"
        warning("Pooled results only available with bootstrap SEs")
    }
    if (se.method == "wfe") 
        stop("wfe is no longer supported. Please specify se.method = 'bootstrap', 'conditional', or 'unconditional'")
    if (!is.null(moderator)) {
        if (attr(sets, "qoi") == "att") {
            s1 = sets[["att"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
        }
        if (attr(sets, "qoi") == "atc") {
            s1 = sets[["atc"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
        }
        if (attr(sets, "qoi") == "ate") {
            s1 <- sets[["att"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
        }
        if (attr(sets, "qoi") == "art") {
            s1 <- sets[["art"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
        }
        if ((attr(sets, "qoi") == "art") || (attr(sets, "qoi") == 
            "att")) {
            att.sets.in <- s1
        }
        else {
            att.sets.in <- NULL
        }
        unit.id <- attr(panel.data, "unit.id")
        time.id <- attr(panel.data, "time.id")
        treatment <- attr(panel.data, "treatment")
        set.list <- handle_moderating_variable(ordered.data = panel.data, 
            att.sets = att.sets.in, atc.sets = sets[["atc"]], 
            moderator = moderator, unit.id = unit.id, time.id = time.id, 
            PM.object = sets, qoi.in = attr(sets, "qoi"))
        res <- lapply(set.list, FUN = panel_estimate, se.method = se.method, 
            number.iterations = number.iterations, df.adjustment = df.adjustment, 
            confidence.level = confidence.level, data = panel.data, 
            pooled = pooled, include.placebo.test = include.placebo.test)
    }
    else {
        
        res = panel_estimate(se.method = se.method, number.iterations = number.iterations, 
            df.adjustment = df.adjustment, confidence.level = confidence.level, 
            sets = sets, data = panel.data, pooled = pooled, 
            include.placebo.test = include.placebo.test, parallel = parallel, 
            num.cores = num.cores)
        
    }
    return(res)
}


panel_estimate <- function (sets, data, se.method = "bootstrap", number.iterations = 1000, 
    df.adjustment = FALSE, confidence.level = 0.95, placebo.test = FALSE, 
    placebo.lead = NULL, pooled = FALSE, include.placebo.test = FALSE, 
    parallel = FALSE, num.cores = 1) {
    lead <- attr(sets, "lead")
    outcome.variable <- attr(sets, "outcome.var")
    if (include.placebo.test == TRUE) {
        old.data <- data
    }
    if (!inherits(sets, "PanelMatch")) 
        stop("sets parameter is not a PanelMatch object")
    qoi <- attr(sets, "qoi")
    if (qoi == "ate") {
        t.sets <- sets[["att"]]
        t.sets.2 <- sets[["atc"]]
        sets.x <- t.sets[sapply(t.sets, length) > 0]
        sets.y <- t.sets.2[sapply(t.sets.2, length) > 0]
        if (length(sets.x) == 0 && length(sets.y) == 0) 
            stop("do not have adequate data to proceed")
    }
    else {
        t.sets <- sets[[qoi]]
    }
    if (length(t.sets) == 0) {
        return(NA)
    }
    lag.in <- attr(t.sets, "lag")
    dependent = outcome.variable
    treatment <- attr(t.sets, "treatment.var")
    unit.id <- attr(t.sets, "id.var")
    time.id <- attr(t.sets, "t.var")
    if (identical(qoi, "ate")) {
        sets.att <- sets[["att"]]
        sets.att <- sets.att[sapply(sets.att, length) > 0]
        sets.atc <- sets[["atc"]]
        sets.atc <- sets.atc[sapply(sets.atc, length) > 0]
    }
    if (identical(qoi, "atc")) {
        sets.atc <- sets[["atc"]]
        sets.atc <- sets.atc[sapply(sets.atc, length) > 0]
        sets.att <- NULL
    }
    if (identical(qoi, "att") || identical(qoi, "art")) {
        sets.att <- sets[[qoi]]
        sets.atc <- NULL
        sets.att <- sets.att[sapply(sets.att, length) > 0]
    }
    if (qoi == "att" || qoi == "art" || qoi == "ate") {
        treated.unit.ids.att <- as.numeric(sub("\\..*", "", names(sets.att)))
        if (identical(qoi, "att") || identical(qoi, "art")) 
            treated.unit.ids.atc <- NULL
    }
    if (qoi == "atc" || qoi == "ate") {
        treated.unit.ids.atc <- as.numeric(sub("\\..*", "", names(sets.atc)))
        if (identical(qoi, "att")) 
            treated.unit.ids.att <- NULL
    }
    if (placebo.test) {
        data <- my_prepare_data(data.in = data, lead = max(placebo.lead):0, 
            sets.att = sets.att, sets.atc = sets.atc, qoi.in = qoi, 
            dependent.variable = dependent)
        pe.results <- calculate_placebo_estimates(qoi.in = qoi, 
            data.in = data, lead = lead, number.iterations = number.iterations, 
            att.treated.unit.ids = treated.unit.ids.att, atc.treated.unit.ids = treated.unit.ids.atc, 
            outcome.variable = dependent, unit.id.variable = unit.id, 
            confidence.level = confidence.level, att.sets = sets.att, 
            atc.sets = sets.atc, lag = lag.in, placebo.lead = placebo.lead, 
            se.method = se.method, parallel = parallel, num.cores = num.cores)
    }
    else {

        ## data <- data[,c("wbcode2", "year", "dem", "y", "tradewb")]
        
        data <- my_prepare_data(data.in = data, lead = lead, sets.att = sets.att, 
            sets.atc = sets.atc, qoi.in = qoi, dependent.variable = dependent)

        pe.results <- calculate_estimates(qoi.in = qoi, data.in = data, 
            lead = lead, number.iterations = number.iterations, 
            att.treated.unit.ids = treated.unit.ids.att, atc.treated.unit.ids = treated.unit.ids.atc, 
            outcome.variable = dependent, unit.id.variable = unit.id, 
            confidence.level = confidence.level, att.sets = sets.att, 
            atc.sets = sets.atc, lag = lag.in, se.method = se.method, 
            pooled = pooled, parallel = parallel, num.cores = num.cores)
        
    }
    if (include.placebo.test == TRUE) {
        pt.results <- placebo_test(sets, old.data, number.iterations = number.iterations, 
            confidence.level = confidence.level, plot = FALSE, 
            se.method = se.method, parallel = parallel, num.cores = num.cores)
        pe.results[["placebo.test"]] <- pt.results
    }
    return(pe.results)
}



calculate_estimates <- function (qoi.in, data.in, lead, number.iterations, att.treated.unit.ids, 
    atc.treated.unit.ids, outcome.variable, unit.id.variable, 
    confidence.level, att.sets, atc.sets, placebo.test = FALSE, 
    lag, se.method, pooled = FALSE, parallel = FALSE, num.cores = 1) {
    pt.estimates <- calculate_point_estimates(qoi.in, data.in, lead, outcome.variable, pooled)
    if (identical(se.method, "bootstrap")) {
        if (parallel) {
            coefs <- handle_bootstrap_parallel(qoi.in = qoi.in, 
                data.in = data.in, lead = lead, number.iterations = number.iterations, 
                att.treated.unit.ids = att.treated.unit.ids, 
                atc.treated.unit.ids = atc.treated.unit.ids, 
                outcome.variable = outcome.variable, unit.id.variable = unit.id.variable, 
                confidence.level = confidence.level, lag = lag, 
                pooled = pooled, num.cores = num.cores)
        }
        else {
            
            coefs <- handle_bootstrap(qoi.in = qoi.in, data.in = data.in, 
                lead = lead, number.iterations = number.iterations, 
                att.treated.unit.ids = att.treated.unit.ids, 
                atc.treated.unit.ids = atc.treated.unit.ids, 
                outcome.variable = outcome.variable, unit.id.variable = unit.id.variable, 
                confidence.level = confidence.level, lag = lag, 
                pooled = pooled)

            ## custom
            coefs <- my_handle_bootstrap(qoi.in = qoi.in, data.in = data.in, 
                lead = lead, number.iterations = number.iterations, 
                att.treated.unit.ids = att.treated.unit.ids, 
                atc.treated.unit.ids = atc.treated.unit.ids, 
                outcome.variable = outcome.variable, unit.id.variable = unit.id.variable, 
                confidence.level = confidence.level, lag = lag, 
                pooled = pooled)
            
        }
        if (identical(qoi.in, "att") || identical(qoi.in, "art")) {
            sets <- att.sets
        }
        else if (identical(qoi.in, "atc")) {
            sets <- atc.sets
        }
        ses <- apply(coefs, 2, sd, na.rm = TRUE)
        if (!pooled) {
            names(ses) <- paste0("t+", lead)
        }
        if (identical(qoi.in, "ate")) {
            z <- list(estimate = pt.estimates, bootstrapped.estimates = coefs, 
                bootstrap.iterations = number.iterations, standard.error = ses, 
                lead = lead, confidence.level = confidence.level, 
                qoi = "ate", matched.sets = list(att = att.sets, 
                  atc = atc.sets), se.method = se.method, pooled = pooled)
        }
        else {
            z <- list(estimate = pt.estimates, bootstrapped.estimates = coefs, 
                bootstrap.iterations = number.iterations, standard.error = ses, 
                lag = lag, lead = lead, confidence.level = confidence.level, 
                qoi = qoi.in, matched.sets = sets, se.method = se.method, 
                pooled = pooled)
        }
        class(z) <- "PanelEstimate"
        return(z)
    }
    else if (identical(se.method, "conditional")) {
        if (identical(qoi.in, "ate")) 
            stop("analytical standard errors not available for ATE")
        if (pooled) 
            stop("Analytical SE's not available for pooled estimates")
        se.s <- handle_conditional_se(qoi.in, data.in, lead, 
            outcome.variable, unit.id.variable)
        if (is.null(atc.sets)) {
            sets <- att.sets
        }
        else if (is.null(att.sets)) {
            sets <- atc.sets
        }
        else {
            stop("missing sets")
        }
        z <- list(estimate = pt.estimates, standard.error = se.s, 
            lag = lag, lead = lead, confidence.level = confidence.level, 
            qoi = qoi.in, matched.sets = sets, se.method = se.method, 
            pooled = pooled)
        class(z) <- "PanelEstimate"
        return(z)
    }
    else if (identical(se.method, "unconditional")) {
        se.s <- handle_unconditional_se(qoi.in, data.in, lead, 
            outcome.variable, unit.id.variable)
        if (is.null(atc.sets)) {
            sets <- att.sets
        }
        else if (is.null(att.sets)) {
            sets <- atc.sets
        }
        else {
            stop("missing sets")
        }
        z <- list(estimate = pt.estimates, standard.error = se.s, 
            lag = lag, lead = lead, confidence.level = confidence.level, 
            qoi = qoi.in, matched.sets = sets, se.method = se.method, 
            pooled = pooled)
        class(z) <- "PanelEstimate"
        return(z)
    }
    else {
        stop("invalid standard error method")
    }
}


my_handle_bootstrap <- function (qoi.in, data.in, lead, number.iterations, att.treated.unit.ids, 
    atc.treated.unit.ids, outcome.variable, unit.id.variable, 
    confidence.level, lag, pooled) {

    ## browser()
    
    coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
    if (qoi.in %in% c("att", "art", "atc")) {

        per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]), 
            FUN = perunitSum, lead.in = lead, dependent.in = outcome.variable, 
            qoi_in = qoi.in)

        per.unit.dit.sums <- by(data.in, as.factor(data.in[, 
                                                           unit.id.variable]), FUN = perunitSum_Dit, qoi_in = qoi.in)
        
        units.id <- names(per.unit.sums)

        ## perunitSum_Dit(data.in[data.in$wbcode2 == 95,], qoi.in) 

        ## perunitSum(data.in[data.in$wbcode2 == 95,], lead.in = lead, dependent.in = outcome.variable, qoi_in = qoi.in) 

        ## per.unit.dit.sums
        
        tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
        colnames(tdf) <- sapply(lead, function(x) paste0("Wit_", 
            qoi.in, x))
        tdf$unit.id <- NA
        tdf$unit.id <- as.character(unlist(units.id))
        tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
        for (k in 1:number.iterations) {
            clusters <- unique(data.in[, unit.id.variable])
            units <- sample(clusters, size = length(clusters), 
                replace = TRUE)
            if (identical(qoi.in, "att") || identical(qoi.in, 
                "art")) {
                treated.unit.ids <- att.treated.unit.ids
            }
            else {
                treated.unit.ids <- atc.treated.unit.ids
            }
            while (all(!units %in% treated.unit.ids)) {
                units <- sample(clusters, size = length(clusters), 
                  replace = TRUE)
            }

            fdf <- as.data.frame(table(units))
            colnames(fdf) <- c("unit.id", "frequency")
            bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
            bdf$frequency[is.na(bdf$frequency)] <- 0
            if (identical(qoi.in, "atc")) {
                multiplier <- -1
            }
            else {
                multiplier <- 1
            }

            dt_bdf <- adt(bdf)[, status := ifelse(unit.id %in% treated.unit.ids, "treat", "control")]
            dt_bdf[status == "control" & frequency != 0]
            


            coefs_sum <- colSums(bdf[, sapply(lead, function(x) paste0("Wit_", qoi.in, x)), drop = FALSE] * # selecting
                    bdf[, "frequency"], # multiplying with bootstrap frequency
                    na.rm = FALSE)

            ## withAutoprint(coefs_sum)

            # first column sums, then divide by something
            coefs[k, ] <- multiplier * coefs_sum/ 
                (sum(bdf[, "frequency"] * bdf[, "Dit"]))
        }
        if (pooled) {
            coefs <- as.matrix(rowMeans(coefs), ncol = 1)
        }
        return(coefs)
    }
    else if (identical(qoi.in, "ate")) {
        qoi.in <- "att"
        per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]), 
            FUN = perunitSum, lead.in = lead, dependent.in = outcome.variable, 
            qoi_in = qoi.in)
        per.unit.dit.sums <- by(data.in, as.factor(data.in[, 
            unit.id.variable]), FUN = perunitSum_Dit, qoi_in = qoi.in)
        units.id <- names(per.unit.sums)
        tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
        colnames(tdf) <- sapply(lead, function(x) paste0("Wit_", 
            "att", x))
        tdf$unit.id <- NA
        tdf$unit.id <- as.character(unlist(units.id))
        tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
        qoi.in <- "atc"
        per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]), 
            FUN = perunitSum, lead.in = lead, dependent.in = outcome.variable, 
            qoi_in = qoi.in)
        per.unit.dit.sums <- by(data.in, as.factor(data.in[, 
            unit.id.variable]), FUN = perunitSum_Dit, qoi_in = qoi.in)
        units.id <- names(per.unit.sums)
        tdf.atc <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
        colnames(tdf.atc) <- sapply(lead, function(x) paste0("Wit_", 
            "atc", x))
        tdf.atc$unit.id <- NA
        tdf.atc$unit.id <- as.character(unlist(units.id))
        tdf.atc$Dit <- as.numeric(unlist(per.unit.dit.sums))
        for (k in 1:number.iterations) {
            clusters <- unique(data.in[, unit.id.variable])
            units <- sample(clusters, size = length(clusters), 
                replace = TRUE)
            while (all(!units %in% att.treated.unit.ids) || all(!units %in% 
                atc.treated.unit.ids)) {
                units <- sample(clusters, size = length(clusters), 
                  replace = TRUE)
            }
            fdf <- as.data.frame(table(units))
            colnames(fdf) <- c("unit.id", "frequency")
            bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
            bdf$frequency[is.na(bdf$frequency)] <- 0
            att_new <- colSums(bdf[, sapply(lead, function(x) paste0("Wit_att", 
                x)), drop = FALSE] * bdf[, "frequency"], na.rm = FALSE)/(sum(bdf[, 
                "frequency"] * bdf[, "Dit"]))
            bdf.atc <- merge(x = tdf.atc, y = fdf, all.x = TRUE)
            bdf.atc$frequency[is.na(bdf.atc$frequency)] <- 0
            atc_new <- -colSums(bdf.atc[, sapply(lead, function(x) paste0("Wit_atc", 
                x)), drop = FALSE] * bdf.atc[, "frequency"], 
                na.rm = FALSE)/(sum(bdf.atc[, "frequency"] * 
                bdf.atc[, "Dit"]))
            coefs[k, ] <- (att_new * sum(bdf$Dit * bdf$frequency) + 
                atc_new * sum(bdf.atc$Dit * bdf.atc$frequency))/(sum(bdf$Dit * 
                bdf$frequency) + sum(bdf.atc$Dit * bdf.atc$frequency))
        }
        if (pooled) {
            coefs <- as.matrix(rowMeans(coefs), ncol = 1)
        }
        return(coefs)
    }
    else {
        stop("qoi misspecified")
    }
}


my_prepare_data <- function (data.in, lead, sets.att = NULL, sets.atc = NULL, qoi.in, 
                             dependent.variable) {

    browser()

    if (identical(qoi.in, "att") || identical(qoi.in, "art") || 
        identical(qoi.in, "ate")) {
        if (!identical(qoi.in, "ate")) 
            qoi.t <- qoi.in
        if (identical(qoi.in, "ate")) 
            qoi.t <- "att"
        for (j in lead) {
            dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.att)
            data.in = merge(x = data.in, y = dense.wits, all.x = TRUE, 
                by.x = colnames(data.in)[1:2], by.y = c("id", 
                  "t"))
            colnames(data.in)[length(data.in)] <- paste0("Wit_", 
                qoi.t, j)
            data.in[is.na(data.in[, length(data.in)]), length(data.in)] <- 0
        }
        data.in[, paste0("dit_", qoi.t)] <- getDits(matched_sets = sets.att, 
            data = data.in)
        colnames(data.in)[length(data.in)] <- paste0("dits_", 
            qoi.t)
        data.in[, paste0("Wit_", qoi.t, "-1")] <- 0
    }
    if (identical(qoi.in, "atc") || identical(qoi.in, "ate")) {
        for (j in lead) {
            dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.atc)
            data.in = merge(x = data.in, y = dense.wits, all.x = TRUE, 
                by.x = colnames(data.in)[1:2], by.y = c("id", 
                  "t"))
            colnames(data.in)[length(data.in)] <- paste0("Wit_atc", 
                j)
            data.in[is.na(data.in[, length(data.in)]), length(data.in)] <- 0
        }
        data.in$dit_atc <- getDits(matched_sets = sets.atc, data = data.in)
        colnames(data.in)[length(data.in)] <- "dits_atc"
        data.in$`Wit_atc-1` <- 0
    }
    data.in[, dependent.variable][is.na(data.in[, dependent.variable])] <- 0
    return(data.in)
}



calculate_point_estimates <- function (qoi.in, data.in, lead, outcome.variable, pooled = FALSE) {

    browser()
    if (identical(qoi.in, "att") || identical(qoi.in, "atc") || 
        identical(qoi.in, "art")) {
        col.idx <- sapply(lead, function(x) paste0("Wit_", qoi.in, 
            x))
        x.in <- data.in[, col.idx, drop = FALSE]
        y.in <- data.in[c(outcome.variable)][, 1]
        z.in <- data.in[, paste0("dits_", qoi.in)]
        o.coefs <- sapply(x.in, equality_four, y = y.in, z = z.in)
        if (identical(qoi.in, "atc")) 
            o.coefs <- -o.coefs
        if (all(lead >= 0)) {
            if (length(lead[lead < 0]) > 1) {
                names(o.coefs)[(length(o.coefs) - max(lead[lead >= 
                  0])):length(o.coefs)] <- sapply(lead[lead >= 
                  0], function(x) paste0("t+", x))
                names(o.coefs)[(length(o.coefs) - length(lead) + 
                  1):length(lead[lead < 0])] <- sapply(lead[lead < 
                  0], function(x) paste0("t", x))
            }
            else {
                names(o.coefs) <- sapply(lead, function(x) paste0("t+", 
                  x))
            }
        }
        if (pooled) {
            o.coefs <- mean(o.coefs, na.rm = TRUE)
            names(o.coefs) <- NULL
        }
    }
    else if (identical(qoi.in, "ate")) {
        o.coefs_att <- sapply(data.in[, sapply(lead, function(x) paste0("Wit_att", 
            x)), drop = FALSE], equality_four, y = data.in[c(outcome.variable)][, 
            1], z = data.in$dits_att)
        o.coefs_atc <- -sapply(data.in[, sapply(lead, function(x) paste0("Wit_atc", 
            x)), drop = FALSE], equality_four, y = data.in[c(outcome.variable)][, 
            1], z = data.in$dits_atc)
        o.coefs <- (o.coefs_att * sum(data.in$dits_att) + o.coefs_atc * 
            sum(data.in$dits_atc))/(sum(data.in$dits_att) + sum(data.in$dits_atc))
        if (length(lead[lead < 0]) > 1) {
            names(o.coefs)[(length(o.coefs) - max(lead[lead >= 
                0])):length(o.coefs)] <- sapply(lead[lead >= 
                0], function(x) paste0("t+", x))
            names(o.coefs)[(length(o.coefs) - length(lead) + 
                1):length(lead[lead < 0])] <- sapply(lead[lead < 
                0], function(x) paste0("t", x))
        }
        else {
            names(o.coefs) <- sapply(lead, function(x) paste0("t+", 
                x))
        }
        if (pooled) {
            o.coefs <- mean(o.coefs, na.rm = TRUE)
            names(o.coefs) <- NULL
        }
    }
    return(o.coefs)
}
