# Incremental analysis of different monitoring strategies
# Varying frequency and age group

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(BCEA)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Functions ----

# Function to plot boxplot whiskers as 95% percentile
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
# Function for cost-efficiency frontier adapted from BCEA package
# Dominance is assessed on median instead of the mean as in original function
ceef.plot.median <- function (he, comparators = NULL, pos = c(1, 1), start.from.origins = TRUE,
                              threshold = NULL, flip = FALSE, dominance = TRUE, relative = FALSE,
                              print.summary = TRUE, graph = c("base", "ggplot2"), return_dominated_only = FALSE,
                              ...) {
  if (is.null(he$c) | is.null(he$e))
    stop("Please use the bcea() function from BCEA version >=2.1-0 or attach the vectors e and c to the bcea object. Please see ?ceef.plot for additional details.")
  if (is.null(threshold))
    threshold <- base::pi/2
  else {
    if (threshold <= 0) {
      warning("The value of the cost-effectiveness threshold should be positive. The argument will be ignored.")
      threshold <- base::pi/2
    }
    else threshold <- atan(threshold)
  }
  exArgs = list(...)
  if (exists("print.plot", exArgs)) {
    print.plot = exArgs$print.plot
  }
  else {
    print.plot = TRUE
  }
  if (!is.null(comparators)) {
    stopifnot(all(comparators %in% 1:he$n.comparators))
    he$comp <- he$comp[comparators]
    he$n.comparators = length(comparators)
    he$n.comparisons = length(comparators) - 1
    he$interventions = he$interventions[comparators]
    he$ref = rank(c(he$ref, he$comp))[1]
    he$comp = rank(c(he$ref, he$comp))[-1]
    he$mod <- TRUE
    he$e <- he$e[, comparators]
    he$c <- he$c[, comparators]
  }
  if (relative) {
    temp <- he
    temp$e <- temp$c <- matrix(NA, he$n.sim, he$n.comparators)
    temp$e[, he$ref] <- temp$c[, he$ref] <- rep(0, he$n.sim)
    temp$e[, -he$ref] <- -he$delta.e
    temp$c[, -he$ref] <- -he$delta.c
    he <- temp
  }
  stopifnot(he$n.comparators >= 2)
  base.graphics <- ifelse(isTRUE(pmatch(graph, c("base",
                                                 "ggplot2")) == 2), FALSE, TRUE)
  c.avg <- e.avg <- x <- y <- e <- e.orig <- c.orig <- NA_real_
  ec.min <- with(he, c(min(apply(e, 2, median)), apply(c, 2,
                                                       median)[which.min(apply(e, 2, median))], which.min(apply(e,
                                                                                                                2, median))))
  e.neg <- ec.min[1] < 0
  c.neg <- any(apply(he$c, 2, median) < 0)
  if (e.neg & !c.neg & start.from.origins) {
    message("Benefits are negative, the frontier will not start from the origins")
    start.from.origins <- FALSE
  }
  if (!e.neg & c.neg & start.from.origins) {
    message("Costs are negative, the frontier will not start from the origins")
    start.from.origins <- FALSE
  }
  if (e.neg & c.neg & start.from.origins) {
    message("Costs and benefits are negative, the frontier will not start from the origins")
    start.from.origins <- FALSE
  }
  e.neg <- ifelse(start.from.origins, e.neg, TRUE)
  data.avg <- data.frame(e.avg = apply(he$e, 2, median) - ifelse(!e.neg,
                                                                 0, ec.min[1]), c.avg = apply(he$c, 2, median) - ifelse(!e.neg,
                                                                                                                        0, ec.min[2]))
  data.avg <- cbind(data.avg, data.avg, as.factor(c(1:dim(data.avg)[1])))
  names(data.avg)[3:5] <- c("e.orig", "c.orig",
                            "comp")
  orig.avg <- data.avg[, 3:5]
  comp <- ifelse(any(apply(data.avg[, 1:2], 1, function(x) isTRUE(sum(x) ==
                                                                    0 & prod(x) == 0))), which(apply(data.avg[, 1:2], 1,
                                                                                                     sum) == 0 & apply(data.avg[, 1:2], 1, prod) == 0), 0)
  ceef.points <- data.frame(x = 0, y = 0, comp = comp)
  repeat {
    if (prod(dim(data.avg)) == 0)
      break
    theta <- with(data.avg, atan(c.avg/e.avg))
    theta.min <- min(theta, na.rm = TRUE)
    if (theta.min > threshold)
      break
    index <- which(theta == theta.min)
    if (length(index) > 1)
      index = index[which.min(data.avg$e.avg[index])]
    ceef.points <- with(data.avg, rbind(ceef.points, c(e.orig[index],
                                                       c.orig[index], comp[index])))
    data.avg[, 1:2] <- data.avg[, 3:4] - matrix(rep(as.numeric(data.avg[index,
                                                                        3:4]), dim(data.avg)[1]), ncol = 2, byrow = TRUE)
    data.avg <- subset(subset(data.avg, c.avg * e.avg > 0),
                       c.avg + e.avg > 0)
  }
  ceef.points$comp <- factor(ceef.points$comp)
  ceef.points$slope <- NA
  for (i in 2:dim(ceef.points)[1]) ceef.points$slope[i] <- with(ceef.points,
                                                                (y[i] - y[i - 1])/(x[i] - x[i - 1]))
  while (dim(ceef.points)[1] > 1 & ceef.points$slope[2] < 0) {
    ceef.points <- ceef.points[-1, ]
    ceef.points$slope[1] <- NA
  }
  scatter.data <- data.frame(e = c(he$e), c = c(he$c), comp = as.factor(sort(rep(1:he$n.comparators,
                                                                                 he$n.sim))))
  ceef.points[, 1] <- ceef.points[, 1] + ifelse(!e.neg, 0,
                                                ec.min[1])
  ceef.points[, 2] <- ceef.points[, 2] + ifelse(!e.neg, 0,
                                                ec.min[2])
  orig.avg[, 1] <- orig.avg[, 1] + ifelse(!e.neg, 0, ec.min[1])
  orig.avg[, 2] <- orig.avg[, 2] + ifelse(!e.neg, 0, ec.min[2])
  ceef.summary <- function(he, ceef.points, orig.avg, include.ICER = FALSE,
                           ...) {
    no.ceef <- which(!1:he$n.comparators %in% ceef.points$comp)
    if (ceef.points$comp[1] == 0)
      ceef.points <- ceef.points[-1, ]
    rownames(ceef.points) <- he$interventions[as.numeric(levels(ceef.points$comp)[ceef.points$comp])]
    if (!include.ICER) {
      ceef.points[, 5] <- atan(ceef.points[, 4]^(1 * ifelse(!flip,
                                                            1, -1)))
      ceef.points <- ceef.points[, -3]
      colnames(ceef.points) <- c("Effectiveness",
                                 "Costs", "Increase slope", "Increase angle")
    }
    else {
      ICERs <- numeric(dim(ceef.points)[1])
      index <- as.numeric(levels(ceef.points$comp)[ceef.points$comp])
      for (i in 1:length(ICERs)) {
        if (ceef.points$comp[i] == he$ref)
          ICERs[i] <- NA_real_
        else ICERs[i] <- he$ICER[index[i] + ifelse(index[i] <
                                                     he$ref, 0, -1)]
      }
      ceef.points[, 3] <- ICERs
      ceef.points[, 5] <- atan(ceef.points[, 4]^(1 * ifelse(!flip,
                                                            1, -1)))
      colnames(ceef.points) <- c("Effectiveness",
                                 "Costs", paste0("ICER ", he$interventions[he$ref],
                                                 " vs."), "Increase slope", "Increase angle")
    }
    if (flip)
      colnames(ceef.points)[1:2] <- colnames(ceef.points[2:1])
    if (length(no.ceef) > 0) {
      noceef.points <- data.frame(matrix(NA_real_, ncol = 4,
                                         nrow = length(no.ceef)))
      noceef.points[, 1:2] <- orig.avg[no.ceef, -3]
      if (!include.ICER) {
        noceef.points <- noceef.points[, -3]
        colnames(noceef.points) <- c("Effectiveness",
                                     "Costs", "Dominance type")
      }
      else {
        ICERs <- numeric(dim(noceef.points)[1])
        for (i in 1:length(ICERs)) {
          if (no.ceef[i] == he$ref)
            ICERs[i] <- NA_real_
          else ICERs[i] <- he$ICER[no.ceef[i] + ifelse(no.ceef[i] <
                                                         he$ref, 0, -1)]
        }
        noceef.points[, 3] <- ICERs
        colnames(noceef.points) <- c("Effectiveness",
                                     "Costs", paste0("ICER ", he$interventions[he$ref],
                                                     " vs."), "Dominance type")
      }
      how.dominated <- rep("Extended dominance",
                           length(no.ceef))
      for (i in 1:length(no.ceef)) for (j in 1:dim(ceef.points)[1]) {
        if ((noceef.points[i, 1] - ceef.points[j, 1]) *
            (noceef.points[i, 2] - ceef.points[j, 2]) <
            0) {
          how.dominated[i] <- "Absolute dominance"
          break
        }
      }

      noceef.points[, ifelse(!include.ICER, 3, 4)] <- how.dominated
      rownames(noceef.points) <- he$interventions[no.ceef]
      if (flip)
        colnames(noceef.points)[1:2] <- colnames(noceef.points)[2:1]
    }

    if (length(no.ceef) > 0) {

      if(return_dominated_only == FALSE) {
        ## Print message starts here ##
        cat("\nCost-effectiveness efficiency frontier summary \n\n")
        cat("Interventions on the efficiency frontier:\n")
        print(ceef.points, quote = F, digits = 5, justify = "center")
        cat("\n")
        cat("Interventions not on the efficiency frontier:\n")
        print(noceef.points, quote = F, digits = 5, justify = "center")
        ## Print message stops here!
      } else if(return_dominated_only == TRUE) {
        ## Addition: return print message of only the dominated strategies
        print(c(rownames(noceef.points)))
      }

    }
  }
  if (print.summary)
    ceef.summary(he, ceef.points, orig.avg, ...)
  colour <- colours()[floor(seq(262, 340, length.out = he$n.comparators))]
  if (base.graphics) {
    if (print.plot) {
      if (is.numeric(pos) & length(pos) == 2) {
        temp <- ""
        if (pos[2] == 0)
          temp <- paste0(temp, "bottom")
        else temp <- paste0(temp, "top")
        if (pos[1] == 0)
          temp <- paste0(temp, "left")
        else temp <- paste0(temp, "right")
        pos <- temp
        if (length(grep("^(bottom|top)(left|right)$",
                        temp)) == 0)
          pos <- FALSE
      }
      if (is.logical(pos)) {
        if (!pos)
          pos = "topright"
        else pos = "topleft"
      }
      if (flip) {
        temp <- scatter.data$e
        scatter.data$e <- scatter.data$c
        scatter.data$c <- temp
        temp <- ceef.points$x
        ceef.points$x <- ceef.points$y
        ceef.points$y <- temp
        temp <- orig.avg$e.orig
        orig.avg$e.orig <- orig.avg$c.orig
        orig.avg$c.orig <- temp
        rm(temp)
      }
      xlab = ifelse((!flip & !relative), "Effectiveness",
                    ifelse((!flip & relative), "Differential effectiveness",
                           ifelse((flip & !relative), "Cost", "Differential cost")))
      ylab = ifelse((!flip & !relative), "Cost",
                    ifelse((!flip & relative), "Differential cost",
                           ifelse((flip & !relative), "Effectiveness",
                                  "Differential effectiveness")))
      plot(NULL, xlim = c(min(range(scatter.data$e)[1],
                              0), max(range(scatter.data$e)[2], 0)), ylim = c(min(range(scatter.data$c)[1],
                                                                                  0), max(range(scatter.data$c)[2], 0)), main = "Cost-effectiveness efficiency frontier",
           xlab = xlab, ylab = ylab)
      if (dominance) {
        for (i in 1:dim(ceef.points)[1]) {
          rect(col = "grey95", border = NA, xleft = ifelse(!flip,
                                                           -1, 1) * 2 * max(abs(range(scatter.data$e))),
               xright = ceef.points$x[i], ybottom = ceef.points$y[i],
               ytop = ifelse(!flip, 1, -1) * 2 * max(abs(range(scatter.data$c))))
        }
        if (dim(ceef.points)[1] > 1)
          for (i in 2:dim(ceef.points)[1]) {
            rect(col = "grey85", border = NA, xleft = ifelse(!flip,
                                                             -1, 1) * 2 * max(abs(range(scatter.data$e))),
                 xright = ceef.points$x[ifelse(!flip, i -
                                                 1, i)], ybottom = ceef.points$y[ifelse(!flip,
                                                                                        i, i - 1)], ytop = ifelse(!flip, 1, -1) *
                   2 * max(abs(range(scatter.data$c))))
          }
      }
      abline(h = 0, col = "grey")
      abline(v = 0, col = "grey")
      for (i in 1:he$n.comparators) with(scatter.data,
                                         points(subset(scatter.data, comp == i)[, -3],
                                                type = "p", pch = 20, cex = 0.35, col = colour[i]))
      points(ceef.points[, 1:2], type = "l", lwd = 2)
      points(orig.avg[, -3], pch = 21, cex = 2, bg = "white",
             col = "black")
      for (i in 1:he$n.comparators) {
        text(orig.avg[i, -3], labels = orig.avg[i, 3],
             col = ifelse(i %in% ceef.points$comp, "black",
                          "grey60"), cex = 0.75)
      }
      text <- paste(1:he$n.comparators, ":", he$interventions)
      legend(pos, text, col = colour, cex = 0.7, bty = "n",
             lty = 1)
      box()
    }
  }
  else {
    if (print.plot) {
      if (!isTRUE(requireNamespace("ggplot2", quietly = TRUE) &
                  requireNamespace("grid", quietly = TRUE))) {
        message("Falling back to base graphics\n")
        ceef.plot(he, flip = flip, comparators = comparators,
                  pos = pos, start.from.origins = start.from.origins,
                  graph = "base")
        return(invisible(NULL))
      }
      opt.theme <- ggplot2::theme()
      exArgs <- list(...)
      if (length(exArgs) >= 1) {
        for (obj in exArgs) if (ggplot2::is.theme(obj))
          opt.theme <- opt.theme + obj
      }
      ceplane <- ggplot2::ggplot(ceef.points, ggplot2::aes(x = x,
                                                           y = y))
      if (dominance) {
        ceplane <- ceplane + ggplot2::geom_rect(data = ceef.points,
                                                ggplot2::aes(xmax = x, ymin = y), ymax = 2 *
                                                  max(abs(range(scatter.data$c))), xmin = -2 *
                                                  max(abs(range(scatter.data$e))), alpha = 0.35,
                                                fill = "grey75")
      }
      ceplane <- ceplane + ggplot2::geom_hline(yintercept = 0,
                                               colour = "grey") + ggplot2::geom_vline(xintercept = 0,
                                                                                      colour = "grey") + ggplot2::geom_point(data = scatter.data,
                                                                                                                             ggplot2::aes(x = e, y = c, colour = comp), size = 1)
      if (dim(ceef.points)[1] > 1)
        ceplane <- ceplane + ggplot2::geom_path()
      xlab = ifelse(!relative, "Effectiveness", "Effectiveness differential")
      ylab = ifelse(!relative, "Cost", "Cost differential")
      ceplane <- ceplane + ggplot2::geom_point(data = orig.avg,
                                               ggplot2::aes(x = e.orig, y = c.orig), size = 5.5,
                                               colour = "black") + ggplot2::geom_point(data = orig.avg,
                                                                                       ggplot2::aes(x = e.orig, y = c.orig), size = 4.5,
                                                                                       colour = "white") + ggplot2::scale_colour_manual("",
                                                                                                                                        labels = paste0(1:he$n.comparators, ": ",
                                                                                                                                                        he$interventions), values = colour, na.value = "black") +
        ggplot2::labs(title = "Cost-effectiveness efficiency frontier",
                      x = xlab, y = ylab) + ggplot2::theme_bw()
      for (i in 1:he$n.comparators) {
        ceplane <- ceplane + ggplot2::geom_text(data = orig.avg[i,
                                                                ], ggplot2::aes(x = e.orig, y = c.orig, label = comp),
                                                size = 3.5, colour = ifelse(i %in% ceef.points$comp,
                                                                            "black", "grey60"))
      }
      jus <- NULL
      if (isTRUE(pos)) {
        pos = "bottom"
        ceplane <- ceplane + ggplot2::theme(legend.direction = "vertical")
      }
      else {
        if (is.character(pos)) {
          choices <- c("left", "right", "bottom",
                       "top")
          pos <- choices[pmatch(pos, choices)]
          jus = "center"
          if (is.na(pos))
            pos = FALSE
        }
        if (length(pos) > 1)
          jus <- pos
        if (length(pos) == 1 & !is.character(pos)) {
          pos <- c(1, 1)
          jus <- pos
        }
      }
      ceplane <- ceplane + ggplot2::theme(legend.position = pos,
                                          legend.justification = jus, legend.title = ggplot2::element_blank(),
                                          legend.background = ggplot2::element_blank(),
                                          text = ggplot2::element_text(size = 11), legend.key.size = grid::unit(0.66,
                                                                                                                "lines"), legend.spacing = grid::unit(-1.25,
                                                                                                                                                      "line"), panel.grid = ggplot2::element_blank(),
                                          legend.key = ggplot2::element_blank(), legend.text.align = 0,
                                          plot.title = ggplot2::element_text(hjust = 0.5,
                                                                             face = "bold", lineheight = 1.05, size = 14.3)) +
        opt.theme
      if (flip)
        ceplane <- ceplane + ggplot2::coord_flip()
      print(ceplane)
    }
  }
}

# Function to assign dominated status to scenario for a single simulation (calls ceef.plot.median)
assign_dominated_strategies <- function(df, exposure, outcome) {

  # Need to artificially repeat value twice in e and c:
  # Effectiveness (outcome) matrix:
  outcome_matrix <- rbind(as.matrix(select(df, scenario, sim, outcome) %>%
                                      spread(key="scenario", value = outcome)%>%
                                      select(-sim)),
                          as.matrix(select(df, scenario, sim,outcome) %>%
                                      spread(key="scenario", value =outcome)%>%
                                      select(-sim)))
  # Cost (exposure) matrix:
  exposure_matrix <- rbind(as.matrix(select(df, scenario, sim, exposure) %>%
                                       spread(key="scenario", value = exposure)%>%
                                       select(-sim)),
                           as.matrix(select(df, scenario, sim,exposure) %>%
                                       spread(key="scenario", value =exposure)%>%
                                       select(-sim)))
  # Assign numbers to each scenario
  interventions_df <- data.frame(scenario = c(colnames(select(df, scenario, sim, outcome) %>%
                                                         spread(key="scenario", value = outcome)%>%
                                                         select(-sim))),
                                 numbers = 1:length(c(colnames(select(df, scenario, sim, outcome) %>%
                                                                 spread(key="scenario", value = outcome)%>%
                                                                 select(-sim)))))


  ceef <- capture.output(ceef.plot.median(bcea(e=outcome_matrix,
                                               c=exposure_matrix,
                                               ref=1,
                                               interventions=interventions_df$numbers,
                                               Kmax=100000000000000,
                                               plot=FALSE),
                                          graph="base", relative = FALSE, return_dominated_only = TRUE))

  # Extracted numbers of dominated interventions
  dominated_numbers <- unlist(regmatches(ceef, gregexpr("[[:digit:]]+", ceef)))[-1]
  # Assign this to scenario
  interventions_df$dominated <- "No"
  interventions_df$dominated[interventions_df$numbers %in% dominated_numbers] <- "Yes"
  # Merge with original df

  dom_df <- left_join(df, interventions_df, by = "scenario") %>%
    select(-numbers)

  return(dom_df)

}

# Function to calculate ICER for specified strategies for a single simulation
# For this, need to manually decide on dominated strategies to exclude and remove these from arguments!
calculate_icer_per_sim <- function(df, exposure, outcome) {

  ranked_strategies <- arrange(df, get(exposure))
  ranked_strategies$diff_exposure <- c(ranked_strategies[,exposure][1], diff(ranked_strategies[,exposure]))
  ranked_strategies$diff_outcome <- c(ranked_strategies[,outcome][1], diff(ranked_strategies[,outcome]))
  ranked_strategies$icer <-ranked_strategies$diff_exposure/ranked_strategies$diff_outcome
  ranked_strategies$icer[is.na(ranked_strategies$icer)] <- 0
  ranked_strategies$comparator <- lag(ranked_strategies$scenario)
  icer <- select(ranked_strategies, -diff_exposure)  # -diff_outcome

  return(icer)

}

# Function to find dominated and extended dominated strategies on cost-effectiveness frontier
find_dominated_strategies <- function(df, exposure, outcome) {
  plot <- ceef.plot.median(bcea(e=cbind(matrix(rep(0, 183)),
                                        as.matrix(select(df, scenario, sim, outcome) %>%
                                                    filter(sim != "X") %>%
                                                    spread(key="scenario", value = outcome)%>%
                                                    select(-sim))),
                                c=cbind(matrix(rep(0, 183)),
                                        as.matrix(select(df, scenario, sim, exposure) %>%
                                                    filter(sim != "X") %>%
                                                    spread(key="scenario", value = exposure)%>%
                                                    select(-sim))),
                                ref=1,
                                interventions=c("No monitoring", colnames(select(df, scenario, sim, outcome) %>%
                                                                            filter(sim != "X") %>%
                                                                            spread(key="scenario", value = outcome)%>%
                                                                            select(-sim))),
                                Kmax=1000000000,
                                plot=FALSE),
                           graph="base", relative = FALSE)
  return(plot)
}

# Function to create dataframe for incremental monitoring strategies plots
create_incremental_plot_df <- function(interactions_df, py_on_treatment_df,
                                       deaths_averted_df, ly_saved_df,
                                       hbsag_test_cost = 10.38,
                                       clinical_assessment_cost = 120,
                                       monitoring_assessment_cost = 15.77,
                                       treatment_py_cost = 66.44, # changed to reflect only 1 monitoring assessment on treatment
                                       scenario_labels_obj = NULL,
                                       ref_label) {
  # interactions_df has columns: scenario, sim, monitoring_assessments, treatment_initiations
  # If reference is status quo, also needs to have hbsag_tests and clinical_assessments column
  # py_on_treatment_df has columns:scenario, sim, py_on_treatment
  # deaths_averted_df and ly_saved_df has columns: scenario, sim, value (sim already in correct format)
  # Ensure all inputs have the same scenarios!
  # If the reference is status quo, ref_label has to be "No treatment"
  interactions <- interactions_df

  if (ref_label == "No treatment") {

    # Previously, included treatment initiation as an interactions rather than monitoring (as represented by PY)
    #   interactions$total_interactions <- interactions$hbsag_tests+
    #      interactions$clinical_assessments+
    #      interactions$monitoring_assessments+
    #      interactions$treatment_initiations
    interactions$sim <- gsub("[^0-9]", "", interactions$sim)

    # Combine interactions and assign costs:
    interactions <- left_join(interactions, py_on_treatment_df,
                              by=c("scenario", "sim")) %>%
      mutate(total_interactions = hbsag_tests+clinical_assessments+monitoring_assessments+
               py_on_treatment*1,   # monitoring on treatment once a year. Includes treatment initiation.
             screening_cost = hbsag_tests*hbsag_test_cost,
             assessment_cost = clinical_assessments*clinical_assessment_cost,
             monitoring_cost = monitoring_assessments*monitoring_assessment_cost,
             treatment_cost = py_on_treatment*treatment_py_cost,
             total_cost = screening_cost+assessment_cost+monitoring_cost+treatment_cost)

  } else {

    #  interactions$total_interactions <- interactions$monitoring_assessments+
    #    interactions$treatment_initiations
    interactions$sim <- gsub("[^0-9]", "", interactions$sim)

    # Combine interactions and assign costs based on Shevanthi's Gambia paper:
    # $15.77 per monitoring assessment and $66.44 per person-year of treatment
    interactions <- left_join(interactions, py_on_treatment_df,
                              by=c("scenario", "sim")) %>%
      mutate(total_interactions = monitoring_assessments+
               py_on_treatment*1, # represents treatment initiation + monitoring on treatment ONCE a year
             monitoring_cost = monitoring_assessments*monitoring_assessment_cost,
             treatment_cost = py_on_treatment*treatment_py_cost,
             total_cost = monitoring_cost+treatment_cost)

  }

  # Combine with deaths averted
  df <- deaths_averted_df

  df <- df %>%
    left_join(interactions, by = c("scenario", "sim")) %>%
    drop_na()
  colnames(df)[colnames(df)=="value"] <- "deaths_averted"

  # Combine with LY saved
  df_ly <- ly_saved_df

  df <- df_ly %>%
    left_join(df, by = c("scenario", "sim")) %>%
    drop_na()
  colnames(df)[colnames(df)=="value"] <- "ly_saved"

  # Scenario labels
  if(length(scenario_labels_obj) > 0L) {
    levels(df$scenario) <- scenario_labels_obj
  }

  # Separate columns for frequency and age group
  df$scenario2 <- df$scenario
  df <- extract(df, col="scenario2", into = c("frequency","age_group"), regex="^(\\S+)\\s+(.*)")

  # Add reference point at 0
  ref_df <- data.frame(matrix(c(ref_label, "X", rep(0,(ncol(df)-4))), nrow = 1))
  ref_df <- ref_df[rep(seq_len(nrow(ref_df)), each = nrow(crossing(df$frequency, df$age_group))),]
  ref_df <- cbind(ref_df, data.frame(crossing(df$frequency, df$age_group)))
  colnames(ref_df) <- colnames(df)

  df$scenario <- as.character(df$scenario)
  df <- rbind(df, ref_df)
  df[,-c(1:2,ncol(df)-1,ncol(df))] <- apply(df[,-c(1:2,ncol(df)-1,ncol(df))], 2, as.numeric)
  df$scenario <- factor(df$scenario)

  return(df)


}

# Function to discount outcomes (deaths averted or LY saved) or interactions
discount_outcome_2020_to_2100 <- function(scenario_object,object_to_subtract,
                                          outcome, interaction_outcome, yearly_discount_rate,
                                          interaction_colname) {
  # This function takes output calculated every 5 years from 2025 to 2100 and
  # sums to outcome between 2020 and 2100
  # comparator_object should be the one with the higher value  # outcome can be "cum_hbv_deaths" or "ly"

  # Check output is every 5 years from 2025 to 2100:
  # Doing this on cum_hbv_deaths because py_on_treatment does not have year attached
  if((all.equal(sapply(object_to_subtract[["cum_hbv_deaths"]], "[[", "by_year"),seq(2025,2100,by=5))==TRUE ||
      is.null(object_to_subtract)) &
     all.equal(sapply(scenario_object[["cum_hbv_deaths"]], "[[", "by_year"),seq(2025,2100,by=5))==TRUE) {

    if(outcome %in% c("cum_hbv_deaths", "ly", "yll", "dalys")) {

      # Transform cumulative outcome into incident deaths averted/life-years saved
      # (from population output, not cohort - though cumulative deaths/LY is not the same,
      # the deaths averted/LY saved between 2 strategies is)
      if(is.null(object_to_subtract)==TRUE) {
        cum_outcome <- do.call("rbind", scenario_object[[outcome]])[,-c(1:3)]
      } else if(is.null(object_to_subtract)==FALSE) {
        cum_outcome <- do.call("rbind", scenario_object[[outcome]])[,-c(1:3)]-
          do.call("rbind", object_to_subtract[[outcome]])[,-c(1:3)]
      }

      inc_outcome <- data.frame(rbind(cum_outcome[1,], apply(cum_outcome,2,diff)))

    } else if (outcome == "py_on_treatment") {

      # Transform cumulative outcome into incident deaths averted/life-years saved
      # (from population output, not cohort - though cumulative deaths/LY is not the same,
      # the deaths averted/LY saved between 2 strategies is)
      if(is.null(object_to_subtract)==TRUE) {
        cum_outcome <- do.call("rbind", scenario_object[[outcome]])
      } else if(is.null(object_to_subtract)==FALSE) {
        cum_outcome <- do.call("rbind", scenario_object[[outcome]])-
          do.call("rbind", object_to_subtract[[outcome]])
      }

      inc_outcome <- data.frame(rbind(cum_outcome[1,], apply(cum_outcome,2,diff)))

    } else if (outcome == "interactions") {

      # Transform cumulative outcome into incidence
      # (from population output, not cohort - though cumulative deaths/LY is not the same,
      # the deaths averted/LY saved between 2 strategies is)
      if(is.null(object_to_subtract)==TRUE) {
        cum_outcome <- do.call("rbind", lapply(scenario_object[[outcome]], "[[", interaction_outcome))[,-c(1:3)]
      } else if(is.null(object_to_subtract)==FALSE) {
        cum_outcome <- do.call("rbind", lapply(scenario_object[[outcome]], "[[", interaction_outcome))[,-c(1:3)]-
          do.call("rbind", lapply(object_to_subtract[[outcome]], "[[", interaction_outcome))[,-c(1:3)]
      }
      inc_outcome <- data.frame(rbind(cum_outcome[1,], apply(cum_outcome,2,diff)))
    }

    # Calculate the yearly discount rate, then take the mean for every 5 years.
    # This assumes that the indicence is constant over those 5 years, though that is not necessarily true.
    total_years <- c(2020:2099)
    n_total_years <- length(total_years) - 1
    discount_df <- data.frame(total_years=total_years,
                              discount_by_year = (1-yearly_discount_rate) ^ (0:n_total_years))
    discount_df$group_5 <- rep(1:(nrow(discount_df)/5), each = 5)
    discount_vec <- group_by(discount_df, group_5) %>%
      summarise(discounts = mean(discount_by_year))
    # Apply this to the incident outcome every 5 years, then sum all incident events between 2020 to 2100
    # to get the cumulative discounted deaths averted/LY saved
    discounted_cum_outcome <- apply(discount_vec$discounts*inc_outcome,2,sum)
    discounted_cum_outcome <- as.data.frame(discounted_cum_outcome)
    colnames(discounted_cum_outcome) <- ifelse(test <- outcome=="interactions", yes=interaction_colname,
                                               no=outcome)
    discounted_cum_outcome$sim <- rownames(discounted_cum_outcome)
    discounted_cum_outcome <- discounted_cum_outcome[,c(2,1)]
    return(discounted_cum_outcome)


  }  else {
    print("Years are not from 2020 to 2100.")
  }


}


# Load files ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# Status quo
out1_it <- readRDS(paste0(out_path, "a1_it_out1_status_quo_cohort_200121.rds"))
out1_it <- out1_it[[1]]
out2 <- readRDS(paste0(out_path, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]

# FILES FOR ANALYSIS ACROSS ALL AGES

# No IT treatment

# New files with correct counting of monitoring events (for analysis of
# cost-effectiveness of IT treatment)
out3n <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_201220.rds"))
out3n <- out3n[[1]]
out5n <- readRDS(paste0(out_path, "a1_out5_screen_2020_monit_5_201220.rds"))
out5n <- out5n[[1]]
out6an <- readRDS(paste0(out_path, "a1_out6a_screen_2020_monit_2_181220.rds"))
out6an <- out6an[[1]]

# No monitoring
out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_011220.rds"))
out3 <- out3[[1]]

# Monitoring all age groups from every 30 to every 1 year (need to update this with IT)
out4d <- readRDS(paste0(out_path, "a1_out4d_screen_2020_monit_30_031220.rds"))
out4d <- out4d[[1]]
out4c <- readRDS(paste0(out_path, "a1_out4c_screen_2020_monit_25_031220.rds"))
out4c <- out4c[[1]]
out4b <- readRDS(paste0(out_path, "a1_out4b_screen_2020_monit_20_021220.rds"))
out4b <- out4b[[1]]
out4a <- readRDS(paste0(out_path, "a1_out4a_screen_2020_monit_15_021220.rds"))
out4a <- out4a[[1]]
out4 <- readRDS(paste0(out_path, "a1_out4_screen_2020_monit_10_301120.rds"))
out4 <- out4[[1]]
out5d <- readRDS(paste0(out_path, "a1_out5d_screen_2020_monit_9_301120.rds"))
out5d <- out5d[[1]]
out5c <- readRDS(paste0(out_path, "a1_out5c_screen_2020_monit_8_301120.rds"))
out5c <- out5c[[1]]
out5b <- readRDS(paste0(out_path, "a1_out5b_screen_2020_monit_7_301120.rds"))
out5b <- out5b[[1]]
out5a <- readRDS(paste0(out_path, "a1_out5a_screen_2020_monit_6_011220.rds"))
out5a <- out5a[[1]]
out5 <- readRDS(paste0(out_path, "a1_out5_screen_2020_monit_5_301120.rds"))
out5 <- out5[[1]]
out6c <- readRDS(paste0(out_path, "a1_out6c_screen_2020_monit_4_301120.rds"))
out6c <- out6c[[1]]
out6b <- readRDS(paste0(out_path, "a1_out6b_screen_2020_monit_3_301120.rds"))
out6b <- out6b[[1]]
out6a <- readRDS(paste0(out_path, "a1_out6a_screen_2020_monit_2_301120.rds"))
out6a <- out6a[[1]]
out6 <- readRDS(paste0(out_path, "a1_out6_screen_2020_monit_1_301120.rds"))
out6 <- out6[[1]]

# Once in a lifetime monitoring at different ages (replaced with new sims with IT)
# These are wrong (modelling of lifetime monitoring not comparable with other frequencies)
#out_lt_monit_30 <- readRDS(paste0(out_path, "a1_out_screen_2020_monit_lifetime_30_101220.rds"))
#out_lt_monit_30 <- out_lt_monit_30[[1]]
#out_lt_monit_45 <- readRDS(paste0(out_path, "a1_out_screen_2020_monit_lifetime_45_101220.rds"))
#out_lt_monit_45 <- out_lt_monit_45[[1]]

# FILES FOR AGE-SPECIFIC ANALYSIS
# IT is treated

# Effect of screening by age (not monitoring)
a2_out3_it <- readRDS(paste0(out_path, "a2_it_out3_screen_2020_monit_0_180121.rds"))
a2_out3_it <- a2_out3_it[[1]]
a4_out3_it <- readRDS(paste0(out_path, "a4_it_out3_screen_2020_monit_0_190121.rds"))
a4_out3_it <- a4_out3_it[[1]]
a5_out3_it <- readRDS(paste0(out_path, "a5_it_out3_screen_2020_monit_0_190121.rds"))
a5_out3_it <- a5_out3_it[[1]]

# No monitoring/All age monitoring
out3_it <- readRDS(paste0(out_path, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]   # No monitoring
out4_it <- readRDS(paste0(out_path, "a1_it_out4_screen_2020_monit_10_140121.rds"))
out4_it <- out4_it[[1]]   # 10 years
out5_it <- readRDS(paste0(out_path, "a1_it_out5_screen_2020_monit_5_161220.rds"))
out5_it <- out5_it[[1]]   # 5 years
out6c_it <- readRDS(paste0(out_path, "a1_it_out6c_screen_2020_monit_4_130121.rds"))
out6c_it <- out6c_it[[1]]   # 4 years
out6b_it <- readRDS(paste0(out_path, "a1_it_out6b_screen_2020_monit_3_130121.rds"))
out6b_it <- out6b_it[[1]]   # 3 years
out6a_it <- readRDS(paste0(out_path, "a1_it_out6a_screen_2020_monit_2_161220.rds"))
out6a_it <- out6a_it[[1]]  # 2 years
out6_it <- readRDS(paste0(out_path, "a1_it_out6_screen_2020_monit_1_130121.rds"))
out6_it <- out6_it[[1]]  # 1 year

# Monitoring different age groups while they are of given age (IT always treated)
monit_out7 <- readRDS(paste0(out_path, "a1_it_monit_out7_161220.rds"))
monit_out7 <- monit_out7[[1]]

monit_out17 <- readRDS(paste0(out_path, "a1_it_monit_out17_161220.rds"))
monit_out17 <- monit_out17[[1]]

monit_out2a <- readRDS(paste0(out_path, "a1_it_monit_out2a_171220.rds"))
monit_out2a <- monit_out2a[[1]]

monit_out9 <- readRDS(paste0(out_path, "a1_it_monit_out9_171220.rds"))
monit_out9 <- monit_out9[[1]]

monit_out18 <- readRDS(paste0(out_path, "a1_it_monit_out18_171220.rds"))
monit_out18 <- monit_out18[[1]]

monit_out4a <- readRDS(paste0(out_path, "a1_it_monit_out4a_171220.rds"))
monit_out4a <- monit_out4a[[1]]

monit_out10 <- readRDS(paste0(out_path, "a1_it_monit_out10_200121.rds"))
monit_out10 <- monit_out10[[1]]

monit_out5a <- readRDS(paste0(out_path, "a1_it_monit_out5a_161220.rds"))
monit_out5a <- monit_out5a[[1]]

monit_out6 <- readRDS(paste0(out_path, "a1_it_monit_out6_200121.rds"))
monit_out6 <- monit_out6[[1]]

monit_out1a <- readRDS(paste0(out_path, "a1_it_monit_out1a_191220.rds"))
monit_out1a <- monit_out1a[[1]]

monit_out8 <- readRDS(paste0(out_path, "a1_it_monit_out8_090121.rds"))
monit_out8 <- monit_out8[[1]]

# 10-yearly frequencies
monit_out7_10 <- readRDS(paste0(out_path, "a1_it_monit_out7_10_150121.rds"))
monit_out7_10 <- monit_out7_10[[1]]

# 4-yearly frequencies
monit_out4c <- readRDS(paste0(out_path, "a1_it_monit_out4c_130121.rds"))
monit_out4c <- monit_out4c[[1]]
monit_out5c <- readRDS(paste0(out_path, "a1_it_monit_out5c_130121.rds"))
monit_out5c <- monit_out5c[[1]]
monit_out2c <- readRDS(paste0(out_path, "a1_it_monit_out2c_140121.rds"))
monit_out2c <- monit_out2c[[1]]
monit_out1c <- readRDS(paste0(out_path, "a1_it_monit_out1c_140121.rds"))
monit_out1c <- monit_out1c[[1]]

# 3-yearly frequencies
monit_out4b <- readRDS(paste0(out_path, "a1_it_monit_out4b_130121.rds"))
monit_out4b <- monit_out4b[[1]]
monit_out2b <- readRDS(paste0(out_path, "a1_it_monit_out2b_130121.rds"))
monit_out2b <- monit_out2b[[1]]
monit_out5b <- readRDS(paste0(out_path, "a1_it_monit_out5b_150121.rds"))
monit_out5b <- monit_out5b[[1]]
monit_out1b <- readRDS(paste0(out_path, "a1_it_monit_out1b_150121.rds"))
monit_out1b <- monit_out1b[[1]]

# Annual frequencies
monit_out4 <- readRDS(paste0(out_path, "a1_it_monit_out4_140121.rds"))
monit_out4 <- monit_out4[[1]]
monit_out5 <- readRDS(paste0(out_path, "a1_it_monit_out5_140121.rds"))
monit_out5 <- monit_out5[[1]]
monit_out2 <- readRDS(paste0(out_path, "a1_it_monit_out2_130121.rds"))
monit_out2 <- monit_out2[[1]]
monit_out1 <- readRDS(paste0(out_path, "a1_it_monit_out1_130121.rds"))
monit_out1 <- monit_out1[[1]]

# Lifetime monitoring events at specific age with IT treatment
# These are wrong/not comparable!!
out_it_lt_monit_30 <- readRDS(paste0(out_path, "a1_it_out_screen_2020_monit_lifetime_30_191220.rds"))
out_it_lt_monit_30 <- out_it_lt_monit_30[[1]]
out_it_lt_monit_45 <- readRDS(paste0(out_path, "a1_it_out_screen_2020_monit_lifetime_45_191220.rds"))
out_it_lt_monit_45 <- out_it_lt_monit_45[[1]]
out_it_lt_monit_30_45 <- readRDS(paste0(out_path, "a1_it_out_screen_2020_monit_lifetime_30_45_181220.rds"))
out_it_lt_monit_30_45 <- out_it_lt_monit_30_45[[1]]

scenario_labels <- list("No treatment" = "status_quo_cohort",
                        "No monitoring" = "screen_2020_monit_0",
                        "5-yearly all ages"="screen_2020_monit_5",
                        "Yearly all ages"="screen_2020_monit_1",
                        "Yearly 15-30"="screen_2020_monit_sim1",
                        "Yearly 15-45"="screen_2020_monit_sim2",
                        "Yearly 30-45"="screen_2020_monit_sim3",
                        "Yearly 30+"="screen_2020_monit_sim4",
                        "Yearly 45+"="screen_2020_monit_sim5",
                        "5-yearly 15-30"="screen_2020_monit_sim6",
                        "5-yearly 15-45"= "screen_2020_monit_sim7",
                        "5-yearly 30-45"="screen_2020_monit_sim8",
                        "5-yearly 30+"="screen_2020_monit_sim9",
                        "5-yearly 45+"="screen_2020_monit_sim10")

## 1) Cost-effectiveness of monitoring frequencies across all ages ----
# Create data frame with all interactions and outcomes of interest (no discounting) ----
freq_interactions <- rbind(
  cbind(scenario = "No monitoring",
        left_join(left_join(left_join(gather(out3$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "At age 30",
        left_join(left_join(left_join(gather(out_lt_monit_30$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out_lt_monit_30$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out_lt_monit_30$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "At age 45",
        left_join(left_join(left_join(gather(out_lt_monit_45$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out_lt_monit_45$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out_lt_monit_45$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 30 years",
        left_join(left_join(left_join(gather(out4d$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4d$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4d$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 25 years",
        left_join(left_join(left_join(gather(out4c$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4c$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4c$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 20 years",
        left_join(left_join(left_join(gather(out4b$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4b$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4b$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 15 years",
        left_join(left_join(left_join(gather(out4a$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4a$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4a$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 10 years",
        left_join(left_join(left_join(gather(out4$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 9 years",
        left_join(left_join(left_join(gather(out5d$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5d$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5d$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 8 years",
        left_join(left_join(left_join(gather(out5c$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5c$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5c$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 7 years",
        left_join(left_join(left_join(gather(out5b$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5b$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5b$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 6 years",
        left_join(left_join(left_join(gather(out5a$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5a$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5a$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 5 years",
        left_join(left_join(left_join(gather(out5$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 4 years",
        left_join(left_join(left_join(gather(out6c$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out6c$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out6c$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 3 years",
        left_join(left_join(left_join(gather(out6b$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out6b$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out6b$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 2 years",
        left_join(left_join(left_join(gather(out6a$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out6a$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out6a$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Every 1 year",
        left_join(left_join(left_join(gather(out6$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out6$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")))


# Extract person-years on treatment
freq_interactions_py_on_treatment <- rbind(
  data.frame(scenario = "No monitoring",
             sim = names(out3$py_on_treatment[[16]]),
             py_on_treatment = out3$py_on_treatment[[16]]),
  data.frame(scenario = "At age 30",
             sim = names(out_lt_monit_30$py_on_treatment[[16]]),
             py_on_treatment = out_lt_monit_30$py_on_treatment[[16]]),
  data.frame(scenario = "At age 45",
             sim = names(out_lt_monit_45$py_on_treatment[[16]]),
             py_on_treatment = out_lt_monit_45$py_on_treatment[[16]]),
  data.frame(scenario = "Every 30 years",
             sim = names(out4d$py_on_treatment[[16]]),
             py_on_treatment = out4d$py_on_treatment[[16]]),
  data.frame(scenario = "Every 25 years",
             sim = names(out4c$py_on_treatment[[16]]),
             py_on_treatment = out4c$py_on_treatment[[16]]),
  data.frame(scenario = "Every 20 years",
             sim = names(out4b$py_on_treatment[[16]]),
             py_on_treatment = out4b$py_on_treatment[[16]]),
  data.frame(scenario = "Every 15 years",
             sim = names(out4a$py_on_treatment[[16]]),
             py_on_treatment = out4a$py_on_treatment[[16]]),
  data.frame(scenario = "Every 10 years",
             sim = names(out4$py_on_treatment[[16]]),
             py_on_treatment = out4$py_on_treatment[[16]]),
  data.frame(scenario = "Every 9 years",
             sim = names(out5d$py_on_treatment[[16]]),
             py_on_treatment = out5d$py_on_treatment[[16]]),
  data.frame(scenario = "Every 8 years",
             sim = names(out5c$py_on_treatment[[16]]),
             py_on_treatment = out5c$py_on_treatment[[16]]),
  data.frame(scenario = "Every 7 years",
             sim = names(out5b$py_on_treatment[[16]]),
             py_on_treatment = out5b$py_on_treatment[[16]]),
  data.frame(scenario = "Every 6 years",
             sim = names(out5a$py_on_treatment[[16]]),
             py_on_treatment = out5a$py_on_treatment[[16]]),
  data.frame(scenario = "Every 5 years",
             sim = names(out5$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]]),
  data.frame(scenario = "Every 4 years",
             sim = names(out6c$py_on_treatment[[16]]),
             py_on_treatment = out6c$py_on_treatment[[16]]),
  data.frame(scenario = "Every 3 years",
             sim = names(out6b$py_on_treatment[[16]]),
             py_on_treatment = out6b$py_on_treatment[[16]]),
  data.frame(scenario = "Every 2 years",
             sim = names(out6a$py_on_treatment[[16]]),
             py_on_treatment = out6a$py_on_treatment[[16]]),
  data.frame(scenario = "Every 1 year",
             sim = names(out6$py_on_treatment[[16]]),
             py_on_treatment = out6$py_on_treatment[[16]])
)

# Outcome 1: DALYs
freq_dalys_averted <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(out3, out_lt_monit_30,
                                                  out_lt_monit_45,out4,out5,out5a,
                                                  out4d, out4c, out4b, out4a,
                                                  out5b,out5c,out5d,
                                                  out6,out6a,out6b,out6c),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
freq_dalys_averted <- subset(freq_dalys_averted, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
freq_dalys_averted$sim <- gsub("[^0-9]", "", freq_dalys_averted$sim)

monitoring_freq_label <- list("No monitoring"="screen_2020_monit_0",
                              "At age 30"="screen_2020_monit_lifetime_30",
                              "At age 45"="screen_2020_monit_lifetime_45",
                              "Every 30 years"="screen_2020_monit_30",
                              "Every 25 years"="screen_2020_monit_25",
                              "Every 20 years"="screen_2020_monit_20",
                              "Every 15 years"="screen_2020_monit_15",
                              "Every 10 years"="screen_2020_monit_10",
                              "Every 9 years"="screen_2020_monit_9",
                              "Every 8 years"="screen_2020_monit_8",
                              "Every 7 years"="screen_2020_monit_7",
                              "Every 6 years"="screen_2020_monit_6",
                              "Every 5 years"="screen_2020_monit_5",
                              "Every 4 years"="screen_2020_monit_4",
                              "Every 3 years"="screen_2020_monit_3",
                              "Every 2 years"="screen_2020_monit_2",
                              "Every 1 year"="screen_2020_monit_1")
freq_dalys_averted$scenario <- factor(freq_dalys_averted$scenario)
levels(freq_dalys_averted$scenario) <- monitoring_freq_label

# Cohort DALYS (to add lifetime monitoring)
freq_dalys_averted_cohort <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
                                 scenario_objects = list(out3,out4,out5,out5a,
                                                         out4d, out4c, out4b, out4a,
                                                         out5b,out5c,out5d,
                                                         out6,out6a,out6b,out6c),
                                 outcome_to_avert = "cohort_dalys",
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no treatment")
freq_dalys_averted_cohort <- subset(freq_dalys_averted_cohort, type == "number_averted") %>%
  select(scenario, sim, value)
freq_dalys_averted_cohort$sim <- gsub("[^0-9]", "", freq_dalys_averted_cohort$sim)
freq_dalys_averted_cohort$scenario <- factor(freq_dalys_averted_cohort$scenario)
levels(freq_dalys_averted_cohort$scenario) <- monitoring_freq_label


# Outcome 2: HBV-related deaths
freq_hbv_deaths_averted <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(out3,out_lt_monit_30, out_lt_monit_45,
                                                  out4,out5,out5a,
                                                  out4d, out4c, out4b, out4a,
                                                  out5b,out5c,out5d,
                                                  out6,out6a,out6b,out6c),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
freq_hbv_deaths_averted <- subset(freq_hbv_deaths_averted, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
freq_hbv_deaths_averted$sim <- gsub("[^0-9]", "", freq_hbv_deaths_averted$sim)
freq_hbv_deaths_averted$scenario <- factor(freq_hbv_deaths_averted$scenario)
levels(freq_hbv_deaths_averted$scenario) <- monitoring_freq_label

# Extra outcome: Proportion of HBV-related deaths in the cohort
# (to add lifetime monitoring)
freq_hbv_deaths_averted_cohort <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
                                 scenario_objects = list(out3,out4,out5,out5a,
                                                         out4d, out4c, out4b, out4a,
                                                         out5b,out5c,out5d,
                                                         out6,out6a,out6b,out6c),
                                 outcome_to_avert = "cohort_cum_hbv_deaths",
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no treatment")
freq_hbv_deaths_averted_cohort <- subset(freq_hbv_deaths_averted_cohort,
                                         type == "proportion_averted") %>%
  select(scenario, sim, value)
freq_hbv_deaths_averted_cohort$sim <- gsub("[^0-9]", "", freq_hbv_deaths_averted_cohort$sim)
freq_hbv_deaths_averted_cohort$scenario <- factor(freq_hbv_deaths_averted_cohort$scenario)
levels(freq_hbv_deaths_averted_cohort$scenario) <- monitoring_freq_label

# Create data frame with discounted interactions and outcomes ----

annual_discounting_rate <- 0.03

# This function is made for this specific analysis:
assemble_discounted_interactions_for_monitoring_frequencies <- function(scenario_object,
                                                                        discount_rate = annual_discounting_rate,
                                                                        no_monitoring_object = out3) {
  # Compares monitoring to out3 and everything else to status quo
  # at given discount rate
  out <- left_join(
    left_join(
      left_join(discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_screened",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "hbsag_tests"),
                discount_outcome_2020_to_2100(scenario_object=no_monitoring_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_assessed",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "clinical_assessments")),
      discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                    object_to_subtract=no_monitoring_object,
                                    outcome="interactions",
                                    interaction_outcome="total_assessed",
                                    yearly_discount_rate=discount_rate,
                                    interaction_colname = "monitoring_assessments")),
    discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                  object_to_subtract=NULL,
                                  outcome="interactions",
                                  interaction_outcome="total_treated",
                                  yearly_discount_rate=discount_rate,
                                  interaction_colname = "treatment_initiations"))

  return(out)
}

freq_interactions_disc <- rbind(
  cbind(scenario = "No monitoring",
        assemble_discounted_interactions_for_monitoring_frequencies(out3)),
  cbind(scenario = "At age 30",
        assemble_discounted_interactions_for_monitoring_frequencies(out_lt_monit_30)),
  cbind(scenario = "At age 45",
        assemble_discounted_interactions_for_monitoring_frequencies(out_lt_monit_45)),
  cbind(scenario = "Every 30 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out4d)),
  cbind(scenario = "Every 25 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out4c)),
  cbind(scenario = "Every 20 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out4b)),
  cbind(scenario = "Every 15 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out4a)),
  cbind(scenario = "Every 10 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out4)),
  cbind(scenario = "Every 9 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out5d)),
  cbind(scenario = "Every 8 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out5c)),
  cbind(scenario = "Every 7 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out5b)),
  cbind(scenario = "Every 6 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out5a)),
  cbind(scenario = "Every 5 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out5)),
  cbind(scenario = "Every 4 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out6c)),
  cbind(scenario = "Every 3 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out6b)),
  cbind(scenario = "Every 2 years",
        assemble_discounted_interactions_for_monitoring_frequencies(out6a)),
  cbind(scenario = "Every 1 year",
        assemble_discounted_interactions_for_monitoring_frequencies(out6))
)

# Discount person-years on treatment
freq_interactions_py_on_treatment_disc <- rbind(
  data.frame(scenario = "No monitoring",
             discount_outcome_2020_to_2100(scenario_object=out3,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "At age 30",
             discount_outcome_2020_to_2100(scenario_object=out_lt_monit_30,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "At age 45",
             discount_outcome_2020_to_2100(scenario_object=out_lt_monit_45,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 30 years",
             discount_outcome_2020_to_2100(scenario_object=out4d,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 25 years",
             discount_outcome_2020_to_2100(scenario_object=out4c,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 20 years",
             discount_outcome_2020_to_2100(scenario_object=out4b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 15 years",
             discount_outcome_2020_to_2100(scenario_object=out4a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 10 years",
             discount_outcome_2020_to_2100(scenario_object=out4,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 9 years",
             discount_outcome_2020_to_2100(scenario_object=out5d,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 8 years",
             discount_outcome_2020_to_2100(scenario_object=out5c,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 7 years",
             discount_outcome_2020_to_2100(scenario_object=out5b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 6 years",
             discount_outcome_2020_to_2100(scenario_object=out5a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 5 years",
             discount_outcome_2020_to_2100(scenario_object=out5,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 4 years",
             discount_outcome_2020_to_2100(scenario_object=out6c,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 3 years",
             discount_outcome_2020_to_2100(scenario_object=out6b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 2 years",
             discount_outcome_2020_to_2100(scenario_object=out6a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "Every 1 year",
             discount_outcome_2020_to_2100(scenario_object=out6,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate))
)
freq_interactions_py_on_treatment_disc$sim <- gsub("[^0-9]", "", freq_interactions_py_on_treatment_disc$sim)


# Outcome 1: HBV related deaths averted
freq_hbv_deaths_averted_disc <- rbind(
  cbind(scenario = "No monitoring",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out3,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "At age 30",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out_lt_monit_30,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "At age 45",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out_lt_monit_45,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 30 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4d,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 25 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4c,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 20 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4b,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 15 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4a,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 10 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 9 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5d,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 8 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5c,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 7 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5b,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 6 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5a,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 5 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 4 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6c,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 3 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6b,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 2 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6a,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 1 year",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)))
freq_hbv_deaths_averted_disc$sim <- gsub("[^0-9]", "", freq_hbv_deaths_averted_disc$sim)
colnames(freq_hbv_deaths_averted_disc)[colnames(freq_hbv_deaths_averted_disc) == "cum_hbv_deaths"] <-
  "value"

# Outcome 2: DALYS averted
freq_dalys_averted_disc <- rbind(
  cbind(scenario = "No monitoring",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out3,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "At age 30",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out_lt_monit_30,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "At age 45",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out_lt_monit_45,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 30 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4d,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 25 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4c,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 20 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4b,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 15 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4a,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 10 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out4,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 9 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5d,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 8 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5c,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 7 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5b,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 6 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5a,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 5 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out5,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 4 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6c,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 3 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6b,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 2 years",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6a,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "Every 1 year",
        discount_outcome_2020_to_2100(scenario_object=out2,
                                      object_to_subtract=out6,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)))
freq_dalys_averted_disc$sim <- gsub("[^0-9]", "", freq_dalys_averted_disc$sim)
colnames(freq_dalys_averted_disc)[colnames(freq_dalys_averted_disc) == "dalys"] <-
  "value"

# Combine into full dataframes with costs ----

# No discounting
freq_df <- create_incremental_plot_df(interactions_df=freq_interactions,
                                      py_on_treatment_df=freq_interactions_py_on_treatment,
                                      deaths_averted_df=freq_hbv_deaths_averted,
                                      ly_saved_df = freq_dalys_averted, # replace LY by DALYs
                                      hbsag_test_cost = 8.3,
                                      clinical_assessment_cost = 84.4,
                                      monitoring_assessment_cost = 40.1,
                                      treatment_py_cost = 60,
                                      scenario_labels_obj = scenario_labels,
                                      ref_label = "No treatment")
colnames(freq_df)[colnames(freq_df)=="ly_saved"] <- "dalys_averted"

# Discount both outcomes and cost
freq_df_disc <- create_incremental_plot_df(interactions_df=freq_interactions_disc,
                                           py_on_treatment_df=freq_interactions_py_on_treatment_disc,
                                           deaths_averted_df=freq_hbv_deaths_averted_disc,
                                           ly_saved_df = freq_dalys_averted_disc, # replace LY by DALYs
                                           hbsag_test_cost = 8.3,
                                           clinical_assessment_cost = 84.4,
                                           monitoring_assessment_cost = 40.1,
                                           treatment_py_cost = 60,
                                           scenario_labels_obj = scenario_labels,
                                           ref_label = "No treatment")
colnames(freq_df_disc)[colnames(freq_df_disc)=="ly_saved"] <- "dalys_averted"

# Discount just costs
freq_df_disc_cost <- create_incremental_plot_df(interactions_df=freq_interactions_disc,
                                                py_on_treatment_df=freq_interactions_py_on_treatment_disc,
                                                deaths_averted_df=freq_hbv_deaths_averted,
                                                ly_saved_df = freq_dalys_averted, # replace LY by DALYs
                                                hbsag_test_cost = 8.3,
                                                clinical_assessment_cost = 84.4,
                                                monitoring_assessment_cost = 40.1,
                                                treatment_py_cost = 60,
                                                scenario_labels_obj = scenario_labels,
                                                ref_label = "No treatment")
colnames(freq_df_disc_cost)[colnames(freq_df_disc_cost)=="ly_saved"] <- "dalys_averted"

# Discount just outcomes
freq_df_disc_outcomes <- create_incremental_plot_df(interactions_df=freq_interactions,
                                                    py_on_treatment_df=freq_interactions_py_on_treatment,
                                                    deaths_averted_df=freq_hbv_deaths_averted_disc,
                                                    ly_saved_df = freq_dalys_averted_disc, # replace LY by DALYs
                                                    hbsag_test_cost = 8.3,
                                                    clinical_assessment_cost = 84.4,
                                                    monitoring_assessment_cost = 40.1,
                                                    treatment_py_cost = 60,
                                                    scenario_labels_obj = scenario_labels,
                                                    ref_label = "No treatment")
colnames(freq_df_disc_outcomes)[colnames(freq_df_disc_outcomes)=="ly_saved"] <- "dalys_averted"

# Analysis: calculate probability of each strategy being non-dominated and ICERs ----
# ICER is calculated for DALYs averted and cost

# A) No discounting ----

# Exclude strategies that have no clinical relevance (monitor > every 10/15 years)
freq_df2 <- subset(freq_df, !(scenario %in% c("Every 30 years", "Every 25 years",
                                              "Every 20 years", "Every 15 years", #"Every 10 years",
                                              "Every 9 years", "Every 8 years", "Every 7 years", "Every 6 years")))

dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- freq_df2[which(freq_df2$sim==
                                               unique(freq_df2$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Lifetime monitoring at age 30 or age 45 are dominated irrespective of
# whether every 15 years is included.
# If monitoring frequencies after lifetime monitoring start at every 5 years, at age 30 is
# not dominated.

# Calculate ICER by simulation
# Previous plot shows no strategies are dominated, so including all
freq_df2 <- subset(freq_df2, !(scenario %in% c("At age 45", "At age 30")))

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- freq_df2[which(freq_df2$sim==
                                     unique(freq_df2$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

# Including lifetime monitoring strategies, every 10 years and then yearly
# from every 5 years, every 10 years would be the most cost-effective strategy under the WTP
# on average and would dominate the lifetime monitoring strategies. Even
# if every 10 years if not considered feasible and would be excluded, 5-yearly monitoring
# would have low probability of being cost-effective, but in this case at age 30 would not
# be dominated and would represent the most cost-effective strategy under the WTP.

# Histograms of ICER
# Manually remove outer quantiles from every 1 year scenario with limits but need to do the same for all
ggplot(icer_df[icer_df$scenario != "Status quo",]) +
  geom_histogram(aes(icer, group = scenario, fill = scenario), bins = 75) +
  facet_wrap(~scenario, ncol = 3, scales = "free_y") +
  xlim(0,38953) +
  theme_bw()
# No monitoring y axis too large
# Shows that ICERS for the larger monitoring frequencies mainly get progressively more uncertain
# However the median is not sensitive to these outliers.

# Forest plot of ICER
# Could make sideways boxplot next to % non_dominated (Alistair thesis p227)
# Only 95 quantiles
ggplot(icer_result[icer_result$scenario != "Status quo",]) +
  geom_pointrange(aes(x = reorder(scenario, icer_median), y = icer_median, ymin = icer_lower, ymax =icer_upper),
                  size=1.2, position = position_dodge(width = 0.3), col = "turquoise") +
  ylim(0,38953) +
  theme_bw() +
  coord_flip()


# B) Discounting of costs and outcomes ----
dominance_prob_list <- list()

freq_df_disc2 <- freq_df_disc
# Exclude strategies that have no clinical relevance (monitor > every 10/15 years)
freq_df_disc2 <- subset(freq_df_disc, !(scenario %in% c("Every 30 years", "Every 25 years",
                                                        "Every 20 years", "Every 15 years", "Every 10 years",
                                                        "Every 9 years", "Every 8 years", "Every 7 years", "Every 6 years")))

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- freq_df_disc2[which(freq_df_disc2$sim==
                                                    unique(freq_df_disc2$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Calculate ICER by simulation
# Including all, lifetime monitoring strategies are dominated
# Excluding all regular frequencies > every 10 years, plus 6-9 years, lifetime monitoring strategies are dominated
# Excluding all regular frequencies > every 5 years, only at age 45 has < 50% prop_non_dominated
freq_df_disc2 <- subset(freq_df_disc2, !(scenario %in% c("At age 30", "At age 45")))
freq_df_disc2 <- subset(freq_df_disc2, !(scenario %in% c("At age 45")))

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- freq_df_disc2[which(freq_df_disc2$sim==
                                          unique(freq_df_disc2$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result_disc <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result_disc

# Cost-effectiveness acceptability curve for different scenarios
acceptability1 <- 0
acceptability2 <- 0
acceptability3 <- 0
acceptability4 <- 0
acceptability5 <- 0
scenario1 <- "No monitoring"
scenario2 <- "Every 30 years"
#scenario2 <- "Every 10 years"
#scenario2 <- "At age 30"
scenario3 <- "Every 15 years"
#scenario3 <- "Every 5 years"
scenario4 <- "Every 10 years"
#scenario4 <- "Every 4 years"
scenario5 <- "Every 5 years"

for (i in 1:752) {
  acceptability1[i] <-
    length(which(icer_df[icer_df$scenario==scenario1,]$icer<=i))/183
  acceptability2[i] <-
    length(which(icer_df[icer_df$scenario==scenario2,]$icer<=i))/183
  acceptability3[i] <-
    length(which(icer_df[icer_df$scenario==scenario3,]$icer<=i))/183
  acceptability4[i] <-
    length(which(icer_df[icer_df$scenario==scenario4,]$icer<=i))/183
  acceptability5[i] <-
    length(which(icer_df[icer_df$scenario==scenario5,]$icer<=i))/183
}
acceptability_curve <- data.frame(scenario = rep(c(scenario1,scenario2,
                                                   scenario3, scenario4,
                                                   scenario5),
                                                 each = 752),
                                  wtp = rep(1:752, times=5),
                                  prob_cost_effective = c(acceptability1,
                                                          acceptability2,
                                                          acceptability3, acceptability4,
                                                          acceptability5))

ggplot(acceptability_curve) +
  geom_line(aes(x=wtp, y = prob_cost_effective*100,
                group = reorder(scenario, -prob_cost_effective),
                colour = reorder(scenario, -prob_cost_effective))) +
  geom_vline(xintercept=391, linetype="dashed") +
  geom_vline(xintercept=518, linetype="dashed") +
  ylab("Probability cost-effective (%)") +
  xlab("Willingness to pay (USD/DALY averted)") +
  labs(colour = "Monitoring scenario") +
  theme_bw() +
  ylim(0,100) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# C) Discount only costs ----
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- freq_df_disc_cost[which(freq_df_disc_cost$sim==
                                                        unique(freq_df_disc_cost$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Calculate ICER by simulation
# Previous plot shows No monitoring is dominated with discounting of costs only

icer_list <- list()

freq_df_disc_cost2 <- subset(freq_df_disc_cost,
                             scenario != "No monitoring")


for(i in 1:183) {
  print(i)
  icer_list[[i]] <- freq_df_disc_cost2[which(freq_df_disc_cost2$sim==
                                               unique(freq_df_disc_cost2$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result_disc_cost <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result_disc_cost

# Analysis: calculate incremental impact (no costs, no discounting) ----

# DALYs averted per different interactions
freq_df$dalys_averted_per_monitoring <- freq_df$dalys_averted/freq_df$monitoring_assessments
freq_df$dalys_averted_per_py_on_treatment <- freq_df$dalys_averted/freq_df$py_on_treatment
freq_df$dalys_averted_per_treatment_initiation <- freq_df$dalys_averted/freq_df$treatment_initiations
freq_df$dalys_averted_per_hbsag_test <-freq_df$dalys_averted/freq_df$hbsag_tests
freq_df$dalys_averted_per_assessment <-freq_df$dalys_averted/freq_df$clinical_assessments

ggplot(subset(freq_df, scenario != "No treatment")) +
  geom_boxplot(aes(x=reorder(scenario,dalys_averted),
                   y = dalys_averted))+
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  ylab("DALYs averted compared to no treatment") +
  xlab("Monitoring frequency") +
  ylim(0,400000)

ggplot(subset(freq_df, scenario != "No treatment")) +
  geom_boxplot(aes(x=reorder(scenario,treatment_initiations),
                   y = treatment_initiations))+
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  ylab("DALYs averted compared to no treatment") +
  xlab("Monitoring frequency") +
  ylim(0,30000)

ggplot(subset(freq_df,scenario != "No monitoring" &
                scenario != "No treatment")) +
  geom_boxplot(aes(x=reorder(scenario,-dalys_averted_per_monitoring),
                   y = dalys_averted_per_monitoring*1000))+
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  ylab("DALYs averted per\n1,000 monitoring assessments") +
  xlab("Monitoring frequency") +
  ylim(0,7000)

ggplot(subset(freq_df,scenario != "No monitoring" &
                scenario != "No treatment")) +
  geom_boxplot(aes(x=reorder(scenario,-dalys_averted_per_treatment_initiation),
                   y = dalys_averted_per_treatment_initiation))+
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  ylab("DALYs averted per\nperson initiated on treatment") +
  xlab("Monitoring frequency") +
  ylim(0,21)


ggplot(subset(freq_df,scenario != "No monitoring" &
                scenario != "No treatment")) +
  geom_boxplot(aes(x=reorder(scenario,-dalys_averted_per_py_on_treatment),
                   y = dalys_averted_per_py_on_treatment*1000)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  ylab("DALYs averted per\n1,000 person-years on treatment") +
  xlab("Monitoring frequency") +
  ylim(0,800)

# Don't think this one is correct!

# Incremental DALYs averted for every 5 years monitoring interval
diff_list <- list()

freq_df_subs <- subset(freq_df, scenario %in% c("No monitoring",
                                                "Every 30 years",
                                                "At age 30", "At age 45",
                                                #"Every 25 years",
                                                #"Every 20 years",
                                                "Every 15 years",
                                                "Every 10 years","Every 5 years",
                                                "Every 1 year"))

for(i in 1:183) {
  print(i)
  diff_list[[i]] <- freq_df_subs[which(freq_df_subs$sim==
                                         unique(freq_df_subs$sim)[i]),]
  #  diff_list[[i]] <- calculate_icer_per_sim(diff_list[[i]],
  #                                           exposure="total_cost",
  #                                           outcome="dalys_averted")
  diff_list[[i]] <- arrange(diff_list[[i]], dalys_averted)
  diff_list[[i]]$diff_outcome <- c(diff_list[[i]][,"dalys_averted"][1],
                                   diff(diff_list[[i]][,"dalys_averted"]))
  diff_list[[i]]$comparator <- lag(diff_list[[i]]$scenario)


}
diff_df <- do.call("rbind", diff_list)

diff_result <- group_by(diff_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(median = median(dalys_averted),
            lower = quantile(dalys_averted, 0.025),
            upper = quantile(dalys_averted, 0.975)) %>%
  arrange(median)
diff_result

# Plot incremental DALYs averted between each 2 strategies
diff_df$scenario <- factor(as.character(diff_df$scenario),
                           levels = c("No monitoring", "At age 30", "Every 30 years",
                                      "At age 45", "Every 15 years", "Every 10 years",
                                      "Every 5 years", "Every 1 year"))

ggplot(subset(diff_df, scenario != "No monitoring"),
       aes(x= scenario, y = diff_outcome)) +
  geom_point(alpha=0.2) +
  #   geom_line(aes(x= scenario, y = diff_outcome, group = sim),alpha=0.05) +
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = "errorbar") +
  stat_summary(fun = median, geom = "point", col = "red", size = 3) +
  stat_summary(fun = median, colour="red", geom="line", aes(group = 1)) +
  ylab("Incremental DALYs averted compared to\nprevious monitoring frequency") +
  xlab("Monitoring frequency") +
  theme_bw()
# Shows that some monitoring frequencies are too infrequent to have much of an impact, as
# people will experience negative outcomes before they have a chance to be monitored and put
# on treatment (<15 years). However monitoring can also be too frequent (>5 years).
# At this point the return is smaller again.

ggplot(subset(freq_df_subs),
       aes(x= reorder(scenario, dalys_averted), y = dalys_averted)) +
  geom_point(alpha=0.2) +
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = "errorbar") +
  stat_summary(fun = median, geom = "point", col = "red", size = 3) +
  #  stat_summary(fun = median, colour="red", geom="line", aes(group = 1)) +
  ylab("DALYs averted in each scenario") +
  theme_bw() +
  ylim(0,361000)
# As before: monitored every 30/25/20 years averts almost the same number of DALYS, then it starts
# to increase a bit. On this plot you cannot see stagnation after 5 years as much.

# Zoom in on yearly increase from every 10 to every 1 year
freq_df_subs2 <- subset(freq_df, scenario %in% c("Every 15 years",
                                                 "Every 10 years","Every 9 years",
                                                 "Every 8 years","Every 7 years",
                                                 "Every 6 years","Every 5 years",
                                                 "Every 4 years",
                                                 "Every 3 years",
                                                 "Every 2 years",
                                                 "Every 1 year"))

for(i in 1:183) {
  print(i)
  diff_list[[i]] <- freq_df_subs2[which(freq_df_subs2$sim==
                                          unique(freq_df_subs2$sim)[i]),]
  diff_list[[i]] <- calculate_icer_per_sim(diff_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
diff_df2 <- do.call("rbind", diff_list)

diff_result2 <- group_by(diff_df2, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(median = median(dalys_averted),
            lower = quantile(dalys_averted, 0.025),
            upper = quantile(dalys_averted, 0.975)) %>%
  arrange(median)
diff_result2

# Plot incremental DALYs averted between each 2 strategies
diff_df2$scenario <- factor(as.character(diff_df2$scenario),
                            levels = c("Every 15 years", "Every 10 years",
                                       "Every 9 years",
                                       "Every 8 years","Every 7 years",
                                       "Every 6 years","Every 5 years",
                                       "Every 4 years",
                                       "Every 3 years",
                                       "Every 2 years",
                                       "Every 1 year"))

ggplot(subset(diff_df2, scenario != "Every 15 years" & scenario != "Every 10 years"),
       aes(x= scenario, y = diff_outcome)) +
  geom_point(alpha=0.2) +
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = "errorbar") +
  #geom_boxplot() +
  stat_summary(fun = median, geom = "point", col = "red", size = 3) +
  stat_summary(fun = median, colour="red", geom="line", aes(group = 1)) +
  ylab("Incremental DALYs averted compared to\nprevious monitoring frequency") +
  theme_bw()
# Shows that some monitoring frequencies are too infrequent to have much of an impact, as
# people will experience negative outcomes before they have a chance to be monitored and put
# on treatment (<15 years). However monitoring can also be too frequent (>5 years).
# At this point the return is smaller again.

# % HBV deaths averted in the cohort!
# This is the proportion of deaths averted that would have occurred WITHOUT treatment
# so in status quo, NOT compared to the previous strategy.

freq_df_cohort <- create_incremental_plot_df(interactions_df=freq_interactions,
                                             py_on_treatment_df=freq_interactions_py_on_treatment,
                                             deaths_averted_df=freq_hbv_deaths_averted_cohort,
                                             ly_saved_df = freq_dalys_averted_cohort, # replace LY by DALYs
                                             hbsag_test_cost = 8.3,
                                             clinical_assessment_cost = 84.4,
                                             monitoring_assessment_cost = 40.1,
                                             treatment_py_cost = 60,
                                             scenario_labels_obj = scenario_labels,
                                             ref_label = "No treatment")
colnames(freq_df_cohort)[colnames(freq_df_cohort)=="ly_saved"] <- "dalys_averted"

freq_df_cohort_subs <- subset(freq_df_cohort, scenario %in% c("No monitoring",
                                                              "Every 30 years", "Every 25 years",
                                                              "Every 20 years","Every 15 years",
                                                              "Every 10 years","Every 5 years",
                                                              "Every 1 year"))

# Doubling frequencies
freq_df_cohort_subs <- subset(freq_df_cohort, scenario %in% c("No monitoring",
                                                              "Every 30 years",
                                                              #"At age 30", "At age 45",
                                                              #"Every 25 years",
                                                              "Every 20 years",
                                                              #"Every 15 years",
                                                              "Every 10 years",
                                                              "Every 5 years",
                                                              "Every 2 years"))


diff_list <- list()
for(i in 1:183) {
  print(i)
  diff_list[[i]] <- freq_df_cohort_subs[which(freq_df_cohort_subs$sim==
                                                unique(freq_df_cohort_subs$sim)[i]),]
  diff_list[[i]] <- calculate_icer_per_sim(diff_list[[i]],
                                           exposure="total_cost",
                                           outcome="deaths_averted")
}
diff_df_cohort_deaths <- do.call("rbind", diff_list)

diff_result_cohort_deaths <- group_by(diff_df_cohort_deaths, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(median = median(diff_outcome),
            lower = quantile(diff_outcome, 0.025),
            upper = quantile(diff_outcome, 0.975)) %>%
  arrange(median)
diff_result_cohort_deaths

# Plot incremental DALYs averted between each 2 strategies
diff_df_cohort_deaths$scenario <- factor(as.character(diff_df_cohort_deaths$scenario),
                                         levels = c("No monitoring", "Every 30 years", "Every 25 years",
                                                    "Every 20 years", "Every 15 years", "Every 10 years",
                                                    "Every 5 years", "Every 1 year"))

diff_df_cohort_deaths$scenario <- factor(as.character(diff_df_cohort_deaths$scenario),
                                         levels = c("No monitoring",
                                                    "Every 30 years",
                                                    #"At age 30", "At age 45",
                                                    #"Every 25 years",
                                                    "Every 20 years",
                                                    #"Every 15 years",
                                                    "Every 10 years",
                                                    "Every 5 years",
                                                    "Every 2 years"))

# Order all:
diff_df_cohort_deaths$scenario <- factor(as.character(diff_df_cohort_deaths$scenario),
                                         levels = c("No monitoring", "At age 30",
                                                    "Every 30 years",
                                                    "At age 45",
                                                    "Every 25 years",
                                                    "Every 20 years",
                                                    "Every 15 years",
                                                    "Every 10 years",
                                                    "Every 9 years",
                                                    "Every 8 years",
                                                    "Every 7 years",
                                                    "Every 6 years",
                                                    "Every 5 years",
                                                    "Every 4 years",
                                                    "Every 3 years",
                                                    "Every 2 years",
                                                    "Every 1 year"))

ggplot(subset(diff_df_cohort_deaths, scenario != "No monitoring" & scenario != "Every 30 years"),
       aes(x= scenario, y = diff_outcome)) +
  geom_point(alpha=0.2) +
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = "errorbar") +
  stat_summary(fun = median, geom = "point", col = "red", size = 3) +
  stat_summary(fun = median, colour="red", geom="line", aes(group = 1)) +
  ylab("Incremental proportion of HBV deaths averted in cohort\ncompared to previous monitoring frequency") +
  theme_bw()

# Lookin at overall proportion averted:
ggplot(subset(freq_df_cohort_subs),
       aes(x= reorder(scenario, deaths_averted), y = deaths_averted)) +
  geom_point(alpha=0.2) +
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = "errorbar") +
  stat_summary(fun = median, geom = "point", col = "red", size = 3) +
  #  stat_summary(fun = median, colour="red", geom="line", aes(group = 1)) +
  ylab("Proportion of deaths averted in each scenario") +
  theme_bw()+
  ylim(0,1)

# Zoom in on yearly increase from every 10 to every 1 year
freq_df_cohort_subs2 <- subset(freq_df_cohort, scenario %in% c("Every 15 years",
                                                               "Every 10 years","Every 9 years",
                                                               "Every 8 years","Every 7 years",
                                                               "Every 6 years","Every 5 years",
                                                               "Every 4 years",
                                                               "Every 3 years",
                                                               "Every 2 years",
                                                               "Every 1 year"))

for(i in 1:183) {
  print(i)
  diff_list[[i]] <- freq_df_cohort_subs2[which(freq_df_cohort_subs2$sim==
                                                 unique(freq_df_cohort_subs2$sim)[i]),]
  diff_list[[i]] <- calculate_icer_per_sim(diff_list[[i]],
                                           exposure="total_cost",
                                           outcome="deaths_averted")
}
diff_df_cohort_deaths2 <- do.call("rbind", diff_list)

diff_result_cohort_deaths2 <- group_by(diff_df_cohort_deaths2, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(median = median(deaths_averted),
            lower = quantile(deaths_averted, 0.025),
            upper = quantile(deaths_averted, 0.975)) %>%
  arrange(median)
diff_result_cohort_deaths2

# Plot incremental DALYs averted between each 2 strategies
diff_df_cohort_deaths2$scenario <- factor(as.character(diff_df_cohort_deaths2$scenario),
                                          levels = c("Every 15 years", "Every 10 years",
                                                     "Every 9 years",
                                                     "Every 8 years","Every 7 years",
                                                     "Every 6 years","Every 5 years",
                                                     "Every 4 years",
                                                     "Every 3 years",
                                                     "Every 2 years",
                                                     "Every 1 year"))

ggplot(subset(diff_df_cohort_deaths2, scenario != "Every 15 years" & scenario != "Every 10 years"),
       aes(x= scenario, y = diff_outcome)) +
  geom_point(alpha=0.2) +
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = "errorbar") +
  stat_summary(fun = median, geom = "point", col = "red", size = 3) +
  stat_summary(fun = median, colour="red", geom="line", aes(group = 1)) +
  ylab("Incremental proportion of HBV deaths averted in cohort\ncompared to previous monitoring frequency") +
  theme_bw()
# Interpret as: monitoring every 9 years averts an extra 1% of deaths that would have occurred in the
# cohort without any treatment compared to monitoring every 10 years

# Forest plot ----
freq_df_daly_summary <- subset(freq_df, scenario %in% c("No monitoring",
                                                        "At age 30", "At age 45",
                                                        "Every 30 years",
                                                        "Every 25 years",
                                                        "Every 20 years", "Every 15 years", "Every 10 years",
                                                        "Every 5 years", "Every 1 year")) %>%
  group_by(scenario) %>%
  summarise(median = median(dalys_averted),
            lower = quantile(dalys_averted, 0.025),
            upper= quantile(dalys_averted, 0.975))

ggplot(freq_df_daly_summary) +
  geom_point(aes(x=reorder(scenario, median), y = median, colour = scenario), size = 3) +
  geom_linerange(aes(x=reorder(scenario, median), ymin = lower, ymax = upper, colour = scenario)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Monitoring\nfrequency") +
  ylab("DALYs averted") +
  ylim(0,360235) +
  coord_flip()

# Example function to find dominated strategies on median (but preferentially look at all sims) ----
# Including status quo
#out3 = 0, out4 = 10, out5d = 9, out5c = 8, out5b= 7, out5a = 6, out5 = 5,
#out6c = 4, out6b = 3, out6a = 2, out6 =1

# Find extended dominated strategies:
# Deaths averted and total interactions
ceef.plot.median(bcea(e=as.matrix(select(freq_df, scenario, sim, deaths_averted) %>%
                                    filter(scenario != "No treatment") %>%
                                    pivot_wider(names_from="scenario", values_from = "deaths_averted")%>%
                                    select(-sim)),
                      c=as.matrix(select(freq_df, scenario, sim, total_interactions) %>%
                                    filter(scenario != "No treatment") %>%
                                    pivot_wider(names_from="scenario", values_from = "total_interactions")%>%
                                    select(-sim)),
                      interventions=c(colnames(select(freq_df, scenario, sim, deaths_averted) %>%
                                                 filter(scenario != "No treatment") %>%
                                                 pivot_wider(names_from="scenario", values_from = "deaths_averted")%>%
                                                 select(-sim))),
                      Kmax=1000000000,
                      plot=FALSE),
                 graph="base", relative = FALSE)

# ICER Plots ----
# TO DO
# Add lines between 0 and the first strategy!
# Ideally here we want to plot the ICER calculated from individual sims as well

freq_df_median <- group_by(freq_df, scenario) %>%
  summarise(deaths_averted = median(deaths_averted),
            total_interactions = median(total_interactions),
            dalys_averted = median(dalys_averted),
            total_cost = median(total_cost))

# Deaths averted vs interactions
ggplot(data = freq_df) +
  geom_line(aes(x = deaths_averted, y = total_interactions, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group = reorder(scenario, total_interactions),
                 colour = reorder(scenario, total_interactions)), alpha = 0.15) +
  stat_ellipse(data = subset(freq_df, scenario != "No treatment"),
               aes(x = deaths_averted, y = total_interactions,
                   group = reorder(scenario, total_interactions), fill= reorder(scenario, total_interactions)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(freq_df_median, scenario != "No monitoring"),   # No monitoring looks to be dominated on average
            aes(x = deaths_averted, y = total_interactions), size = 1) +
  geom_point(data = freq_df_median,
             aes(x = deaths_averted, y = total_interactions,
                 group = reorder(scenario, total_interactions),
                 colour = reorder(scenario, total_interactions)),
             size = 5) +
  geom_point(data = freq_df_median,
             aes(x = deaths_averted, y = total_interactions,
                 group = reorder(scenario, total_interactions), colour = reorder(scenario, total_interactions)),
             size = 5, shape = 1, colour = "black") +
  # scale_fill_manual(values = rev(brewer.pal(13,"RdYlBu"))) +
  # scale_colour_manual("Monitoring frequency",
  #                      values = c("black", rev(brewer.pal(13,"RdYlBu")))) +
  guides(fill=FALSE) +
  xlab("Incremental total clinical interactions") +
  ylab("Incremental HBV-related deaths averted") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# Calculate median ICER compared to SQ for all strategies
freq_df %>%
  mutate(icer_sq = total_cost/dalys_averted) %>%
  group_by(scenario) %>%
  summarise(median_icer = median(icer_sq)) %>%
  arrange(median_icer)

# DALYS averted vs costs
ggplot(data = freq_df) +
  geom_line(aes(x = dalys_averted, y = total_cost, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), alpha = 0.15) +
  stat_ellipse(data = subset(freq_df, scenario != "No treatment"),
               aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost), fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = freq_df_median,
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = freq_df_median,
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost), colour = reorder(scenario, total_cost)),
             size = 5) +
  geom_point(data = freq_df_median,
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost), colour = reorder(scenario, total_cost)),
             size = 5, shape = 1, colour = "black") +
  #  scale_fill_manual(values = rev(brewer.pal(11,"RdYlBu"))) +
  #  scale_colour_manual("Monitoring frequency",
  #                      values = c("black", rev(brewer.pal(11,"RdYlBu")))) +
  guides(fill=FALSE) +
  geom_abline(slope=391, intercept = 0, linetype = "dashed") +
  geom_abline(slope=518, intercept = 0, linetype = "dashed") +
  xlab("Incremental total cost (USD 2019)") +
  ylab("Incremental DALYs averted") +
  labs(colour="Monitoring frequencies") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# DALYS averted vs costs
ggplot(data = subset(freq_df, scenario %in%
                       c("No treatment", "No monitoring", "Every 30 years", "Every 25 years",
                         "Every 20 years", "Every 15 years", "Every 10 years",
                         "Every 5 years", "Every 1 year"))) +
  geom_line(aes(x = dalys_averted, y = total_cost, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), alpha = 0.15) +
  stat_ellipse(data = subset(freq_df, scenario %in%
                               c("No monitoring", "Every 30 years", "Every 25 years",
                                 "Every 20 years", "Every 15 years", "Every 10 years",
                                 "Every 5 years", "Every 1 year")),
               aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost), fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(freq_df_median, scenario %in% c("No treatment", "No monitoring", "Every 30 years", "Every 25 years",
                                                          "Every 20 years", "Every 15 years", "Every 10 years",
                                                          "Every 5 years", "Every 1 year")),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = subset(freq_df_median, scenario %in% c("No treatment", "No monitoring", "Every 30 years", "Every 25 years",
                                                           "Every 20 years", "Every 15 years", "Every 10 years",
                                                           "Every 5 years", "Every 1 year")),
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost), colour = reorder(scenario, total_cost)),
             size = 5) +
  geom_point(data = subset(freq_df_median, scenario %in% c("No treatment", "No monitoring", "Every 30 years", "Every 25 years",
                                                           "Every 20 years", "Every 15 years", "Every 10 years",
                                                           "Every 5 years", "Every 1 year")),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost), colour = reorder(scenario, total_cost)),
             size = 5, shape = 1, colour = "black") +
  #  scale_fill_manual(values = rev(brewer.pal(11,"RdYlBu"))) +
  #  scale_colour_manual("Monitoring frequency",
  #                      values = c("black", rev(brewer.pal(11,"RdYlBu")))) +
  guides(fill=FALSE) +
  geom_abline(slope=391, intercept = 0, linetype = "dashed") +
  geom_abline(slope=518, intercept = 0, linetype = "dashed") +
  xlab("Incremental total cost (USD 2019)") +
  ylab("Incremental DALYs averted") +
  labs(colour="Monitoring frequencies") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# Plot with discounting
freq_df_disc_median <- group_by(freq_df_disc, scenario) %>%
  summarise(deaths_averted = median(deaths_averted),
            total_interactions = median(total_interactions),
            dalys_averted = median(dalys_averted),
            total_cost = median(total_cost))

ggplot(data = freq_df_disc) +
  geom_line(aes(x = dalys_averted, y = total_cost, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), alpha = 0.15) +
  stat_ellipse(data = subset(freq_df_disc, !(scenario %in%
                                               c("No treatment"))),
               aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost), fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(freq_df_disc_median, !(scenario %in% c("At age 30", "At age 45"))),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = freq_df_disc_median,
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost), colour = reorder(scenario, total_cost)),
             size = 5) +
  geom_point(data = freq_df_disc_median,
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost), colour = reorder(scenario, total_cost)),
             size = 5, shape = 1, colour = "black") +
  #  scale_fill_manual(values = rev(brewer.pal(11,"RdYlBu"))) +
  #  scale_colour_manual("Monitoring frequency",
  #                      values = c("black", rev(brewer.pal(11,"RdYlBu")))) +
  guides(fill=FALSE) +
  geom_abline(slope=391, intercept = 0, linetype = "dashed") +
  geom_abline(slope=518, intercept = 0, linetype = "dashed") +
  ylab("Incremental total cost (USD 2019)") +
  xlab("Incremental DALYs averted") +
  labs(colour="Monitoring frequencies") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))


# Minimal deaths averted vs interactions plot for IVHEM ----
ggplot(data = subset(incremental_df_freq, scenario %in% c("Status quo", "No monitoring",
                                                          "Every 10 years", "Every 5 years", "Every 3 years",
                                                          "Every 2 years",
                                                          "Every 1 year"))) +
  #  geom_line(aes(y = deaths_averted, x = interactions, group = sim),
  #            colour = "grey", alpha = 0.3) +
  #  geom_point(aes(y = deaths_averted, x= interactions,
  #                 group = scenario, colour = scenario), alpha = 0.15) +
  stat_ellipse(data = subset(incremental_df_freq, scenario %in% c("No monitoring",
                                                                  "Every 10 years", "Every 5 years", "Every 3 years",
                                                                  "Every 2 years",
                                                                  "Every 1 year")),
               aes(y=deaths_averted,x=interactions/1000,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(incremental_df_freq_summary, scenario %in% c("Status quo",
                                                                       "Every 10 years", "Every 5 years", "Every 3 years",
                                                                       "Every 2 years",
                                                                       "Every 1 year")),
            aes(x = median_interactions/1000,
                y = median_deaths_averted), size = 1) +
  geom_point(data = subset(incremental_df_freq_summary, scenario %in% c("Status quo", "No monitoring",
                                                                        "Every 10 years", "Every 5 years", "Every 3 years",
                                                                        "Every 2 years",
                                                                        "Every 1 year")),
             aes(x = median_interactions/1000,
                 y = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)),
             size = 7) +
  geom_point(data = subset(incremental_df_freq_summary, scenario %in% c("Status quo", "No monitoring",
                                                                        "Every 10 years", "Every 5 years", "Every 3 years",
                                                                        "Every 2 years",
                                                                        "Every 1 year")),
             aes(x = median_interactions/1000,
                 y = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 7, shape = 1, colour = "black") +
  scale_fill_manual(values = viridis(6, option = "D")) +  # rev(brewer.pal(6,"RdYlBu"))
  scale_colour_manual("Monitoring frequency",
                      values = c("black", viridis(6, option = "D")),
                      labels = c("Status quo" = "No screening and\ntreatment")) +
  guides(fill=FALSE) +
  xlab("Incremental total resources utilised\n(thousands)") +
  ylab("Incremental HBV-related deaths averted") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  coord_flip()

# Calculate probability of each strategy being non-dominated (no discounting) ----
# For LY saved and cost:
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df_ly_freq[which(incremental_df_ly_freq$sim==
                                                             unique(incremental_df_ly_freq$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="ly_saved")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Calculate ICER by simulation
# Previous plot shows no strategies are dominated, so including all

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- incremental_df_ly_freq[which(incremental_df_ly_freq$sim==
                                                   unique(incremental_df_ly_freq$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="ly_saved")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)

# Manually remove outer quantiles from every 1 year scenario with limits but need to do the same for all
ggplot(icer_df[icer_df$scenario != "Status quo",]) +
  geom_histogram(aes(icer, group = scenario, fill = scenario), bins = 75) +
  facet_wrap(~scenario, ncol = 3, scales = "free_y") +
  xlim(0,38953) +
  theme_bw()
# No monitoring y axis too large
# Shows that ICERS for the larger monitoring frequencies mainly get progressively more uncertain
# However the median is not sensitive to these outliers.

# Could make sideways boxplot next to % non_dominated (Alistair thesis p227)
# Only 95 quantiles
ggplot(icer_result[icer_result$scenario != "Status quo",]) +
  geom_pointrange(aes(x = reorder(scenario, icer_median), y = icer_median, ymin = icer_lower, ymax =icer_upper),
                  size=1.2, position = position_dodge(width = 0.3), col = "turquoise") +
  ylim(0,38953) +
  theme_bw() +
  coord_flip()

# For LY saved and interactions:
dominance_prob_list_int <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list_int[[i]] <- incremental_df_ly_freq[which(incremental_df_ly_freq$sim==
                                                                 unique(incremental_df_ly_freq$sim)[i]),]
  dominance_prob_list_int[[i]] <- assign_dominated_strategies(dominance_prob_list_int[[i]],
                                                              exposure="interactions",
                                                              outcome="ly_saved")
}
dominance_prob_int_df <- do.call("rbind", dominance_prob_list_int)
dominance_prob_int_result <- group_by(dominance_prob_int_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_int_result$prob_non_dominated <- dominance_prob_int_result$No/
  (dominance_prob_int_result$Yes+dominance_prob_int_result$No)

ggplot(dominance_prob_int_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Calculate "ICER" by simulation => Incremental interactions per effectiveness ratio
# All but no monitoring have >50% prob of being non-dominated => exclude

iier_list <- list()

for(i in 1:183) {
  print(i)
  iier_list[[i]] <- incremental_df_ly_freq[which(incremental_df_ly_freq$scenario != "No monitoring" &
                                                   incremental_df_ly_freq$sim==
                                                   unique(incremental_df_ly_freq$sim)[i]),]
  iier_list[[i]] <- calculate_icer_per_sim(iier_list[[i]],
                                           exposure="interactions",
                                           outcome="ly_saved")
}
iier_df <- do.call("rbind", iier_list)
iier_result <- group_by(iier_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(iier_median = median(icer),
            iier_lower = quantile(icer, 0.025),
            iier_upper = quantile(icer, 0.975)) %>%
  arrange(iier_median)

ggplot(iier_result[iier_result$scenario != "Status quo",]) +
  geom_pointrange(aes(x = reorder(scenario, iier_median), y = iier_median, ymin = iier_lower, ymax =iier_upper),
                  size=1.2, position = position_dodge(width = 0.3), col = "turquoise") +
  ylim(0,1075) +
  theme_bw() +
  coord_flip()

## 2) Compare cost-effectiveness of adding IT as a category to treat ----
interactions_test_it <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        left_join(left_join(left_join(gather(out3n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3n$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_5",
        left_join(left_join(left_join(gather(out5n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5n$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_2",
        left_join(left_join(left_join(gather(out6an$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out6an$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out6an$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_0_it",
        left_join(left_join(left_join(gather(out3_it$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3_it$interactions[[16]]$total_assessed[-c(1:3)]-out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_5_it",
        left_join(left_join(left_join(gather(out5_it$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3_it$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5_it$interactions[[16]]$total_assessed[-c(1:3)]-out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_2_it",
        left_join(left_join(left_join(gather(out6a_it$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3_it$interactions[[16]]$total_assessed[-c(1:3)],   # stays the same
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out6a_it$interactions[[16]]$total_assessed[-c(1:3)]-out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out6a_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")))

ggplot(interactions_test_it) +
  geom_boxplot(aes(x=scenario, y = monitoring_assessments))
ggplot(interactions_test_it) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations))
# Monitoring assessments almost the same but treatment initiations are substantially
# higher if treating IT!

# Extract person-years on treatment
py_on_treatment_test_it <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             sim = names(out3n$py_on_treatment[[16]]),
             py_on_treatment = out3n$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_5",
             sim = names(out5n$py_on_treatment[[16]]),
             py_on_treatment = out5n$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_2",
             sim = names(out6an$py_on_treatment[[16]]),
             py_on_treatment = out6an$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_0_it",
             sim = names(out3_it$py_on_treatment[[16]]),
             py_on_treatment = out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_5_it",
             sim = names(out5_it$py_on_treatment[[16]]),
             py_on_treatment = out5_it$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_2_it",
             sim = names(out6a_it$py_on_treatment[[16]]),
             py_on_treatment = out6a_it$py_on_treatment[[16]])
)
py_on_treatment_test_it$sim <- gsub("[^0-9]", "", py_on_treatment_test_it$sim)

ggplot(py_on_treatment_test_it) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment))
# PY on treatment also substantially higher if treating IT

dalys_averted_no_it <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(out3n,out5n,out6an),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted_with_it <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(out3_it,out5_it,out6a_it),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted_with_it$scenario <- paste0(dalys_averted_with_it$scenario, "_it")
dalys_averted_test_it <- rbind(dalys_averted_no_it,
                               dalys_averted_with_it)
dalys_averted_test_it$sim <- gsub("[^0-9]", "", dalys_averted_test_it$sim)
colnames(dalys_averted_test_it)[colnames(dalys_averted_test_it) == "dalys"] <-
  "value"

ggplot(subset(dalys_averted_test_it, by_year==2100 & type == "proportion_averted")) +
  geom_boxplot(aes(x=scenario, y = value))
# Only tiny bit higher with IT treatment

# HBV-related deaths
deaths_averted_no_it <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(out3n,out5n,out6an),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
deaths_averted_with_it <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(out3_it,out5_it,out6a_it),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
deaths_averted_with_it$scenario <- paste0(deaths_averted_with_it$scenario, "_it")
deaths_averted_test_it <- rbind(deaths_averted_no_it,
                                deaths_averted_with_it)
deaths_averted_test_it$sim <- gsub("[^0-9]", "", deaths_averted_test_it$sim)
colnames(deaths_averted_test_it)[colnames(deaths_averted_test_it) == "cum_hbv_deaths"] <-
  "value"

ggplot(subset(deaths_averted_test_it, by_year==2100 & type == "proportion_averted")) +
  geom_boxplot(aes(x=scenario, y = value))
# Same for deaths

deaths_averted_test_it<- subset(deaths_averted_test_it, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
dalys_averted_test_it<- subset(dalys_averted_test_it, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)

# Combine into dataframe
test_it_df <- create_incremental_plot_df(interactions_df=interactions_test_it,
                                      py_on_treatment_df=py_on_treatment_test_it,
                                      deaths_averted_df=deaths_averted_test_it,
                                      ly_saved_df = dalys_averted_test_it, # replace LY by DALYs
                                      hbsag_test_cost = 8.3,
                                      clinical_assessment_cost = 84.4,
                                      monitoring_assessment_cost = 40.1,
                                      treatment_py_cost = 60,
                                      #scenario_labels_obj = scenario_labels,
                                      ref_label = "No treatment")
colnames(test_it_df)[colnames(test_it_df)=="ly_saved"] <- "dalys_averted"

dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- test_it_df[which(test_it_df$sim==
                                               unique(test_it_df$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Dominated strategies are screen_2020_monit_0_it and screen_2020_monit_2
# screen_2020_monit_5 has 25% probability of being dominated

# That means treating IT becomes better if there is also monitoring - this makes sense

# Calculate ICER by simulation
# Previous plot shows no strategies are dominated, so including all
test_it_df2 <- subset(test_it_df, !(scenario %in% c("screen_2020_monit_0_it",
                                                    "screen_2020_monit_2")))

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- test_it_df2[which(test_it_df2$sim==
                                     unique(test_it_df2$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result
# Most cost-effective is no monitoring, no IT treatment followed by 5-yearly monitoring
# and no IT treatment.
# Addition of IT treatment to that is then probably not cost-effective anymore.

# Changed this to corrct count of monitoring when there is no IT treatment
# but this did not change conclusion!

# That IT treatment without monitoring is dominated is a very clear result (in almost all
# simulations). Actually the % of simulations cost-effective is also quite
# considerably lower when adding IT treatment to 5-yearly monitoring.

# Conclusion: Cost-effectiveness of treating IT already depends on how often you
# monitor, and if we want to keep it general with regards to treatment criteria
# we may not want to present this analysis explicitly.

test_it_df_median <- group_by(test_it_df, scenario) %>%
  summarise(deaths_averted = median(deaths_averted),
            total_interactions = median(total_interactions),
            dalys_averted = median(dalys_averted),
            total_cost = median(total_cost))

# DALYS averted vs costs
# Add different shapes for dominated vs non-dominated
ggplot(data = test_it_df) +
  #  geom_line(aes(x = dalys_averted, y = total_cost, group = sim),
  #            colour = "grey", alpha = 0.3) +
  geom_point(aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), alpha = 0.15) +
  stat_ellipse(data = subset(test_it_df, scenario != "No treatment"),
               aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(test_it_df_median, !(scenario %in% c("screen_2020_monit_0_it",
                                                             "screen_2020_monit_2"))),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = test_it_df_median,
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5) +
  geom_point(data = test_it_df_median,
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual(values = rev(brewer.pal(6,"RdYlBu"))) +
  scale_colour_manual("Monitoring frequency",
                      values = c("black", rev(brewer.pal(6,"RdYlBu")))) +
  guides(fill=FALSE) +
  geom_abline(slope=391, intercept = 0, linetype = "dashed") +
  geom_abline(slope=518, intercept = 0, linetype = "dashed") +
  xlab("Incremental total cost (USD 2019)") +
  ylab("Incremental DALYs averted") +
  labs(colour="Monitoring frequencies") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# Cost-effectiveness acceptability curve for different scenarios
acceptability1 <- 0
acceptability2 <- 0
acceptability3 <- 0
acceptability4 <- 0
scenario1 <- "screen_2020_monit_0"
scenario2 <- "screen_2020_monit_5"
scenario3 <- "screen_2020_monit_5_it"
scenario4 <- "screen_2020_monit_2_it"

for (i in 1:752) {
  acceptability1[i] <-
    length(which(icer_df[icer_df$scenario==scenario1,]$icer<=i))/183
  acceptability2[i] <-
    length(which(icer_df[icer_df$scenario==scenario2,]$icer<=i))/183
  acceptability3[i] <-
    length(which(icer_df[icer_df$scenario==scenario3,]$icer<=i))/183
  acceptability4[i] <-
    length(which(icer_df[icer_df$scenario==scenario4,]$icer<=i))/183

}
acceptability_curve <- data.frame(scenario = rep(c(scenario1,scenario2,
                                                   scenario3, scenario4),
                                                 each = 752),
                                  wtp = rep(1:752, times=4),
                                  prob_cost_effective = c(acceptability1,
                                                          acceptability2,
                                                          acceptability3,
                                                          acceptability4))

ggplot(acceptability_curve) +
  geom_line(aes(x=wtp, y = prob_cost_effective*100,
                group = reorder(scenario, -prob_cost_effective),
                colour = reorder(scenario, -prob_cost_effective))) +
  geom_vline(xintercept=391, linetype="dashed") +
  geom_vline(xintercept=518, linetype="dashed") +
  ylab("Probability cost-effective (%)") +
  xlab("Willingness to pay (USD/DALY averted)") +
  labs(colour = "Monitoring scenario") +
  theme_bw() +
  ylim(0,100) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

## 3) Cost-effectiveness of different monitoring frequencies by age (with IT treatment) ----
# Create data frame with all interactions and outcomes of interest (discounting switch) ----

# Need to load function
annual_discounting_rate <- 0.03  # can be 0

# Newly added objects:
# out4_it, out6c_it, out6b_it, out6_it
# monit_out4c, monit_out5c, monit_out2c, monit_out1c
# monit_out4b, monit_out2b
# monit_out4, monit_out5, monit_out2, monit_out1,
# monit_out5b, monit_out1b, monit_out7_10

# Objects to include in analysis
# For interactions in this case, the no_monitoring_object is out3_it for all of these
# And there is no object_to_substract for calculation of person years of treatment
# For HBV deaths and DALYS averted, the scenario_object (reference scenario)
# is always out2
object_list <- list(out3_it, out4_it, out5_it, out6c_it, out6b_it,
                    out6a_it, out6_it, monit_out7, monit_out17, monit_out2a,
                 monit_out9, monit_out18, monit_out4a, monit_out10, monit_out5a,
                 monit_out6,monit_out1a, monit_out8, monit_out7_10,
                 monit_out4c, monit_out5c, monit_out2c, monit_out1c,
                 monit_out4b, monit_out2b,monit_out5b, monit_out1b,
                 monit_out4, monit_out5, monit_out2, monit_out1,
                 out_it_lt_monit_30,out_it_lt_monit_45, out_it_lt_monit_30_45)

# Extract interactions and person-years on treatment, and HBV-related deaths
# ad DALYs averted, in a for loop
age_interactions <- list()
age_interactions_py_on_treatment <- list()
age_hbv_deaths_averted <- list()
age_dalys_averted <- list()

for (i in 1:length(object_list)) {
  age_interactions[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          assemble_discounted_interactions_for_monitoring_frequencies(object_list[[i]],
          no_monitoring_object = out3_it))
  age_interactions_py_on_treatment[[i]] <-
    data.frame(scenario = object_list[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate))
  age_hbv_deaths_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="cum_hbv_deaths",
                                        yearly_discount_rate=annual_discounting_rate))
  age_dalys_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="dalys",
                                        yearly_discount_rate=annual_discounting_rate))
}
age_interactions <- do.call("rbind", age_interactions)
age_interactions_py_on_treatment <- do.call("rbind", age_interactions_py_on_treatment)
age_hbv_deaths_averted <- do.call("rbind", age_hbv_deaths_averted)
age_dalys_averted <- do.call("rbind", age_dalys_averted)

age_interactions_py_on_treatment$sim <- gsub("[^0-9]", "", age_interactions_py_on_treatment$sim)
age_hbv_deaths_averted$sim <- gsub("[^0-9]", "", age_hbv_deaths_averted$sim)
colnames(age_hbv_deaths_averted)[colnames(age_hbv_deaths_averted) == "cum_hbv_deaths"] <-
  "value"

age_dalys_averted$sim <- gsub("[^0-9]", "", age_dalys_averted$sim)
colnames(age_dalys_averted)[colnames(age_dalys_averted) == "dalys"] <-
  "value"

# Combine into full dataframe with costs ----
age_df <- create_incremental_plot_df(interactions_df=age_interactions,
                                      py_on_treatment_df=age_interactions_py_on_treatment,
                                      deaths_averted_df=age_hbv_deaths_averted,
                                      ly_saved_df = age_dalys_averted, # replace LY by DALYs
                                      hbsag_test_cost = 8.3,
                                      clinical_assessment_cost = 84.4,
                                      monitoring_assessment_cost = 40.1,
                                      treatment_py_cost = 60,
                                      #scenario_labels_obj = scenario_labels,
                                      ref_label = "No treatment")
colnames(age_df)[colnames(age_df)=="ly_saved"] <- "dalys_averted"

# Lifetime monitoring simulations are wrong - drop these
age_df <- subset(age_df, !(scenario %in% c("screen_2020_monit_lifetime_30",
                                           "screen_2020_monit_lifetime_45",
                                           "screen_2020_monit_lifetime_30_45")))


# Analysis: calculate probability of each strategy being non-dominated and ICERs ----

dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- age_df[which(age_df$sim==
                                               unique(age_df$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust =1))

# Find dominated strategies (<50% non-dominated) with 3% discounting and monit_sim7_10 included:
# Only non-dominated strategies are:
# screen_2020_monit_1, screen_2020_monit_2, screen_2020_monit_3,
# screen_2020_monit_4, screen_2020_monit_5,
# screen_2020_monit_sim7,screen_2020_monit_sim7_10,
# screen_2020_monit_0, screen_2020_monit_10
# The latter 2 have >25% probability of being dominated
# With monit_sim7_10 excluded:
# Add screen_2020_monit_sim8, screen_2020_monit_sim6

# Previously looking only at 5 and 2-yearly frequencies:
# screen_2020_monit_sim7, screen_2020_monit_2,
# screen_2020_monit_5, screen_2020_monit_sim8,
# screen_2020_monit_sim6, and screen_2020_monit_0
# screen_2020_monit_0 has > 35% probability of being dominated though.
# Clear distinction among the others

# Note: if looking at deaths averted instead of DALYs averted, the non-dominated strategies are:
# screen_2020_monit_5, screen_2020_monit_2, and screen_2020_monit_0 (same as DALYS),
# screen_2020_monit_sim9, screen_2020_monit_sim10
# Whereas sim7 and 6 have a pretty high probability of being dominated.
# Sim9 and sim10 are focused on the older (30+/45+ ages)!


# Calculate ICER by simulation on non-dominated strategies
age_df2 <- subset(age_df, scenario %in% c("screen_2020_monit_0", "screen_2020_monit_1",
                                          "screen_2020_monit_2", "screen_2020_monit_3",
                                          "screen_2020_monit_4", "screen_2020_monit_5",
                                          "screen_2020_monit_10",
                                          "screen_2020_monit_sim7",
                                          "screen_2020_monit_sim7_10"))

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- age_df2[which(age_df2$sim==
                                     unique(age_df2$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

# Monitoring every 5 years among the younger age groups is more cost-effective than
# every 10 years across all ages. Note the result of cost-effectiveness of 5-yearly
# monitoring in all age groups depends on inclusion 10-yearly frequency (without it might just be,
# with 10-yearly included it definitely is not)
# If 10-yearly monitoring among <45 year olds is included as an option, then this is
# the most cost-effective option and th 5-yearly frequency in this age group is just not
# cost-effective anymore on the median (ICER of 526).But still sufficient to
# potentially argue for this since it is more feasible in reality.

# Monitoring every 10 years in 15-45 year olds would correspond at most to 3 monitoring
# assessments (in 15 year olds)
# Need to simulate the population to monitor to count this exacty (essentially the
# screened compartments or the cohort - the treated people)
# Cohort size - initial treatment initiations
# Approximated pop to monitor:
(unlist(monit_out7_10$cohort_size[,-1])-
    age_df[age_df$scenario=="screen_2020_monit_0",]$treatment_initiations)

quantile(age_df[age_df$scenario=="screen_2020_monit_sim7_10",]$monitoring_assessments/
  (unlist(monit_out7_10$cohort_size[,-1])-
     age_df[age_df$scenario=="screen_2020_monit_0",]$treatment_initiations),
  c(0.5,0.025,0.975))
# 10-yearly monitoring would correspond to less than 1 monitoring assessment per person
# at 80% monitoring prob (but this approximation does not account for age groups being treated)

# Original persons to screen (15-45 year olds at entry into cohort):

# Monitoring assessments per person:
# This is monitoring assessments per number of persons originally eligible for that
# monitoring (all 15-45 year olds at entry into cohort in principle, though technically
# one could consider only <36 year olds being at risk for the 10-yearly)
quantile(age_df[age_df$scenario=="screen_2020_monit_sim7_10",]$monitoring_assessments/
  rowSums(out3_it$cohort_size_screened[,which(ages==15):which(ages==45-da)]),
  c(0.5,0.025,0.975))
# 10-yearly monitoring would correspond to less than 1 monitoring assessment per person
# at 80% monitoring prob
quantile(age_df[age_df$scenario=="screen_2020_monit_sim7_10",]$monitoring_assessments/
           rowSums(out3_it$cohort_size_screened[,which(ages==15):which(ages==36-da)]),
         c(0.5,0.025,0.975))
# If excluding >35 year olds it is 1.3 per person
quantile(age_df[age_df$scenario=="screen_2020_monit_sim7",]$monitoring_assessments/
           rowSums(out3_it$cohort_size_screened[,which(ages==15):which(ages==45-da)]),
         c(0.5,0.025,0.975))
# About 1.55 per person at 80% probability for the 5-yearly frequency in the same age group

# ICER plots ----
# Need to do this on a subset

age_df$frontier <- "Dominated"
age_df$frontier[age_df$scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                "screen_2020_monit_2", "screen_2020_monit_sim7",
                                "screen_2020_monit_sim8",
                                "screen_2020_monit_sim6", "No treatment")] <- "Non-dominated"



age_df_median <- group_by(age_df, scenario, frontier) %>%
  summarise(deaths_averted = median(deaths_averted),
            total_interactions = median(total_interactions),
            dalys_averted = median(dalys_averted),
            total_cost = median(total_cost))

# DALYS averted vs costs
ggplot(data = age_df) +
#  geom_line(aes(x = dalys_averted, y = total_cost, group = sim),
#            colour = "grey", alpha = 0.3) +
#  geom_point(aes(x = dalys_averted, y = total_cost,
#                 group = reorder(scenario, total_cost),
#                 colour = reorder(scenario, total_cost)), alpha = 0.15) +
  stat_ellipse(data = subset(age_df, scenario != "No treatment"),
               aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(age_df_median, frontier == "Non-dominated" |
                            scenario == "No treatment"),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = age_df_median,
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost),
                 shape = frontier),
             size = 5) +
  geom_point(data = subset(age_df_median, frontier == "Dominated"),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 1,
             colour = "black") +
  geom_point(data = subset(age_df_median, frontier == "Non-dominated"),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 2,
             colour = "black") +
  #scale_fill_manual(values = rev(brewer.pal(11,"RdYlBu"))) +
  scale_fill_manual(values = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))) +
  scale_colour_manual("Monitoring strategy\n[age group (frequency)]",
                      #values = c("black", rev(brewer.pal(14,"RdYlBu"))),
                      values = c("black",  brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6)),
                      labels = c("screen_2020_monit_0" = "15+ (0)",
                                 "screen_2020_monit_5" = "15+ (5)",
                                 "screen_2020_monit_2" = "15+ (2)",
                                 "screen_2020_monit_sim7" = "15-45 (5), 45+ (0)",
                                 "screen_2020_monit_sim17" = "15-30 (2), 30-45 (5), 45+ (0)",
                                 "screen_2020_monit_sim2a" = "15-45 (2), 45+ (0)",
                                 "screen_2020_monit_sim9" = "15-30 (0), 30+ (5)",
                                 "screen_2020_monit_sim18" = "15-30 (0), 30-45 (5), 45+ (2)",
                                 "screen_2020_monit_sim4a" = "15-30 (0), 30+ (2)",
                                 "screen_2020_monit_sim10" = "15-45 (0), 45+ (5)",
                                 "screen_2020_monit_sim5a" = "15-45 (0), 45+ (2)",
                                 "screen_2020_monit_sim6" = "15-30 (5), 30+ (0)",
                                 "screen_2020_monit_sim1a" = "15-30 (2), 30+ (0)",
                                 "screen_2020_monit_sim8" = "15-30 (0), 30-45 (5), 45+ (0)",
                                 "No treatment" = "No treatment")) +
  guides(fill=FALSE) +
  guides(shape=FALSE) +
#  geom_abline(slope=391, intercept = 0, linetype = "dashed") +
#  geom_abline(slope=518, intercept = 0, linetype = "dashed") +
  ylab("Incremental total cost (USD 2019)") +
  xlab("Incremental DALYs averted") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# Notation: age group (frequency)

# Plot with only non-dominated strategies
ggplot(data = subset(age_df, frontier == "Non-dominated" | scenario == "No treatment")) +
  geom_line(aes(x = dalys_averted, y = total_cost, group = sim),
              colour = "grey", alpha = 0.2) +
  geom_segment(aes(x = 0, y = 0, xend = dalys_averted, yend = total_cost),
              data = subset(age_df, scenario == "screen_2020_monit_0"),
              colour = "grey", alpha = 0.2) +
  geom_point(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   colour = reorder(scenario, total_cost)), alpha = 0.25) +
  stat_ellipse(data = subset(age_df, frontier == "Non-dominated" & scenario != "No treatment"),
               aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.1) +
  geom_line(data = subset(age_df_median, frontier == "Non-dominated" |
                            scenario == "No treatment"),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = subset(age_df_median, frontier == "Non-dominated" |
                             scenario == "No treatment"),
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             shape = 17,
             size = 5) +
  geom_point(data = subset(age_df_median, frontier == "Non-dominated" |
                             scenario == "No treatment"),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 2,
             colour = "black") +
  #scale_fill_manual(values = rev(brewer.pal(11,"RdYlBu"))) +
  scale_fill_manual(values = brewer.pal(name="Dark2", n = 6)) +
  scale_colour_manual("Monitoring strategy\n[age group (frequency)]",
                      #values = c("black", rev(brewer.pal(14,"RdYlBu"))),
                      values = c("black",  brewer.pal(name="Dark2", n = 6)),
                      labels = c("screen_2020_monit_0" = "15+ (0)",
                                 "screen_2020_monit_5" = "15+ (5)",
                                 "screen_2020_monit_2" = "15+ (2)",
                                 "screen_2020_monit_sim7" = "15-45 (5), 45+ (0)",
                                 "screen_2020_monit_sim17" = "15-30 (2), 30-45 (5), 45+ (0)",
                                 "screen_2020_monit_sim2a" = "15-45 (2), 45+ (0)",
                                 "screen_2020_monit_sim9" = "15-30 (0), 30+ (5)",
                                 "screen_2020_monit_sim18" = "15-30 (0), 30-45 (5), 45+ (2)",
                                 "screen_2020_monit_sim4a" = "15-30 (0), 30+ (2)",
                                 "screen_2020_monit_sim10" = "15-45 (0), 45+ (5)",
                                 "screen_2020_monit_sim5a" = "15-45 (0), 45+ (2)",
                                 "screen_2020_monit_sim6" = "15-30 (5), 30+ (0)",
                                 "screen_2020_monit_sim1a" = "15-30 (2), 30+ (0)",
                                 "screen_2020_monit_sim8" = "15-30 (0), 30-45 (5), 45+ (0)",
                                 "No treatment" = "No treatment")) +
  guides(fill=FALSE) +
  #  geom_abline(slope=391, intercept = 0, linetype = "dashed") +
  #  geom_abline(slope=518, intercept = 0, linetype = "dashed") +
  ylab("Incremental total cost (USD 2019)") +
  xlab("Incremental DALYs averted") +
  labs(colour="Monitoring frequencies") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# Notation: age group (frequency)







# Cost-acceptability curve to do ----

# For 1 screen_2020_monit_0,screen_2020_monit_sim7_10,screen_2020_monit_sim7,
# screen_2020_monit_10  and screen_2020_monit_5
# if 10-yearly frequencies are included
acceptability1 <- 0
acceptability2 <- 0
acceptability3 <- 0
acceptability4 <- 0
acceptability5 <- 0
scenario1 <- "screen_2020_monit_0"
scenario2 <- "screen_2020_monit_sim7_10"
scenario3 <- "screen_2020_monit_sim7"
scenario4 <- "screen_2020_monit_10"
scenario5 <- "screen_2020_monit_5"

for (i in 1:752) {   # 725 is the max WTP
  acceptability1[i] <-
    length(which(icer_df[icer_df$scenario==scenario1,]$icer<=i&
                   icer_df[icer_df$scenario==scenario1,]$icer>=0))/183
  acceptability2[i] <-
    length(which(icer_df[icer_df$scenario==scenario2,]$icer<=i&
                   icer_df[icer_df$scenario==scenario2,]$icer>=0))/183
  acceptability3[i] <-
    length(which(icer_df[icer_df$scenario==scenario3,]$icer<=i&
                   icer_df[icer_df$scenario==scenario3,]$icer>=0))/183
  acceptability4[i] <-
    length(which(icer_df[icer_df$scenario==scenario4,]$icer<=i&
                   icer_df[icer_df$scenario==scenario4,]$icer>=0))/183
  acceptability5[i] <-
    length(which(icer_df[icer_df$scenario==scenario5,]$icer<=i&
                   icer_df[icer_df$scenario==scenario5,]$icer>=0))/183
}
acceptability_curve <- data.frame(scenario = rep(c(scenario1,scenario2,
                                                   scenario3, scenario4,
                                                   scenario5),
                                                 each = 752),
                                  wtp = rep(1:752, times=5),
                                  prob_cost_effective = c(acceptability1,
                                                          acceptability2,
                                                          acceptability3,
                                                          acceptability4,
                                                          acceptability5))

ggplot(acceptability_curve) +
  geom_line(aes(x=wtp, y = prob_cost_effective*100,
                group = reorder(scenario, -prob_cost_effective),
                colour = reorder(scenario, -prob_cost_effective))) +
  geom_vline(xintercept=391, linetype="dashed") +
  geom_vline(xintercept=518, linetype="dashed") +
  ylab("Probability cost-effective (%)") +
  xlab("Willingness to pay (USD/DALY averted)") +
  labs(colour = "Monitoring scenario") +
  theme_bw() +
  #ylim(0,100) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
## 4) Plots of breakdown of impact of 5-yearly monitoring and screening by separate age group ----

# EFFECT OF MONITORING BY AGE

# DALYs averted by monitoring
dalys_averted_by_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out3_it,
                                 scenario_objects = list(monit_out6,
                                                         monit_out8,
                                                         monit_out10),
                                 outcome_to_avert = "dalys",
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
dalys_averted_by_age_group <- subset(dalys_averted_by_age_group, type == "number_averted" &
                                       by_year==2100) %>%
  select(scenario, sim, value)
dalys_averted_by_age_group$sim <- gsub("[^0-9]", "", dalys_averted_by_age_group$sim)
colnames(dalys_averted_by_age_group)[3] <- "dalys_averted"

deaths_averted_by_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out3_it,
                          scenario_objects = list(monit_out6,
                                                  monit_out8,
                                                  monit_out10),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "treatment programme without monitoring")
deaths_averted_by_age_group <- subset(deaths_averted_by_age_group, type == "number_averted" &
                                       by_year==2100) %>%
  select(scenario, sim, value)
# ISSUE WITH SIM NAMES HERE, USE THOSE FROM DALYS
#deaths_averted_by_age_group$sim <- gsub("[^0-9]", "", deaths_averted_by_age_group$sim)
deaths_averted_by_age_group$sim <-dalys_averted_by_age_group$sim
colnames(deaths_averted_by_age_group)[3] <- "deaths_averted"

# In terms of interactions, want the extra monitoring assessments, treatment initiations and
# person years of treatment due to monitoring

interactions_by_age_group <- rbind(
  cbind(scenario = "screen_2020_monit_sim6",
        left_join(gather(monit_out6$interactions[[16]]$total_assessed[-c(1:3)]-
                                     out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"),
                  gather(monit_out6$interactions[[16]]$total_treated[-c(1:3)]-
                           out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim8",
        left_join(gather(monit_out8$interactions[[16]]$total_assessed[-c(1:3)]-
                           out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out8$interactions[[16]]$total_treated[-c(1:3)]-
                           out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim10",
        left_join(gather(monit_out10$interactions[[16]]$total_assessed[-c(1:3)]-
                           out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out10$interactions[[16]]$total_treated[-c(1:3)]-
                           out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")))
interactions_by_age_group$sim <- gsub("[^0-9]", "", interactions_by_age_group$sim)

py_on_treatment_by_age_group <-rbind(
  data.frame(scenario = "screen_2020_monit_sim6",
             sim = names(monit_out6$py_on_treatment[[16]]),
             py_on_treatment = monit_out6$py_on_treatment[[16]]-
               out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim8",
             sim = names(monit_out8$py_on_treatment[[16]]),
             py_on_treatment = monit_out8$py_on_treatment[[16]]-
               out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim10",
             sim = names(monit_out10$py_on_treatment[[16]]),
             py_on_treatment = monit_out10$py_on_treatment[[16]]-
               out3_it$py_on_treatment[[16]]))

df_by_age_group <- left_join(left_join(
  left_join(interactions_by_age_group, py_on_treatment_by_age_group,
                             by = c("scenario", "sim")),
                             dalys_averted_by_age_group, by = c("scenario", "sim")),
  deaths_averted_by_age_group, by = c("scenario", "sim"))
df_by_age_group$scenario <- factor(df_by_age_group$scenario)
levels(df_by_age_group$scenario) <- list("15-30 years" = "screen_2020_monit_sim6",
                                        "30-45 years" = "screen_2020_monit_sim8",
                                        "45+ years" = "screen_2020_monit_sim10")



# Incremental DALYs averted per monitoring assessment
p1 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/monitoring_assessments)) +
  ylab("Incremental DALYs averted\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,1)

# Incremental DALYs averted per treatment initiations
p2 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/treatment_initiations)) +
  ylab("Incremental DALYs averted\nper treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,30)

# Incremental DALYs averted per PY on treatment
p3 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/py_on_treatment)) +
  ylab("Incremental DALYs averted\nper PY on treatment") +
  xlab("Monitored age group") +
  ylim(0,0.75)

# Incremental deaths averted per monitoring assessment
p4 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/monitoring_assessments)) +
  ylab("Incremental HBV-related deaths averted\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,0.04)

# Incremental deaths  averted per treatment initiations
p5 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/treatment_initiations)) +
  ylab("Incremental HBV-related deaths averted\nper treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,0.7)

# Incremental deaths  averted per PY on treatment
p6 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/py_on_treatment)) +
  ylab("Incremental HBV-related deaths averted\nper PY on treatment") +
  xlab("Monitored age group") +
  ylim(0,0.05)

# Incremental treatment initiations per monitoring assessment
p7 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations/monitoring_assessments)) +
  ylab("New treatment initiations\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,0.2)

# Incremental PY on treatment per monitoring assessment
p8 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/monitoring_assessments)) +
  ylab("Incremental person-years on treatment\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,1.75)

# Incremental PY on treatment per treatment initiation
p9 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/treatment_initiations)) +
  ylab("Incremental person-years on treatment\nper treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,35)

grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3)
grid.arrange(p7,p8, ncol = 1)

# Compare this with previous result!! Totally different because of change in
# initial treatment initiations (out3_it) - this would be case in all scenarios including
#

# Additional plots per carrier (especially monitoring assessments per carrier)
# would be helpful (this would be basically the cohort size over time in the different age groups)
# But monitoring assessments is kind of a proxy for carriers over time (except
# that age groups are not equally broad)
# Carriers at entry would also be interesting.
# Also per remaining life expectancy/average age at death? But this would
# need to be calculated in group of given age only (not entire cohort)

# Lifetime risk of HBV death
# 15-30 year olds (monit_sim6), 30-45 year olds (monit_sim8) and 45+ year olds (monit_sim10)
# Extract incremental impact and interactions compared to no monitoring

quantile(monit_out6$cohort_cum_hbv_deaths[,-1]/monit_out6$cohort_size[,-1],
         c(0.5,0.025,0.975))
# Lifetime risk of HBV death in 15-30 year old screened+treated+monitored carriers =
# 7% (3-13%)
quantile(monit_out10$cohort_cum_hbv_deaths[,-1]/monit_out10$cohort_size[,-1],
         c(0.5,0.025,0.975))
# Lifetime risk of HBV death in 45+ year old screened+treated+monitored carriers =
# 5% (2-9%)


# But need cohort simulations of a cohort of this age without treatment/monitoring
# to look at reduction in lifetime risk! Have simulations of screening+treatment without
# monitoring.

# Average age at death in screened+treated cohort
ggplot(gather(rbind(out3_it$cohort_age_at_death, out6_it$cohort_age_at_death,
                    out1_it$cohort_age_at_death),
       key = "sim", value = "age", -scenario))+
  geom_boxplot(aes(x=reorder(scenario, age),y=age))

quantile(out6_it$cohort_age_at_death[,-1]-out3_it$cohort_age_at_death[,-1],
         c(0.5,0.025,0.975))*12
# Average age at death postponed by 8 (3-20) months by monitoring every year
quantile(out3_it$cohort_age_at_death[,-1]-out1_it$cohort_age_at_death[,-1],
         c(0.5,0.025,0.975))
# Average age at death postponed by 1.5 (0.7-3) years by treating without monitoring
quantile(((out3_it$cohort_cum_hbv_deaths[,-1]/out3_it$cohort_size[,-1])-
           (out6_it$cohort_cum_hbv_deaths[,-1]/out6_it$cohort_size[,-1]))/
           (out3_it$cohort_cum_hbv_deaths[,-1]/out3_it$cohort_size[,-1]),
         c(0.5,0.025,0.975))
# Lifetime risk of dying from HBV reduced by 4% (1-8%) by monitoring every year
# in absolute terms or 52% (29-74%) relative to lifetime risk with no monitoring!

quantile(((out1_it$cohort_cum_hbv_deaths[,-1]/out1_it$cohort_size[,-1])-
            (out3_it$cohort_cum_hbv_deaths[,-1]/out3_it$cohort_size[,-1]))/
           (out1_it$cohort_cum_hbv_deaths[,-1]/out1_it$cohort_size[,-1]),
         c(0.5,0.025,0.975))
# Lifetime risk of dying from HBV reduced by 6% (3-12) by treating without monitoring
# in absolute terms or 45% (31-59%) relative to lifetime risk with no treatment!


## EFFECT OF SCREENING WITHOUT MONITORING BY AGE

# DALYs averted
dalys_averted_by_screened_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a2_out3_it,
                                                  a4_out3_it,
                                                  a5_out3_it),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "No treatment")
dalys_averted_by_screened_age_group <- subset(dalys_averted_by_screened_age_group,
                                              type == "number_averted" &
                                       by_year==2100) %>%
  select(scenario, sim, value)
dalys_averted_by_screened_age_group$sim <- gsub("[^0-9]", "",
                                                dalys_averted_by_screened_age_group$sim)
colnames(dalys_averted_by_screened_age_group)[3] <- "dalys_averted"

deaths_averted_by_screened_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a2_out3_it,
                                                  a4_out3_it,
                                                  a5_out3_it),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "No treatment")
deaths_averted_by_screened_age_group <- subset(deaths_averted_by_screened_age_group,
                                               type == "number_averted" &
                                        by_year==2100) %>%
  select(scenario, sim, value)
deaths_averted_by_screened_age_group$sim <- gsub("[^0-9]", "", deaths_averted_by_screened_age_group$sim)
colnames(deaths_averted_by_screened_age_group)[3] <- "deaths_averted"

# In terms of interactions, want the HBsAg tests, initial clinical assessments,
# initial treatment initiations and person years of treatment due to initial treatment

interactions_by_screened_age_group <- rbind(
  cbind(scenario = "a2_screen_2020_monit_0",
        left_join(left_join(
          gather(a2_out3_it$interactions[[16]]$total_screened[-c(1:3)],
                 key = "sim", value = "hbsag_tests"),
          gather(a2_out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                 key = "sim", value = "clinical_assessments"), by = "sim"),
          gather(a2_out3_it$interactions[[16]]$total_treated[-c(1:3)],
                 key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "a4_screen_2020_monit_0",
        left_join(left_join(
          gather(a4_out3_it$interactions[[16]]$total_screened[-c(1:3)],
                 key = "sim", value = "hbsag_tests"),
          gather(a4_out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                 key = "sim", value = "clinical_assessments"), by = "sim"),
          gather(a4_out3_it$interactions[[16]]$total_treated[-c(1:3)],
                 key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "a5_screen_2020_monit_0",
        left_join(left_join(
          gather(a5_out3_it$interactions[[16]]$total_screened[-c(1:3)],
                 key = "sim", value = "hbsag_tests"),
          gather(a5_out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                 key = "sim", value = "clinical_assessments"), by = "sim"),
          gather(a5_out3_it$interactions[[16]]$total_treated[-c(1:3)],
                 key = "sim", value = "treatment_initiations"), by = "sim")))
interactions_by_screened_age_group$sim <- gsub("[^0-9]", "", interactions_by_screened_age_group$sim)

py_on_treatment_by_screened_age_group <-rbind(
  data.frame(scenario = "a2_screen_2020_monit_0",
             sim = names(a2_out3_it$py_on_treatment[[16]]),
             py_on_treatment = a2_out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "a4_screen_2020_monit_0",
             sim = names(a4_out3_it$py_on_treatment[[16]]),
             py_on_treatment = a4_out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "a5_screen_2020_monit_0",
             sim = names(a5_out3_it$py_on_treatment[[16]]),
             py_on_treatment = a5_out3_it$py_on_treatment[[16]]))

df_by_screened_age_group <- left_join(left_join(
  left_join(interactions_by_screened_age_group, py_on_treatment_by_screened_age_group,
            by = c("scenario", "sim")),
  dalys_averted_by_screened_age_group, by = c("scenario", "sim")),
  deaths_averted_by_screened_age_group, by = c("scenario", "sim"))
df_by_screened_age_group$scenario <- factor(df_by_screened_age_group$scenario)
levels(df_by_screened_age_group$scenario) <- list("15-30 years" = "a4_screen_2020_monit_0",
                                         "30-45 years" = "a5_screen_2020_monit_0",
                                         "45+ years" = "a2_screen_2020_monit_0")

# Add per sAg test

# Incremental DALYs averted per clinical assessment
ps1 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/clinical_assessments)) +
  ylab("Incremental DALYs averted\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,4.5)

# Incremental DALYs averted per treatment initiations
ps2 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/treatment_initiations)) +
  ylab("Incremental DALYs averted\nper treatment initiation") +
  xlab("Screened age group") +
  ylim(0,30)

# Incremental DALYs averted per PY on treatment
ps3 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/py_on_treatment)) +
  ylab("Incremental DALYs averted\nper PY on treatment") +
  xlab("Screened age group") +
  ylim(0,1)

# Incremental deaths averted per monitoring assessment
ps4 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/clinical_assessments)) +
  ylab("Incremental HBV-related deaths averted\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,0.2)

# Incremental deaths  averted per treatment initiations
ps5 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/treatment_initiations)) +
  ylab("Incremental HBV-related deaths averted\nper treatment initiation") +
  xlab("Screened age group") +
  ylim(0,0.8)

# Incremental deaths  averted per PY on treatment
ps6 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/py_on_treatment)) +
  ylab("Incremental HBV-related deaths averted\nper PY on treatment") +
  xlab("Screened age group") +
  ylim(0,0.04)

# Incremental treatment initiations per clinical assessment
ps7 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations/clinical_assessments)) +
  ylab("New treatment initiations\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,0.4)

# Incremental PY on treatment per clinical assessment
ps8 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/clinical_assessments)) +
  ylab("Incremental person-years on treatment\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,7.5)

# Incremental PY on treatment per treatment initiation
ps9 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/treatment_initiations)) +
  ylab("Incremental person-years on treatment\nper treatment initiation") +
  xlab("Screened age group") +
  ylim(0,45)

grid.arrange(ps1,ps2,ps3,ps4,ps5,ps6, ncol = 3)

grid.arrange(ps7,ps9,ps8,p7,p9, p8, ncol = 3)

# Might want to show something like: % of total treatment initiations
# after first vs monitoring assessments by age

## 5) Screening by age: cost-effectiveness of including <30 year olds in one-time screen ----

dalys_averted <- cbind(rep(0,183),
                       t((out2$dalys[[16]][,-c(1:3)]-a4_out3_it$dalys[[16]][,-c(1:3)])+
                           (out2$dalys[[16]][,-c(1:3)]-a5_out3_it$dalys[[16]][,-c(1:3)])),
                       t(out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)]))

cost <- cbind(rep(0,183),
              t((a4_out3_it$interactions[[16]]$total_screened[,-c(1:3)]+
                   a5_out3_it$interactions[[16]]$total_screened[,-c(1:3)])*8.3+
                  (a4_out3_it$interactions[[16]]$total_assessed[,-c(1:3)]+
                     a5_out3_it$interactions[[16]]$total_assessed[,-c(1:3)])*84.4+
                  (a4_out3_it$py_on_treatment[[16]]+a5_out3_it$py_on_treatment[[16]])*60),
              t((a4_out3_it$interactions[[16]]$total_screened[,-c(1:3)]+
                   a5_out3_it$interactions[[16]]$total_screened[,-c(1:3)]+
                   a2_out3_it$interactions[[16]]$total_screened[,-c(1:3)])*8.3+
                  (a4_out3_it$interactions[[16]]$total_assessed[,-c(1:3)]+
                     a5_out3_it$interactions[[16]]$total_assessed[,-c(1:3)]+
                     a2_out3_it$interactions[[16]]$total_assessed[,-c(1:3)])*84.4+
                  (a4_out3_it$py_on_treatment[[16]]+
                     a5_out3_it$py_on_treatment[[16]]+a2_out3_it$py_on_treatment[[16]])*60))

ceef.plot(bcea(e=dalys_averted,
               c=cost,
               ref=1,
               interventions=c("No treatment", "Screen 30+", "Screen all"),
               Kmax=100000000000000),
          graph="base", relative = FALSE)
# Extending the screen to all ages is still below WTP (on the mean)
