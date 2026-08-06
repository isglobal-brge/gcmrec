# Tests for the graphics layer (ggplot2) and the interpretation helpers.

fit_graphics <- local({
  data(readmission, envir = environment())
  suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
                          data = readmission, s = 3000))
})

test_that("plot.gcmrec returns a ggplot with a confidence band", {
  p <- plot(fit_graphics)
  expect_s3_class(p, "ggplot")

  # The band exists and contains the estimate
  d <- p$data
  expect_true(all(c("lower", "upper") %in% names(d)))
  expect_true(all(d$lower <= d$estimate + 1e-8))
  expect_true(all(d$upper >= d$estimate - 1e-8))
  # Survival stays within [0, 1]
  expect_true(all(d$lower >= 0 & d$upper <= 1))

  expect_s3_class(plot(fit_graphics, type.plot = "hazard"), "ggplot")
  expect_s3_class(plot(fit_graphics, conf.int = FALSE), "ggplot")
  expect_error(plot(fit_graphics, type.plot = "nope"), "hazard or survival")
})

test_that("the band widens as the confidence level increases", {
  d95 <- plot(fit_graphics)$data
  d99 <- plot(fit_graphics, level = 0.99)$data
  width95 <- d95$upper - d95$lower
  width99 <- d99$upper - d99$lower
  expect_true(mean(width99) > mean(width95))
})

test_that("plotBaseline draws several models together", {
  data(readmission)
  mod.min <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
           data = readmission, s = 3000, typeEffage = "minimal")
  )
  p <- plotBaseline(list(perfect = fit_graphics, minimal = mod.min))
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$model), c("perfect", "minimal"))
  expect_error(plotBaseline(list(1, 2)), "gcmrec objects")
})

test_that("plotForest shows the HRs with their intervals", {
  p <- plotForest(fit_graphics)
  expect_s3_class(p, "ggplot")
  # One row per covariate (alpha is left out)
  expect_equal(nrow(p$data), length(coef(fit_graphics)) - 1)
  expect_equal(p$data$estimate, unname(exp(coef(fit_graphics)[-1])),
               tolerance = 1e-8)

  # Custom labels are applied
  p2 <- plotForest(fit_graphics, labels = c(sex = "Female vs male"))
  expect_true("Female vs male" %in% levels(p2$data$term))
})

test_that("graph.caltimes returns a ggplot and sorts the subjects", {
  data(readmission)
  sub <- readmission[readmission$id %in% 1:30, ]

  expect_s3_class(graph.caltimes(sub), "ggplot")
  expect_s3_class(graph.caltimes(sub, sortevents = "events"), "ggplot")
  expect_s3_class(graph.caltimes(sub, sortevents = "none"), "ggplot")
  expect_s3_class(graph.caltimes(sub, var = sub$sex), "ggplot")

  # sortevents = "followup" actually sorts (it used to be ignored)
  p <- graph.caltimes(sub, sortevents = "followup")
  seg <- p$layers[[1]]$data
  expect_false(is.unsorted(seg$tau[order(seg$row)]))

  # Legacy list format
  data(hydraulic)
  expect_s3_class(graph.caltimes(hydraulic), "ggplot")
})

test_that("mcf estimates the mean cumulative function", {
  data(readmission)
  m <- mcf(readmission)
  expect_s3_class(m, "mcf")

  # It is monotone non-decreasing and ends at the total events per subject
  expect_false(is.unsorted(m$estimate))
  expect_true(all(m$lower <= m$estimate & m$upper >= m$estimate))
  expect_true(all(m$n.risk > 0))
  expect_equal(sum(m$n.event), sum(readmission$event == 1))

  expect_s3_class(plot(m), "ggplot")
  expect_output(print(m), "Mean cumulative function")
})

test_that("mcf by group separates the cumulative risk", {
  data(readmission)
  m <- mcf(readmission, group = readmission$dukes)
  expect_true("group" %in% names(m))
  expect_equal(nlevels(m$group), 3L)

  # Stage D accumulates more events than stage A-B
  final <- tapply(m$estimate, m$group, max)
  expect_gt(final[["3"]], final[["1"]])

  expect_s3_class(plot(m), "ggplot")
  expect_error(mcf(readmission, group = 1:3), "one value per row")
})

test_that("predict computes risks and per-profile curves", {
  profiles <- data.frame(dukes = c(1, 3), sex = c(1, 1))

  risk <- predict(fit_graphics, profiles, type = "risk")
  lp <- predict(fit_graphics, profiles, type = "lp")
  expect_equal(risk, exp(lp))
  # Stage D carries more risk than stage A-B
  expect_gt(risk[[2]], risk[[1]])

  surv <- predict(fit_graphics, profiles, type = "survival")
  expect_equal(nlevels(surv$profile), 2L)
  expect_true(all(surv$estimate >= 0 & surv$estimate <= 1))

  # The curve for the higher-risk profile lies below
  s1 <- surv$estimate[surv$profile == levels(surv$profile)[1]]
  s2 <- surv$estimate[surv$profile == levels(surv$profile)[2]]
  expect_true(all(s2 <= s1 + 1e-8))

  haz <- predict(fit_graphics, profiles, type = "hazard")
  expect_false(is.unsorted(haz$estimate[haz$profile ==
                                          levels(haz$profile)[1]]))

  expect_s3_class(plotPredict(fit_graphics, profiles), "ggplot")
  expect_error(predict(fit_graphics), "newdata")
  expect_error(predict(fit_graphics, data.frame(nope = 1)), "does not")
})

test_that("anova tests alpha = 1 via a likelihood ratio test", {
  a <- anova(fit_graphics)
  expect_s3_class(a, "anova.gcmrec")
  expect_equal(nrow(a), 2L)
  expect_equal(a$Df[2], 1)
  # The model with free alpha cannot have a lower likelihood
  expect_gte(a$logLik[2], a$logLik[1])
  expect_true(a$Chisq[2] >= 0)
  expect_true(a$`Pr(>Chisq)`[2] >= 0 && a$`Pr(>Chisq)`[2] <= 1)
  expect_output(print(a), "Likelihood ratio test")
})

test_that("anova compares models nested in the covariates", {
  data(readmission)
  small <- suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                                   data = readmission, s = 3000))
  a <- anova(small, fit_graphics)
  expect_equal(nrow(a), 2L)
  expect_equal(a$Df[2], 1)
  expect_gte(a$Chisq[2], 0)
})

test_that("anova refuses to compare models with and without frailties", {
  skip_on_cran()
  data(readmission)
  # The likelihoods are not comparable: the one from the frailty model
  # is conditional on the estimated frailties
  mod.fra <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
           data = readmission, s = 3000, Frailty = TRUE)
  )
  expect_error(anova(fit_graphics, mod.fra), "same scale")
})

test_that("anova(fit) refits alpha = 1 from any environment", {
  data(readmission)
  # Fit done inside another function: the refit cannot rely on
  # re-evaluating the original call
  ajusta <- function(...) {
    suppressWarnings(gcmrec(..., data = readmission, s = 3000))
  }
  m <- ajusta(Survr(id, time, event) ~ as.factor(dukes))
  a <- anova(m)
  expect_s3_class(a, "anova.gcmrec")

  # And it matches fitting the reduced model by hand
  m0 <- suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                                data = readmission, s = 3000,
                                rhoFunc = "Identity"))
  expect_equal(a$logLik[1], m0$loglik, tolerance = 1e-8)
})

test_that("keep.data = FALSE slims the object and anova explains it", {
  data(readmission)
  light <- suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                                   data = readmission, s = 3000,
                                   keep.data = FALSE))
  expect_null(light$gcmrec.data)
  expect_lt(as.numeric(object.size(light)),
            as.numeric(object.size(fit_graphics)))
  expect_error(anova(light), "keep.data = FALSE")

  # Comparing two models still works without the stored data
  light2 <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
           data = readmission, s = 3000, keep.data = FALSE)
  )
  expect_s3_class(anova(light, light2), "anova.gcmrec")
})

test_that("anova rejects comparisons that are not a valid LRT", {
  data(readmission)
  small <- suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                                   data = readmission, s = 3000))

  # Different effective age scale: the models are not nested
  mod.min <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ as.factor(dukes),
           data = readmission, s = 3000, typeEffage = "minimal")
  )
  expect_error(anova(small, mod.min), "effective age scale")

  # Different data
  sub <- readmission[readmission$id %in% unique(readmission$id)[1:150], ]
  mod.sub <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ as.factor(dukes), data = sub, s = 3000)
  )
  expect_error(anova(mod.sub, fit_graphics), "same data")

  # Non-nested covariates
  mod.other <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ sex + chemo, data = readmission,
           s = 3000)
  )
  expect_error(anova(mod.other, fit_graphics), "not nested")

  # Decreasing order of complexity
  expect_error(anova(fit_graphics, small), "not nested")
})

test_that("anova warns when comparing two frailty models", {
  skip_on_cran()
  data(readmission)
  f1 <- suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                                data = readmission, s = 3000,
                                Frailty = TRUE))
  f2 <- suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes) +
                                  sex, data = readmission, s = 3000,
                                Frailty = TRUE))
  expect_warning(anova(f1, f2), "only.*approximate")
})

test_that("theme_gcmrec can be composed with the plots", {
  p <- plot(fit_graphics) + theme_gcmrec(base_size = 14)
  expect_s3_class(p, "ggplot")
})
