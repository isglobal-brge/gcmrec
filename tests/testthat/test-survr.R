# Unit tests for the Survr response object.

test_that("Survr builds a valid object from correct data", {
  x <- Survr(id = c(1, 1, 2), time = c(5, 3, 7), event = c(1, 0, 0))
  expect_true(is.Survr(x))
  expect_equal(dim(x), c(3L, 3L))
  expect_equal(colnames(x), c("id", "time", "event"))
})

test_that("Survr requires one censored time per subject", {
  # Subject 1 has no censored record (event == 0)
  expect_error(
    Survr(id = c(1, 1, 2), time = c(5, 3, 7), event = c(1, 1, 0)),
    "censored"
  )
})

test_that("Survr requires events coded 0-1", {
  expect_error(
    Survr(id = c(1, 1, 2), time = c(5, 3, 7), event = c(2, 0, 0)),
    "0-1"
  )
})

test_that("is.Survr tells Survr objects apart from plain matrices", {
  expect_false(is.Survr(matrix(1:6, ncol = 3)))
})
