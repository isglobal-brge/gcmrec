# Unit tests for addCenTime: it adds a censored record with time 0
# to subjects whose last record is an event.

test_that("addCenTime adds a censored record only where one is missing", {
  datin <- data.frame(
    id    = c(1, 1, 2),
    time  = c(5, 3, 7),
    event = c(1, 1, 0)
  )
  out <- addCenTime(datin)

  # Subject 1 gains a censored record with time 0; subject 2 is unchanged
  expect_equal(nrow(out), 4L)
  added <- out[out$id == 1, ][3, ]
  expect_equal(added$time, 0)
  expect_equal(added$event, 0)
  expect_equal(nrow(out[out$id == 2, ]), 1L)
})

test_that("addCenTime does not modify data that are already complete", {
  datin <- data.frame(
    id    = c(1, 1, 2),
    time  = c(5, 3, 7),
    event = c(1, 0, 0)
  )
  expect_equal(addCenTime(datin), datin)
})

test_that("addCenTime honors non-standard column positions", {
  datin <- data.frame(
    evento = c(1, 0),
    sujeto = c(1, 1),
    tiempo = c(5, 3)
  )
  out <- addCenTime(datin, id = 2, time = 3, event = 1)
  expect_equal(out, datin)
})
