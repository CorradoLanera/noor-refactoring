# Packages --------------------------------------------------------
library(testthat)
library(bench)
library(profvis)
library(readr)
library(here)
library(geepack)
library(janitor)




# Setup -----------------------------------------------------------

dat <- here("data/reduced 1.rds") |>
  read_rds() |>
  clean_names()

clusz <- table(dat[["patient"]])
waves <- dat[["time"]]
corstrv <- "unstructured"




# Functions -------------------------------------------------------
compose_waves <- function(x, wave) {
  n <- length(wave)

  stringi::stri_c(
    wave[[x]],
    wave[seq(from = x + 1L, to = n)],
    sep = ":"
  )
}

crossutri <- function(wave) {
  n <- length(wave)

  if (n == 1) {
    return(NULL)
  }

  wave[] <- paste(wave)
  seq_len(n - 1L) |>
    lapply(compose_waves, wave)
}


generate_design <- function(clusz, waves, corstrv) {
  if (corstrv == 1) {
    return(matrix(0, 0, 0))
  }
  crs <- clusz * (clusz - 1) / 2
  if (corstrv == 2 || corstrv == 3) {
    ans <- matrix(1, length(clusz), 1)
    colnames(ans) <- c("alpha")
  } else {
    id <- rep(seq_along(clusz), clusz)
    z1 <- unlist(lapply(split(waves, id), crossutri), use.names = FALSE)
    z2 <- unlist(crossutri(seq_len(max(clusz))), use.names = FALSE)
    z <- factor(z1, levels = unique.default(z2))

    # Avoid dispatch for model.matrix
    ans <- model.matrix.default(~ z - 1)

    # stringi::stri_c much faster then paste
    znames <- stringi::stri_c("alpha", z2, sep = ".")

    # faster than
    # > colnames(ans) <- dn (much less gc()s) # Original, and slowest!
    # and faster then
    # > attributes(ans)[["dimnames"]][[2]] <- znames
    dn <- dimnames(ans)
    dn[[2L]] <- znames
    dimnames(ans) <- dn
  }
  ans
}





# Test ------------------------------------------------------------

with_reporter(default_reporter(), {

  test_that("refactoring works", {
    # setup
    expected <- genZcor(clusz, waves, corstrv)

    # evaluation
    res <- generate_design(clusz, waves, corstrv)

    # tests
    res |>
      expect_equal(expected)
  })

})


# Benchmark -------------------------------------------------------

bench::mark(
  refactored = generate_design(clusz, waves, corstrv),
  original = genZcor(clusz, waves, corstrv),
  min_iterations = 2,
  check = FALSE,
  memory = FALSE,
  filter_gc = FALSE
) |>
  plot()



# Profile ---------------------------------------------------------

prof <- profvis({
  refactored <- generate_design(clusz, waves, corstrv)
  original <- genZcor(clusz, waves, corstrv)
})
prof
