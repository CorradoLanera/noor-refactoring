refactored_crossutri <- function(wave) {
  n <- length(wave)

  if (n == 1) {
    return(NULL)
  }

  ans <- rep(0, n * (n - 1) / 2)

  k <- 1
  for (i in seq_len(n - 1)) {
    for (j in seq(from = i + 1, to = n)) {
      ans[k] <- paste(wave[i], wave[j], sep = ":")
      k <- k + 1
    }
  }
  ans
}

refactored_genZcor <- function(clusz, waves, corstrv) {
  if (corstrv == 1) {
    return(matrix(0, 0, 0))
  }
  crs <- clusz * (clusz - 1) / 2
  if (corstrv == 2 || corstrv == 3) {
    ans <- matrix(1, length(clusz), 1)
    colnames(ans) <- c("alpha")
  } else {
    id <- rep(seq_along(clusz), clusz)
    z1 <- unlist(lapply(split(waves, id), refactored_crossutri))
    z2 <- unlist(refactored_crossutri(seq_len(max(clusz))))
    z <- factor(z1, levels = unique.default(z2))
    ans <- model.matrix(~ z - 1)
    znames <- paste("alpha", z2, sep = ".")
    colnames(ans) <- znames
  }
  ans
}
