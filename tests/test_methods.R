# In-memory mathematical checks only. No study records or derived data are used.
set.seed(42)
w <- c(.35, .65); p0 <- c(.8, .2); p1 <- c(.6, .4)
r0 <- c(.1, .3); r1 <- c(.15, .25)
total <- sum(w * (p1 * r1 - p0 * r0))
composition <- sum(w * (p1 - p0) * (r0 + r1) / 2)
within <- sum(w * (p1 + p0) * (r1 - r0) / 2)
stopifnot(abs(total - composition - within) < 1e-14)
# A zero Poisson plug-in count is fixed; a Jeffreys posterior is not degenerate.
stopifnot(all(rpois(1000, 0) == 0))
posterior <- rgamma(50000, shape = .5, rate = 1000)
stopifnot(all(posterior > 0), var(posterior) > 0,
          abs(mean(posterior) - .0005) < .00003)
# Equal-weight Fay-Feuer limits reduce to exact scaled Poisson limits.
count <- 3; a <- .002; rate <- count * a; v <- count * a^2
lo <- qgamma(.025, rate^2 / v, scale = v / rate)
hi <- qgamma(.975, (rate + a)^2 / (v + a^2), scale = (v + a^2) / (rate + a))
stopifnot(abs(lo - a * qchisq(.025, 2 * count) / 2) < 1e-14,
          abs(hi - a * qchisq(.975, 2 * (count + 1)) / 2) < 1e-14)
# With no events the conservative gamma upper limit stays positive.
stopifnot(qgamma(.975, shape = 1, scale = a) > 0)
cat('PASS: decomposition identity, zero-count uncertainty and gamma limits.\n')
