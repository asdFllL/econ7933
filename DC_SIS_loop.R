
# implement DC-SIS proposed in [Li, Zhong and Zhu (2012, JASA)]

# calculating distance correlation: dcorr()
dcov <- function(u, v, U2ind, U3ind) { # nolint

  S1 <- mean(abs(u[U2ind[, 1]] - u[U2ind[, 2]]) * abs(v[U2ind[, 1]] - v[U2ind[, 2]])) # nolint
  S2 <- mean(abs(u[U2ind[, 1]] - u[U2ind[, 2]])) * mean(abs(v[U2ind[, 1]] - v[U2ind[, 2]])) # nolint
  S3 <- mean(abs(u[U3ind[, 1]] - u[U3ind[, 3]]) * abs(v[U3ind[, 2]] - v[U3ind[, 3]])) # nolint

  sqrt(S1 + S2 - 2 * S3) # nolint
}

dcorr <- function(u, v, U2ind, U3ind) { # nolint

  dcov(u, v, U2ind, U3ind) / sqrt(dcov(u, u, U2ind, U3ind) * dcov(v, v, U2ind, U3ind)) # nolint
}

# select active predictors with given cutoff/threshold
DC_SIS_cutoff <- function(dataset, U2ind, U3ind, cutoff) { # nolint

  data <- dataset
  V <- data[, 1]
  X <- data[, (2:ncol(data))] # nolint
  omega <- c()
  p <- ncol(X) # nolint

  for (k in 1:p) {
    omega[k] <- dcorr(X[, k], V, U2ind, U3ind)^2 # nolint
  }

  # which((omega >= threshold) == TRUE)                # screening using threshold # nolint
  head(order(omega, decreasing = TRUE), n = cutoff) # screening using cutoff # nolint
  # omega
}

DC_SIS_threshold <- function(dataset, U2ind, U3ind, threshold) { # nolint

  data <- dataset
  V <- data[, 1]
  X <- data[, (2:ncol(data))] # nolint
  omega <- c()
  p <- ncol(X) # nolint

  for (k in 1:p) {
    omega[k] <- dcorr(X[, k], V, U2ind, U3ind)^2 # nolint
  }

  which((omega >= threshold) == TRUE) # screening using threshold # nolint
  # head(order(omega, decreasing = TRUE), n = cutoff)    # screening using cutoff # nolint
  # omega
}

DC_SIS_check <- function(dataset, U2ind, U3ind, cutoff) { # nolint

  data <- dataset
  V <- data[, 1]
  X <- data[, (2:ncol(data))] # nolint
  omega <- c()
  p <- ncol(X) # nolint

  for (k in 1:p) {
    omega[k] <- dcorr(X[, k], V, U2ind, U3ind)^2 # nolint
  }

  # which((omega >= threshold) == TRUE)                # screening using threshold # nolint
  # head(order(omega, decreasing = TRUE), n = cutoff)    # screening using cutoff # nolint
  omega
}


# 要用threshold再算一次，以防包括了特别小的被涵盖 R Function只输出最后一行

# U2ind/U3ind(2阶，3阶 U indecient) gtools permutations()
# logic probit 对照组设置


# U2ind = permutations(n, 2, repeats.allowed = TRUE)
# U3ind = permutations(n, 3, repeats.allowed = TRUE)

# first try
library("gtools")
library("haven")

U2 <- permutations(500, 2, repeats.allowed = TRUE)
U3 <- permutations(500, 3, repeats.allowed = TRUE)

data1 <- read_dta("hPHI_data_wave17.dta")

# 将PHI17变成第一列
cols <- colnames(data1)
new_cols <- c(cols[3], cols[1:2], cols[4:47])
data1 <- data1[, new_cols]

# 转化为numeric类型
data1 <- as.data.frame(lapply(data1, as.numeric))

View(data1)

data_test <- data1[1:500, ]

col_names <- colnames(data_test)
col_names[c(29, 30, 31, 46)]

t1 <- proc.time()

# DC_SIS_cutoff(data_test, U2, U3, 5)
col_names[DC_SIS_threshold(data_test, U2, U3, 0.1)]
# DC_SIS_check(data_test, U2, U3)



t2 <- proc.time()
t <- t2 - t1
print(paste0("执行时间：", t[3][1], "秒"))

loop_data <- function(dataset, n, U2, U3, threshold) {
  k <- nrow(dataset)
  m <- floor(k / n)
  output_list <- c()
  for (i in 1:m) {
    j <- (i - 1) * n + 1
    data_run <- dataset[j:i * n, ]
    result_1 <- DC_SIS_threshold(data_run, U2, U3, threshold)
    output_list <- append(output_list, result_1)
  }
  data_run_2 <- dataset[m * n + 1:k, ]
  result_2 <- DC_SIS_threshold(data_run_2, U2, U3, threshold)
  output_list <- append(output_list, result_2)
  output_list <- output_list[!duplicated(output_list)]
  output_list
}


data2 <- read_dta("PHI_wave17.dta")
View(data2)


test <- c(2, 2, 2, 4, 5, 6, 7, 7)
duplicated(test)
test <- test[!duplicated(test)]
test
