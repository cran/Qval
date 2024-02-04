################################################################
#                           Example 1                          #
#             The GDI method to validate Q-matrix              #
################################################################
set.seed(123)

library(Qval)

## generate Q-matrix and data
K <- 4
I <- 20
example.Q <- sim.Q(K, I)
IQ <- list(
  P0 = runif(I, 0.0, 0.2),
  P1 = runif(I, 0.8, 1.0)
)
example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")

## simulate random mis-specifications
example.MQ <- sim.MQ(example.Q, 0.1)

## using MMLE/EM to fit CDM model first
example.CDM.obj <- CDM(example.data$dat, example.MQ)

## using the fitted CDM.obj to avoid extra parameter estimation.
Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "GDI")


## also can validate the Q-matrix directly
Q.GDI.obj <- validation(example.data$dat, example.MQ)

## item level iteration
Q.GDI.obj <- validation(example.data$dat, example.MQ, method = "GDI", iter.level = "item", maxitr = 150)

## search method
Q.GDI.obj <- validation(example.data$dat, example.MQ, method = "GDI", search.method = "ESA")

## cut-off point
Q.GDI.obj <- validation(example.data$dat, example.MQ, method = "GDI", eps = 0.90)

## check QRR
print(getQRR(example.Q, Q.GDI.obj$Q.sug))
