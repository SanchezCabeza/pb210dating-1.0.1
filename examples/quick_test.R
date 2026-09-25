# Quick test for pb210dating version 1.0.1

library(pb210dating)

result <- Pb210Dating(
  LayerDating = 17,
  Path = tempdir(),
  MonteCarloRunsCF = 100,
  MonteCarloRunsCFCS = 100
)
