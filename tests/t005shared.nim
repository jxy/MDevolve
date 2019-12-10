import mdevolve
import sys2d, utils

let
  (VFhl01, K) = newIntegratorPair(updateVFhl01, updateK)
  VFh = VFhl01[0]
  VFl0 = VFhl01[1]
  VFl1 = VFHl01[2]

runtest("3*LeapFrog", 2, @[-1.326510932742053e-05, -3.326948818882514e-06, -8.324050204677746e-07]):
  # Exactly as in t001rleapfrog.nim
  var
    VK = mkLeapFrog(VFh, K, 4)
    Vl0K = mkLeapFrog(VFl0, K, 2)
    Vl1K = mkLeapFrog(VFl1, K, 1)
  newParallelEvolution(VK, Vl0K, Vl1K)

runtest("3*Omelyan2MN", 2, @[-5.924369554932696e-07, -1.493390566764674e-07, -3.741168219661972e-08]):
  # Same amount of integration steps in t002recursive.nim.
  # Since the steps are not equal size in recursive integrators, this is different.
  var
    VK = mkOmelyan2MN(K, VFh, 36)
    Vl0K = mkOmelyan2MN(K, VFl0, 6)
    Vl1K = mkOmelyan2MN(K, VFl1, 1)
  newParallelEvolution(VK, Vl0K, Vl1K)

runtest("3*Omelyan4MN5FV", 1, @[-6.775154013372031e-07, -1.698951566098117e-07, -4.25064610176662e-08]):
  # Same amount of integration steps in t002recursive.nim.
  # Since the steps are not equal size in recursive integrators, this is different.
  # Because of how we arrange the parallel evolution, this is no longer order 4.
  var
    VK = mkOmelyan4MN5FV(K, VFh, 100)
    Vl0K = mkOmelyan4MN5FV(K, VFl0, 10)
    Vl1K = mkOmelyan4MN5FV(K, VFl1, 1)
  newParallelEvolution(VK, Vl0K, Vl1K)
