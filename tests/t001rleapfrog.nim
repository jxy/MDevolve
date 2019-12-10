import mdevolve
import sys2d, utils

runtest("Recursive Leapfrog", 2, @[-1.326510932742053e-05, -3.326948818882514e-06, -8.324050204677746e-07]):
  let
    (V, K) = newIntegratorPair((updateVFh, updateVFl0, updateVFl1), updateK)
    (VFh, VFl0, VFl1) = V
  let
    VK = mkLeapFrog(VFh, K, 2)
    Vl0VK = mkLeapFrog(VFl0, VK, 2)
  mkLeapFrog(VFl1, Vl0VK)

runtest("Recursive Leapfrog", 2, @[-1.326510932742053e-05, -3.326948818882514e-06, -8.324050204677746e-07]):
  let (VFhl01, K) = newIntegratorPair(updateVFhl01, updateK)
  let
    VFh = VFhl01[0]
    VFl0 = VFhl01[1]
    VFl1 = VFHl01[2]
  let
    VK = mkLeapFrog(VFh, K, 2)
    Vl0VK = mkLeapFrog(VFl0, VK, 2)
  mkLeapFrog(VFl1, Vl0VK)
