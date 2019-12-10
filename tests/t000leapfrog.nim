import mdevolve
import sys2d, utils

runtest("Leapfrog", 2, @[7.544205908427415e-06, 1.863113559208429e-06, 4.643416811056511e-07]):
  let (V, K) = newIntegratorPair(updateVF, updateK)
  mkLeapFrog(V, K)
