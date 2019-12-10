import mdevolve/v0
import sys2d, utils

runtest("3*LeapFrog", 2):
  var
    VK = mkLeapFrog(updateVFh, updateK, 3, shared = 1)
    Vl0K = mkLeapFrog(updateVFl0, updateK, 2, shared = 1)
    Vl1K = mkLeapFrog(updateVFl1, updateK, 1, shared = 1)
  mkSharedEvolution(VK, Vl0K, Vl1K)
