import mdevolve/v0
import sys2d, utils

runtest("Leapfrog", 2):
  let
    VK = mkLeapFrog(updateVFh, updateK, 2)
    Vl0VK = mkLeapFrog(updateVFl0, VK.evolve, 2)
  mkLeapFrog(updateVFl1, Vl0VK.evolve)
runtest("SW92", 2):
  let
    VK = mkSW92(updateK, updateVFh, 2)
    Vl0VK = mkSW92(VK.evolve, updateVFl0, 2)
  mkSW92(Vl0VK.evolve, updateVFl1)
runtest("Omelyan2MN", 2):
  let
    VK = mkOmelyan2MN(updateK, updateVFh, 2)
    Vl0VK = mkOmelyan2MN(VK.evolve, updateVFl0, 2)
  mkOmelyan2MN(Vl0VK.evolve, updateVFl1)
runtest("Omelyan4MN4FP", 1):
  let
    VK = mkOmelyan4MN4FP(updateK, updateVFh, 2)
    Vl0VK = mkOmelyan4MN4FP(VK.evolve, updateVFl0, 2)
  mkOmelyan4MN4FP(Vl0VK.evolve, updateVFl1)
runtest("Omelyan4MN5FV", 1):
  let
    VK = mkOmelyan4MN5FV(updateVFh, updateK, 2)
    Vl0VK = mkOmelyan4MN5FV(updateVFl0, VK.evolve, 2)
  mkOmelyan4MN5FV(updateVFl1, Vl0VK.evolve)
runtest("LF/4MN4FP/4MN5FV", 2):
  let
    VK = mkLeapFrog(updateVFh, updateK, 64)
    Vl0VK = mkOmelyan4MN4FP(VK.evolve, updateVFl0, 1)
  mkOmelyan4MN5FV(updateVFl1, Vl0VK.evolve)
runtest("FGYin11", 1):
  let
    VK = mkFGYin11(updateVFh, updateK, updateVfgFh, (xsave.save s.x), (s.x.load xsave), 2)
    Vl0VK = mkFGYin11(updateVFl0, VK.evolve, updateVfgFl0, (xsave.save s.x), (s.x.load xsave), 2)
  mkFGYin11(updateVFl1, Vl0VK.evolve, updateVfgFl1, (xsave.save s.x), (s.x.load xsave))
runtest("3*LeapFrog", 2):
  var
    VK = mkLeapFrog(updateVFh, updateK, 3, shared = 1)
    Vl0K = mkLeapFrog(updateVFl0, updateK, 2, shared = 1)
    Vl1K = mkLeapFrog(updateVFl1, updateK, 1, shared = 1)
  mkSharedEvolution(VK, Vl0K, Vl1K)
runtest("3*Omelyan2MN", 2):
  var
    VK = mkOmelyan2MN(updateK, updateVFh, 3, shared = 0)
    Vl0K = mkOmelyan2MN(updateK, updateVFl0, 1, shared = 0)
    Vl1K = mkOmelyan2MN(updateK, updateVFl1, 1, shared = 0)
  mkSharedEvolution(VK, Vl0K, Vl1K)
runtest("3*Omelyan4MN4FP", 2):
  var
    VK = mkOmelyan4MN4FP(updateK, updateVFh, 3, shared = 0)
    Vl0K = mkOmelyan4MN4FP(updateK, updateVFl0, 1, shared = 0)
    Vl1K = mkOmelyan4MN4FP(updateK, updateVFl1, 1, shared = 0)
  mkSharedEvolution(VK, Vl0K, Vl1K)
runtest("3*Omelyan4MN5FV", 2):
  var
    VK = mkOmelyan4MN5FV(updateVFh, updateK, 3, shared = 1)
    Vl0K = mkOmelyan4MN5FV(updateVFl0, updateK, 1, shared = 1)
    Vl1K = mkOmelyan4MN5FV(updateVFl1, updateK, 1, shared = 1)
  mkSharedEvolution(VK, Vl0K, Vl1K)
runtest("3*FGYin11", 2):
  var
    VK = mkFGYin11(updateVFh, updateK, updateVfgFh, (xsave.save s.x), (s.x.load xsave), 36, shared = 1)
    Vl0K = mkFGYin11(updateVFl0, updateK, updateVfgFl0, (xsave.save s.x), (s.x.load xsave), 6, shared = 1)
    Vl1K = mkFGYin11(updateVFl1, updateK, updateVfgFl1, (xsave.save s.x), (s.x.load xsave), 1, shared = 1)
  mkSharedEvolution(VK, Vl0K, Vl1K)
runtest("2MN+2*4MN4FP", 2):
  var
    VK = mkOmelyan2MN(updateK, updateVFh, 8, shared = 0)
    Vl0K = mkOmelyan4MN4FP(updateK, updateVFl0, 1, shared = 0)
    Vl1K = mkOmelyan4MN4FP(updateK, updateVFl1, 1, shared = 0)
  mkSharedEvolution(VK, Vl0K, Vl1K)
