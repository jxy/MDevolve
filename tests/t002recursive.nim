import mdevolve
import sys2d, utils

template runs(setup:untyped):auto {.dirty.} =
  block:
    setup
    runtest("Recursive SW92", 2, @[-2.496588332734362e-06, -6.263639766856954e-07, -1.567297263083134e-07]):
      let
        VK = mkSW92(K, VFh, 2)
        Vl0VK = mkSW92(VK, VFl0, 2)
      mkSW92(Vl0VK, VFl1)
    runtest("Recursive Omelyan2MN", 2, @[-5.198453085775157e-07, -1.312348767434912e-07, -3.288839067749905e-08]):
      let
        VK = mkOmelyan2MN(K, VFh, 2)
        Vl0VK = mkOmelyan2MN(VK, VFl0, 2)
      mkOmelyan2MN(Vl0VK, VFl1)
    runtest("Recursive Omelyan4MN4FP", 1, @[8.033465626056113e-08, 5.035945438436329e-09, 3.149835947624524e-10]):
      let
        VK = mkOmelyan4MN4FP(K, VFh, 2)
        Vl0VK = mkOmelyan4MN4FP(VK, VFl0, 2)
      mkOmelyan4MN4FP(Vl0VK, VFl1)
    runtest("Recursive Omelyan4MN5FV", 1, @[2.510426089230577e-09, 1.586677456089092e-10, 9.971357073368381e-12]):
      let
        VK = mkOmelyan4MN5FV(K, VFh, 2)
        Vl0VK = mkOmelyan4MN5FV(VK, VFl0, 2)
      mkOmelyan4MN5FV(Vl0VK, VFl1)
    runtest("Recursive Omelyan4MN5FP", 1, @[1.231617319241707e-08, 7.717484429292654e-10, 4.824718402574035e-11]):
      let
        VK = mkOmelyan4MN5FP(K, VFh, 2)
        Vl0VK = mkOmelyan4MN5FP(VK, VFl0, 2)
      mkOmelyan4MN5FP(Vl0VK, VFl1)
    runtest("Recursive LF/4MN4FP/4MN5FV", 2, @[5.050797557970554e-10, 3.280709037767338e-11, 2.283506717049022e-12]):
      let
        VK = mkLeapFrog(K, VFh, 64)
        Vl0VK = mkOmelyan4MN4FP(VK, VFl0, 1)
      mkOmelyan4MN5FV(Vl0VK, VFl1)

runs:
  let
    (V, K) = newIntegratorPair((updateVFh, updateVFl0, updateVFl1), updateK)
    (VFh, VFl0, VFl1) = V

runs:
  let (VFhl01, K) = newIntegratorPair(updateVFhl01, updateK)
  let
    VFh = VFhl01[0]
    VFl0 = VFhl01[1]
    VFl1 = VFHl01[2]
