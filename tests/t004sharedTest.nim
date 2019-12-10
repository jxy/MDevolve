import mdevolve
import sys2d, utils

proc mkTest0(T,V:Integrator, steps = 1):Integrator =
  newSerialEvolution("Test0",
    @[(T, 0.375), (V, 0.4), (T, 0.125), (V, 0.2), (T, 0.125), (V, 0.4), (T, 0.375)],
    steps, 2)

proc mkTest1(T,V:Integrator, steps = 1):Integrator =
  newSerialEvolution("Test1",
    @[(T, 0.75), (V, 0.5), (T, -0.5), (V, 0.5), (T, 0.75)],
    steps, 2)

proc mkTest2(T,V:Integrator, steps = 1):Integrator =
  newSerialEvolution("Test2",
    @[(T, 0.5), (V, 0.25), (T, -0.25), (V, 0.25), (T, 0.5), (V, 0.25), (T, -0.25), (V, 0.25), (T, 0.5)],
    steps, 2)

proc mkTest3(T,V:Integrator, steps = 1):Integrator =
  newSerialEvolution("Test3",
    @[(T, 0.625), (V, 0.2), (T, -0.5), (V, 0.1), (T, 0.625), (V, 0.2), (T, -0.5), (V, 0.2), (T, 0.625), (V, 0.1), (T, -0.5), (V, 0.2), (T, 0.625)],
    steps, 2)

block:
  let
    (VFhl, K) = newIntegratorPair(updateVFhl, updateK)
    VFh = VFhl[0]
    VFl = VFhl[1]

  runtest("2*Test0", 1):
    var
      VK = mkTest0(K, VFh, 2)
      Vl = mkTest0(K, VFl, 1)
    newParallelEvolution(VK, Vl)

  runtest("2*Test1", 1):
    var
      VK = mkTest1(K, VFh, 2)
      Vl = mkTest1(K, VFl, 1)
    newParallelEvolution(VK, Vl)

  runtest("Test0+Test2", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl = mkTest0(K, VFl, 1)
    newParallelEvolution(VK, Vl)

  runtest("Test1+Test2", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl = mkTest1(K, VFl, 1)
    newParallelEvolution(VK, Vl)

  runtest("Test2+Test3", 1):
    var
      VK = mkTest3(K, VFh, 1)
      Vl = mkTest2(K, VFl, 1)
    newParallelEvolution(VK, Vl)

block:
  let
    (VFhl01, K) = newIntegratorPair(updateVFhl01, updateK)
    VFh = VFhl01[0]
    VFl0 = VFhl01[1]
    VFl1 = VFhl01[2]

  runtest("Test0+Test1+Test2", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest0(K, VFl0, 1)
      Vl1K = mkTest1(K, VFl1, 1)
    newParallelEvolution(VK, Vl0K, Vl1K)

  runtest("Test1+Test2+Test3", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest3(K, VFl0, 1)
      Vl1K = mkTest1(K, VFl1, 1)
    newParallelEvolution(VK, Vl0K, Vl1K)

  runtest("Test0^4+Test1+Test2", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest0(K, VFl0, 4)
      Vl1K = mkTest1(K, VFl1, 1)
    newParallelEvolution(VK, Vl0K, Vl1K)

  runtest("Test0*Test2+Test1", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest0(VK, VFl0, 1)
      Vl1K = mkTest1(K, VFl1, 1)
    newParallelEvolution(Vl0K, Vl1K)

  runtest("Test0*Test2+Test1", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest0(VK, VFl0, 1)
      Vl1K = mkTest1(K, VFl1, 8)
    newParallelEvolution(Vl0K, Vl1K)

  runtest("Test3^4+Test1+Test2", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest3(K, VFl0, 4)
      Vl1K = mkTest1(K, VFl1, 1)
    newParallelEvolution(VK, Vl0K, Vl1K)

  runtest("Test3*Test2+Test1", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest3(VK, VFl0, 1)
      Vl1K = mkTest1(K, VFl1, 1)
    newParallelEvolution(Vl0K, Vl1K)

  runtest("Test3*Test2+Test1", 1):
    var
      VK = mkTest2(K, VFh, 1)
      Vl0K = mkTest3(VK, VFl0, 1)
      Vl1K = mkTest1(K, VFl1, 8)
    newParallelEvolution(Vl0K, Vl1K)
