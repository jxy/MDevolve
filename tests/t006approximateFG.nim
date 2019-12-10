import mdevolve
import sys2d, utils
import math

block:
  let
    (VF, K) = newIntegratorPair(updateVFGa, updateK)
    V = VF[0]

  # cf. t003forcegradient.nim

  runtest("Omelyan4MN2F1GV:a", 1, @[-3.006083422185668e-07, -1.901979396201625e-08, -1.192570264763049e-09]):
    mkOmelyan4MN2F1GV(T=K, V=V)

  runtest("Omelyan4MN3F1GP:a", 1, @[-3.301239726027916e-08, -2.11963513407909e-09, -1.334152788246001e-10]):
    mkOmelyan4MN3F1GP(T=K, V=V)

  runtest("Omelyan4MN3F1GV:a", 1, @[-5.439543726559748e-08, -3.350519861555767e-09, -2.085482897484781e-10]):
    mkOmelyan4MN3F1GV(T=K, V=V)

  runtest("Omelyan4MN5F2GP:a", 1, @[-2.17127864132749e-07, 2.423135025964029e-10, 1.433519969396002e-11], 2):
    mkOmelyan4MN5F2GP(T=K, V=V)

  runtest("Omelyan6MN5F3GP:a", 1, @[-1.350265969257691e-06, -2.746947558307511e-09, -1.847129116328006e-10], 2):
    mkOmelyan6MN5F3GP(T=K, V=V)

  runtest("Omelyan6MN5F5GP:a", 1, @[-1.131957490407842e-06, 1.556658180135173e-08, 9.752727514467097e-10], 2):
    mkOmelyan6MN5F5GP(T=K, V=V)

  runtest("Omelyan8S11F11GP:a", 1, @[-0.001140020893141358, 7.290049643948748e-07, 1.28624142448075e-08], 4):
    mkOmelyan8S11F11GP(T=K, V=V)

block:
  let
    (VF, K) = newIntegratorPair(updateVFGa2, updateK)
    V = VF[0]

  # cf. t003forcegradient.nim

  runtest("Omelyan4MN2F1GV:a2", 1, @[-3.099900278691337e-07, -1.96140153008173e-08, -1.229832236049333e-09]):
    mkOmelyan4MN2F1GV(T=K, V=V)

  runtest("Omelyan4MN3F1GP:a2", 1, @[-3.609571064266959e-08, -2.315034830502327e-09, -1.456721410164619e-10]):
    mkOmelyan4MN3F1GP(T=K, V=V)

  runtest("Omelyan4MN3F1GV:a2", 1, @[-8.138008666946916e-08, -4.991609792170948e-09, -3.103952650462816e-10]):
    mkOmelyan4MN3F1GV(T=K, V=V)

  runtest("Omelyan4MN5F2GP:a2", 1, @[-2.215046071007976e-07, -3.967937090010309e-11, -3.362199407774824e-12], 2):
    mkOmelyan4MN5F2GP(T=K, V=V)

  runtest("Omelyan6MN5F3GP:a2", 1, @[-1.326181828353867e-06, 2.249442854207473e-10, 3.72613051524695e-12], 2):
    mkOmelyan6MN5F3GP(T=K, V=V)

  runtest("Omelyan6MN5F5GP:a2", 1, @[-1.278005976601548e-06, 2.218860650771148e-10, 3.682609772681644e-12], 2):
    mkOmelyan6MN5F5GP(T=K, V=V)

  runtest("Omelyan8S11F11GP:a2", 1, @[-7.098181613462984e-05, 1.018743727154714e-07, -8.415890206947552e-11], 4):
    mkOmelyan8S11F11GP(T=K, V=V)
