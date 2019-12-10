import mdevolve
import sys2d, utils
import math

block:
  let
    (VF, K) = newIntegratorPair(updateVFGa, updateK)

  runtest("CreutzGocksch:Omelyan4MN2F1GV:a", 1, @[8.269949169914526e-10, 1.367173041444403e-11, 2.136069099378801e-13]):
    mkCreutzGocksch mkOmelyan4MN2F1GV(T=K, V=VF[0])

  runtest("CreutzGocksch:Omelyan4MN2F2GV:a", 1, @[8.520317784643794e-11, -2.704347856763434e-11, -2.243760732767441e-12]):
    mkCreutzGocksch mkOmelyan4MN2F2GV(T=K, V=VF[0])

  runtest("CreutzGocksch:Omelyan4MN3F1GP:a", 1, @[-4.040123791071437e-10, -6.544098596350523e-12, -1.008082506359642e-13]):
    mkCreutzGocksch mkOmelyan4MN3F1GP(T=K, V=VF[0])

  runtest("CreutzGocksch:Omelyan4MN5F2GP:a", 1, @[-7.361622422763503e-11, -1.173283692423865e-12, -1.554312234475219e-14]):
    mkCreutzGocksch mkOmelyan4MN5F2GP(T=K, V=VF[0])

  runtest("Omelyan8CK4P4:Omelyan4MN2F1GV:a", 1, @[3.946776827823406e-05, -1.270631826422175e-07, -9.163780845256042e-13], 4):
    mkOmelyan8CK4P4 mkOmelyan4MN2F1GV(T=K, V=VF[0])

  runtest("Omelyan8CK4P4:Omelyan4MN2F2GV:a", 1, @[2.743373372715574e-05, -9.56785606298638e-08, 8.520784078314136e-11], 4):
    mkOmelyan8CK4P4 mkOmelyan4MN2F2GV(T=K, V=VF[0])

  runtest("Omelyan8CK4P4:Omelyan4MN3F1GP:a", 1, @[-1.278117947078883e-05, 4.891433658116284e-08, 4.529709940470639e-13], 4):
    mkOmelyan8CK4P4 mkOmelyan4MN3F1GP(T=K, V=VF[0])

  runtest("Omelyan8CK4P4:Omelyan4MN3F1GV:a", 1, @[8.658901908731309e-05, -9.466855366291327e-08, 6.22449469744879e-09], 4):
    mkOmelyan8CK4P4 mkOmelyan4MN3F1GV(T=K, V=VF[0])

  runtest("Omelyan8CK4P4:Omelyan4MN5F2GP:a", 1, @[-3.430390906222769e-06, 1.451429954002492e-08, 9.747758156208874e-14], 4):
    mkOmelyan8CK4P4 mkOmelyan4MN5F2GP(T=K, V=VF[0])

  runtest("Omelyan10CK4P6:Omelyan4MN2F1GV:a", 1, @[5.335483779589367e-05, -1.371238589342738e-08, -4.85500528668581e-11], 6):
    mkOmelyan10CK4P7 mkOmelyan4MN2F1GV(T=K, V=VF[0])

  runtest("Omelyan10CK4P6:Omelyan4MN2F2GV:a", 1, @[3.9174368965611e-05, -4.639084694169071e-08, -5.618319143252393e-10], 6):
    mkOmelyan10CK4P7 mkOmelyan4MN2F2GV(T=K, V=VF[0])

  runtest("Omelyan10CK4P6:Omelyan4MN3F1GP:a", 1, @[-2.677615667101563e-05, 7.084718212091445e-08, 8.63524807215299e-11], 6):
    mkOmelyan10CK4P7 mkOmelyan4MN3F1GP(T=K, V=VF[0])

  runtest("Omelyan10CK4P6:Omelyan4MN3F1GV:a", 1, @[-0.0001108941197380808, -2.381718461563764e-06, -3.792836045235504e-08], 6):
    mkOmelyan10CK4P7 mkOmelyan4MN3F1GV(T=K, V=VF[0])

  runtest("Omelyan10CK4P6:Omelyan4MN5F2GP:a", 1, @[-4.272498320023743e-06, 1.239832503330263e-08, 1.030375784694115e-11], 6):
    mkOmelyan10CK4P7 mkOmelyan4MN5F2GP(T=K, V=VF[0])

  runtest("Omelyan12CK4P12:Omelyan4MN2F1GV:a", 1, @[-1.081068243458105e-05, -2.380192976048079e-08, -3.328448627826219e-13], 6):
    mkOmelyan12CK4P12 mkOmelyan4MN2F1GV(T=K, V=VF[0])

  runtest("Omelyan12CK4P12:Omelyan4MN2F2GV:a", 1, @[-7.561781960596647e-06, -1.786506942380583e-08, -1.332738364112629e-10], 6):
    mkOmelyan12CK4P12 mkOmelyan4MN2F2GV(T=K, V=VF[0])

  runtest("Omelyan12CK4P12:Omelyan4MN3F1GP:a", 1, @[-2.881535328125295e-07, -1.196932331026801e-09, 7.749356711883593e-13], 6):
    mkOmelyan12CK4P12 mkOmelyan4MN3F1GP(T=K, V=VF[0])

  runtest("Omelyan12CK4P12:Omelyan4MN3F1GV:a", 1, @[2.232111008138027e-05, 7.701530568837711e-08, -9.776101705938345e-09], 6):
    mkOmelyan12CK4P12 mkOmelyan4MN3F1GV(T=K, V=VF[0])

  runtest("Omelyan12CK4P12:Omelyan4MN5F2GP:a", 1, @[-7.941368340702581e-08, -1.910205327249059e-11, 1.170175067954915e-13], 6):
    mkOmelyan12CK4P12 mkOmelyan4MN5F2GP(T=K, V=VF[0])

block:
  let
    (VF, K) = newIntegratorPair(updateVFGa2, updateK)

  runtest("Omelyan10CK6P4:Omelyan6MN5F3GP:a2", 1, @[5.222074905164575e-06, -1.85372432959241e-08, -2.122857445385762e-11], 6):
    mkOmelyan10CK6P4 mkOmelyan6MN5F3GP(T=K, V=VF[0])

  runtest("Omelyan12CK6P7:Omelyan6MN5F3GP:a2", 1, @[-6.326127325095854e-06, -8.995630906838414e-09, 1.813038608133866e-11], 6):
    mkOmelyan12CK6P7 mkOmelyan6MN5F3GP(T=K, V=VF[0])

block:
  let
    (VF, K) = newIntegratorPair(updateVFG, updateK)

  runtest("CreutzGocksch:Omelyan6MN5F3GP", 1, @[0.0002423680431111741, -3.302853980091669e-06, -2.71259192885509e-09], 6):
    mkCreutzGocksch mkOmelyan6MN5F3GP(T=K, V=VF[0])

  runtest("CreutzGocksch:Omelyan8S11F11GP", 1, @[-1.088491866219066e-05, 4.575905632187016e-08, 6.981792921578744e-11], 6):
    mkCreutzGocksch mkOmelyan8S11F11GP(T=K, V=VF[0])

  runtest("Omelyan10CK6P4:Omelyan6MN5F3GP", 1, @[5.228193666440006e-06, -1.852199327245785e-08, -2.123057285530194e-11], 6):
    mkOmelyan10CK6P4 mkOmelyan6MN5F3GP(T=K, V=VF[0])

  runtest("Omelyan12CK6P7:Omelyan6MN5F3GP", 1, @[-6.325788849848024e-06, -8.995665101707573e-09, 1.813127425975836e-11], 6):
    mkOmelyan12CK6P7 mkOmelyan6MN5F3GP(T=K, V=VF[0])

# no gradient versions
block:
  let
    (V, K) = newIntegratorPair(updateVF, updateK)

  runtest("CreutzGocksch:Omelyan6MN7FV", 1, @[-0.0003826625796956584, 5.974412358611403e-06, 5.107910983070951e-09], 6):
    mkCreutzGocksch mkOmelyan6MN7FV(T=K, V=V)

  runtest("Omelyan10CK6P4:Omelyan6MN7FV", 1, @[-2.096292517528298e-05, 4.324073721306831e-08, 5.864109198228107e-11], 6):
    mkOmelyan10CK6P4 mkOmelyan6MN7FV(T=K, V=V)

  runtest("Omelyan12CK6P7:Omelyan6MN7FV", 1, @[1.340377017866601e-05, 9.883400320376268e-09, -2.583933067512589e-11], 6):
    mkOmelyan12CK6P7 mkOmelyan6MN7FV(T=K, V=V)
