import mdevolve
import math, os

let Verb = existsenv("VerbTest")
# Simple 2D system
type A = array[2,float]
type S = object
  x:A
  rs:A  # Share some calculations
  old:bool  # We update rs if old
# State save/load.
proc save(x:var A, y:S) = x = y.x
proc load(y:var S, x:A) =
  y.x = x
  y.old = true
# Works like an array.  Need to remember to update `old`.
template `[]`(s:S,i:untyped):untyped = s.x[i]
proc updateR(s:var S) =
  if s.old:
    let
      x0 = s[0]
      x1 = s[1]
      r2 = x0*x0+x1*x1
    s.rs[0] = r2
    s.rs[1] = sqrt r2
    s.old = false
proc r(s:var S):float =
  s.updateR
  s.rs[1]
proc r2(s:var S):float =
  s.updateR
  s.rs[0]
# Harmonic potential and its derivative.
proc Vh(s:var S):float = 0.01 * s.r2 / 2.0
proc Fh(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
  result[0] = -0.01 * x0
  result[1] = -0.01 * x1
# Separating the log potential for this demo
proc Vl0(s:var S):float = -ln(s.r/(s.r+1))
proc Fl0(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    d = 1/(s.r2*(1+s.r))
  result[0] = x0 * d
  result[1] = x1 * d
proc Vl1(s:var S):float = -ln(s.r+1)
proc Fl1(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    d = 1/(s.r2+s.r)
  result[0] = x0 * d
  result[1] = x1 * d
proc K(p:A):float =
  let
    p0 = p[0]
    p1 = p[1]
    k2 = p0*p0+p1*p1
  k2 / 2
proc dK(p:A):A = p

# Prepare a system for our MD evolution routines.
# Assemble it as TestSys
type TestSys = object
  x:S
  p:A
proc Vh(s:var TestSys):float = Vh s.x
proc Fh(s:var TestSys):A = Fh s.x
proc Vl0(s:var TestSys):float = Vl0 s.x
proc Fl0(s:var TestSys):A = Fl0 s.x
proc Vl1(s:var TestSys):float = Vl1 s.x
proc Fl1(s:var TestSys):A = Fl1 s.x
proc K(s:TestSys):float = K s.p
proc dK(s:TestSys):A = dK s.p
proc status(s:var TestSys):string =
  result = "x: " & $s.x[0] & " " & $s.x[1]
  result &= "  p: " & $s.p[0] & " " & $s.p[1]
  let
    vc = s.Vh
    vl0 = s.Vl0
    vl1 = s.Vl1
    k = s.K
    h = vc + vl0 + vl1 + k
  result &= "  Vh: " & $vc & "  Vl0: " & $vl0 & "  Vl1: " & $vl1 & "  K: " & $k & "  H: " & $h

# Declare the system and prepare for the evolution.
var s:TestSys
# Basic symplectic update
proc updateK(t:float) =
  let d = dK s
  s.x[0] += t*d[0]
  s.x[1] += t*d[1]
  s.x.old = true
template mkUpdateV(force:untyped):untyped =
  proc `updateV force`(t:float) =
    let d = force s
    s.p[0] += t*d[0]
    s.p[1] += t*d[1]
  proc `updateVfg force`(t:float) =
    let d = force s
    s.x[0] += t*d[0]
    s.x[1] += t*d[1]
    s.x.old = true
mkUpdateV Fh
mkUpdateV Fl0
mkUpdateV Fl1
var xsave:A  # location to save x

let
  t = 0.1
  n = 12
template runtest0(label, ss, genHMCbody:untyped) =
  var H = genHMCbody
  echo "######## BEGIN ",label," ########"
  var dHs,dt,er = newseq[float]()
  for k in 0..2:
    let steps = ss*2^k
    if Verb: echo "# steps: ",steps
    H.steps = steps
    s.x.x = [1.0,0]
    s.x.old = true
    s.p = [-1.5,0.5]
    let h0 = s.K + s.Vh + s.Vl0 + s.Vl1
    for i in 0..<n:
      if Verb: echo (t*i.float)," ",s.status
      H.evolve t
    if Verb: echo (t*n.float)," ",s.status
    dt.add(t / steps.float)
    dHs.add(s.K + s.Vh + s.Vl0 + s.Vl1 - h0)
    s.p[0] = -s.p[0]
    s.p[1] = -s.p[1]
    for i in 0..<n: H.evolve t
    if Verb: echo "# reversed to origin: ",s.status
    let h1 = s.K + s.Vh + s.Vl0 + s.Vl1
    er.add((h1-h0)/max(abs(h0),abs(h1)))
    if Verb: echo "\n"
  var m,e = 0.0
  for i in 0..<dt.len:
    if i > 0:  # naive estimate of the error scaling
      let
        a = ln(abs(dHs[i]/dHs[i-1]))/ln(dt[i]/dt[i-1])
        d = a - m
        di = d / i.float
      m += di
      e += d * di * float(i-1)
      echo "# naive integrator error scale factor: ",a
    echo "# dt,dH,er: ",dt[i]," ",dHs[i]," ",er[i]
    if abs(er[i]) > 1e-12:
      echo "# Error: failed reversibility test."
      quit 1
  echo "# Estimate naive error scale factor: ",m," +/- ",sqrt(e/float((dt.len-1)*(dt.len-2)))
  echo "######## END ",label," ########"
template runtest(label, s, genHMCbody:untyped) =
  proc run {.gensym.} = runtest0(label, s, genHMCbody)
  run()
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
