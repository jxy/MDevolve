import macros, os
import math

var Verbose* = existsenv("VerbTest")

macro verbose(xs: varargs[untyped]):untyped =
  var echo = newcall(bindsym"echo")
  for x in xs: echo.add x
  result = quote do:
    if Verbose: `echo`

let
  tDef = 0.1
  nDef = 12
  xt = [ 0.0979503986562138, 1.104054119512884]
  pt = [-0.3096105574933478, 1.614829451689578]
template runtest0(label, ss, deltaHs, n, t, genHMCbody:untyped) =
  var H = genHMCbody
  let delHs:seq[float] = deltaHs
  echo "######## BEGIN ",label," ########"
  when mdevolveVersion > 0: echo H
  var dHs,dt,er = newseq[float]()
  var failed = 0
  for k in 0..2:
    record.setlen(0)
    let steps = ss*2^k
    verbose "# steps: ",steps
    H.steps = steps
    s.x.x = [1.0,0]
    s.x.old = true
    s.p = [-1.5,0.5]
    let h0 = s.K + s.Vh + s.Vl
    for i in 0..<n:
      verbose (t*i.float)," ",s.status
      H.evolve t
    when mdevolveVersion > 0: H.finish
    verbose (t*n.float)," ",s.status
    # echo "# x: ",s.x.x,"  p: ",s.p
    if k == 2:
      var dx,dp:array[2,float]
      let dtol = 1E-3
      for i in 0..1:
        dx[i] = xt[i] - s.x.x[i]
        dp[i] = pt[i] - s.p[i]
        if abs(dx[i]) > dtol or abs(dp[i]) > dtol:
          failed += 1
          echo "Error: large error."
      echo "# dx: ",dx,"  dp: ",dp
    dt.add(t / steps.float)
    dHs.add(s.K + s.Vh + s.Vl0 + s.Vl1 - h0)
    if not record.statusCheck:
      when mdevolveVersion > 0:
        failed += 1
    record.setlen(0)
    s.p[0] = -s.p[0]
    s.p[1] = -s.p[1]
    for i in 0..<n: H.evolve t
    when mdevolveVersion > 0: H.finish
    verbose "# reversed to origin: ",s.status
    let h1 = s.K + s.Vh + s.Vl
    er.add((h1-h0)/max(abs(h0),abs(h1)))
    verbose "\n"
    echo "# dt,dH,er: ",dt[^1]," ",dHs[^1]," ",er[^1]
    if not record.statusCheck:
      when mdevolveVersion > 0:
        failed += 1
  var m,e = 0.0
  for i in 0..<dt.len:
    if i > 0:  # naive estimate of the error scaling
      let
        a = ln(abs(dHs[i]/dHs[i-1]))/ln(dt[i]/dt[i-1])
        d = a - m
        di = d / i.float
      m += di
      e += d * di * float(i-1)
      echo "# naive integrator error scale factor (", i-1, "->", i, "): ",a
    if abs(er[i]) > 1e-12:
      echo "# Error: failed reversibility test."
      failed += 1
  echo "# Estimate naive error scale factor: ",m," +/- ",sqrt(e/float((dt.len-1)*(dt.len-2)))
  for i in 0..<delHs.len:
    if abs(dHs[i]-delHs[i])/max(abs(dHs[i]),abs(delHs[i])) > 1e-9:
      echo "# Error: expecting delHs[", i, "]=",delHs[i]
      echo "#        Actual dHs[",i,"]=",dHs[i]
      echo "#        Failed reproducing dH"
      failed += 1
  echo "######## END ",label," ########"
  if failed > 0: quit failed
template runtest*(label, s, genHMCbody:untyped) =
  proc run {.gensym.} =
    let
      n = nDef
      t = tDef
    runtest0(label, s, @[], n, t, genHMCbody)
  run()
template runtest*(label, s, dH, genHMCbody:untyped) =
  proc run {.gensym.} =
    let
      n = nDef
      t = tDef
    runtest0(label, s, dH, n, t, genHMCbody)
  run()
template runtest*(label, s, dH, factor, genHMCbody:untyped) =
  proc run {.gensym.} =
    let
      n = nDef div factor
      t = nDef.float / n.float * tDef
    runtest0(label, s, dH, n, t, genHMCbody)
  run()
