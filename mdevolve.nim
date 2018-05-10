type Leapfrog[X,Y] = object
  x:X
  y:Y
  n:int
template mkLeapfrog*[X,Y](steps:int, algoX:X, algoY:Y):untyped =
  Leapfrog[type(algoX),type(algoY)](n:steps, x:algoX, y:algoY)
proc evolve*(S:LeapFrog, t:float) =
  mixin evolve
  let
    h = t / S.n.float
    hh = 0.5*h
  S.x.evolve hh
  S.y.evolve h
  for i in 1..<S.n:
    S.x.evolve h
    S.y.evolve h
  S.x.evolve hh

type S5[X,Y] = object
  x:X
  y:Y
  n:int
  l:float
template mkSW92*[X,Y](steps:int, algoX:X, algoY:Y):untyped =
  ## Sexton & Weingarten (1992)
  S5[type(algoX),type(algoY)](n:steps, x:algoX, y:algoY, l:1.0/6.0)
template mkOmelyan2MN*[X,Y](steps:int, T:X, V:Y,
    lambda = 0.1931833275037836):untyped =
  ## Omelyan et. al. (2002)
  S5[type(T),type(V)](n:steps, x:T, y:V, l:lambda)
proc evolve*(S:S5, t:float) =
  mixin evolve
  let
    h = t / S.n.float
    h2 = 0.5 * h
    h6 = S.l * h
    h3 = 2.0 * h6
    h23 = h - h3
  S.x.evolve h6
  S.y.evolve h2
  S.x.evolve h23
  S.y.evolve h2
  for i in 1..<S.n:
    S.x.evolve h3
    S.y.evolve h2
    S.x.evolve h23
    S.y.evolve h2
  S.x.evolve h6

type S9[X,Y] = object
  x:X
  y:Y
  n:int
  r,t,l:float
template mkOmelyan4MN4FP*[X,Y](steps:int, T:X, V:Y,
    rho = 0.1786178958448091,
    theta = -0.06626458266981843,
    lambda = 0.7123418310626056):untyped =
  ## Omelyan et. al. (2003)
  S9[type(T),type(V)](n:steps, x:T, y:V, r:rho, t:theta, l:lambda)
proc evolve*(S:S9, t:float) =
  mixin evolve
  let
    h = t / S.n.float
    a = S.r * h
    b = S.l * h
    c = S.t * h
    d = (0.5 - S.l) * h
    e = (1.0 - 2*(S.t+S.r)) * h
    a2 = 2 * a
  S.x.evolve a
  S.y.evolve b
  S.x.evolve c
  S.y.evolve d
  S.x.evolve e
  S.y.evolve d
  S.x.evolve c
  S.y.evolve b
  for i in 1..<S.n:
    S.x.evolve a2
    S.y.evolve b
    S.x.evolve c
    S.y.evolve d
    S.x.evolve e
    S.y.evolve d
    S.x.evolve c
    S.y.evolve b
  S.x.evolve a

type S11[X,Y] = object
  x:X
  y:Y
  n:int
  r,t,l,m:float
template mkOmelyan4MN5FV*[X,Y](steps:int, V:X, T:Y,
    theta = 0.08398315262876693,
    rho = 0.2539785108410595,
    lambda = 0.6822365335719091,
    mu = -0.03230286765269967):untyped =
  ## Omelyan et. al. (2003)
  S11[type(V),type(T)](n:steps, x:V, y:T, r:rho, t:theta, l:lambda, m:mu)
proc evolve*(S:S11, t:float) =
  mixin evolve
  let
    h = t / S.n.float
    a = S.t * h
    b = S.r * h
    c = S.l * h
    d = S.m * h
    e = (0.5 - (S.l+S.t)) * h
    f = (1.0 - 2*(S.m+S.r)) * h
    a2 = 2 * a
  S.x.evolve a
  S.y.evolve b
  S.x.evolve c
  S.y.evolve d
  S.x.evolve e
  S.y.evolve f
  S.x.evolve e
  S.y.evolve d
  S.x.evolve c
  S.y.evolve b
  for i in 1..<S.n:
    S.x.evolve a2
    S.y.evolve b
    S.x.evolve c
    S.y.evolve d
    S.x.evolve e
    S.y.evolve f
    S.x.evolve e
    S.y.evolve d
    S.x.evolve c
    S.y.evolve b
  S.x.evolve a

type FGYin11[X,Y] = object
  ## Force Gradient Integrator, H. Yin (2011)
  x:X
  y:Y
  n:int
template mkFGYin11*[X,Y](steps:int, algoX:X, algoY:Y):untyped =
  FGYin11[type(algoX),type(algoY)](n:steps, x:algoX, y:algoY)
type FGupdate[X] = object
  s:X
proc evolve*(S:FGYin11, t:float) =
  mixin evolve
  ## Force Gradient Integrator, H. Yin (2011)
  let
    h = t / S.n.float
    h2 = 0.5*h
    h3 = h/3.0
    h23 = (2.0/3.0)*h
    h6 = h/6.0
  S.x.evolve h6
  S.y.evolve h2
  S.x.fgupdate.evolve h23
  S.y.evolve h2
  for i in 1..<S.n:
    S.x.evolve h3
    S.y.evolve h2
    S.x.fgupdate.evolve h23
    S.y.evolve h2
  S.x.evolve h6

type Concur[X,Y] = object
  x:X
  y:Y
proc evolve*(S:Concur, t:float) =
  mixin evolve
  # Calculate forces at the same time
  S.x.evolve t
  S.y.evolve t

type Combine[X,Y] = object
  x:X
  y:Y
template sharedEvolver(X,Y:untyped):untyped = discard  # FIXME
iterator combinedSteps(S,T:untyped):float = discard  # FIXME
proc evolve*(S:Combine, t:float) =
  mixin evolve
  # Combine two integrators.
  # Requires both share exactly one base updater,
  # i.e. exp(t T)), that updates X with P,
  # and other updaters unique to each that updat P with X,
  # i.e. exp(t Sx), exp(t Sy), with [Sx,Sy]=0,
  # such that the combined evolution integrate the Hamiltonian
  # system of T+Sx+Sy.
  # Only works for positive time steps in the shared updater, T.

  # Discover which type is the shared one
  type T = sharedEvolver(S.X, S.Y)

  for d in combinedSteps(S, T):
    S.x.evolveStep(t,d)
    S.y.evolveStep(t,d)

when isMainModule:
  import math, os

  let Verb = existsenv("VerbTest")
  # Simple 2D system
  type A = array[2,float]
  type S = object
    x:A
    rs:A  # Share some calculations
    old:bool  # We update rs if old
  # To keep our old code
  template `[]`(s:S,i:untyped):untyped = s.x[i]
  proc updateR(s:var S) =
    if s.old:
      let
        x0 = s.x[0]
        x1 = s.x[1]
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
  # Coulomb and Harmonic potential
  proc Vc(s:var S):float = -1 / s.r + s.r2 / 2.0
  proc Fc(s:var S):A =
    let
      x0 = s[0]
      x1 = s[1]
      r3 = s.r*s.r2
    result[0] = -(x0 + x0 / r3)
    result[1] = -(x1 + x1 / r3)
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

  # Assemble it as TestSys
  type TestSys = object
    x:S
    p:A
  proc Vc(s:var TestSys):float = Vc s.x
  proc Fc(s:var TestSys):A = Fc s.x
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
      vc = s.Vc
      vl0 = s.Vl0
      vl1 = s.Vl1
      k = s.K
      h = vc + vl0 + vl1 + k
    result &= "  Vc: " & $vc & "  Vl0: " & $vl0 & "  Vl1: " & $vl1 & "  K: " & $k & "  H: " & $h

  # Basic symplectic update
  type
    MDK = object
      s:ptr TestSys
    MDVc = object
      s:ptr TestSys
    MDVl0 = object
      s:ptr TestSys
    MDVl1 = object
      s:ptr TestSys
  # For this simple model we overload delta (-dV/dx or dK/dp)
  proc delta(s:MDK):A = dK s.s[]
  proc delta(s:MDVc):A = Fc s.s[]
  proc delta(s:MDVl0):A = Fl0 s.s[]
  proc delta(s:MDVl1):A = Fl1 s.s[]
  proc evolve(s:MDK, t:float) =
    let d = delta s
    s.s.x[0] += t*d[0]
    s.s.x[1] += t*d[1]
    s.s.x.old = true
  proc evolve(s:MDVc|MDVl0|MDVl1, t:float) =
    let d = delta s
    s.s.p[0] += t*d[0]
    s.s.p[1] += t*d[1]
    s.s.x.old = true

  # Force gradient update
  proc fgupdate(s:MDVc|MDVl0|MDVl1):auto = FGupdate[type(s)](s:s)
  proc evolve(s:FGupdate, t:float) =
    # t = 2/3 tau (Yin, 2011)
    let
      s = s.s  # Nim allows us to do this.
      d = delta s
      h = (3.0/32.0)*t*t  # 1r24=3r32**:2r3
    # Approximation is valid when h*d is small.
    let x = s.s.x
    s.s.x[0] += h*d[0]
    s.s.x[1] += h*d[1]
    s.s.x.old = true
    let f = delta s
    s.s.p[0] += t*f[0]
    s.s.p[1] += t*f[1]
    s.s.x = x
    s.s.x.old = true

  var s:TestSys
  let
    t = 1.0
    n = 10
  template runtest(algo0, algo1, algo2, label, ss:untyped) =
    # algo0/1/2 are the outer, middle, and inner integrators
    echo "######## BEGIN ",label," ########"
    var dHs,dt,er = newseq[float]()
    for steps in ss:
      if Verb: echo "# steps: ",steps
      let
        VK = algo2[MDVc,MDK](n:4, x:MDVc(s:s.addr), y:MDK(s:s.addr))
        Vl0VK = algo1[MDVl0,type(VK)](n:4, x:MDVl0(s:s.addr), y:VK)
        H = algo0[MDVl1,type(Vl0VK)](n:steps, x:MDVl1(s:s.addr), y:Vl0VK)
      s.x.x = [1.0,0]
      s.x.old = true
      s.p = [0.0,0.7]
      let h0 = s.K + s.Vc + s.Vl0 + s.Vl1
      for i in 0..<n:
        if Verb: echo (t*i.float)," ",s.status
        H.evolve t
      if Verb: echo (t*n.float)," ",s.status
      dt.add(t / steps.float)
      dHs.add(s.K + s.Vc + s.Vl0 + s.Vl1 - h0)
      for i in 0..<n:
        H.evolve(-t)
      if Verb: echo "# reversed to origin: ",s.status
      let h1 = s.K + s.Vc + s.Vl0 + s.Vl1
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
  runtest(Leapfrog, Leapfrog, Leapfrog, "Leapfrog", [20,40,80,160,320,640])
  runtest(SW92, SW92, SW92, "SW92", [20,40,80,160,320,640])
  runtest(FGYin11, SW92, Leapfrog, "FG/SW/LF", [20,40,80,160,320,640])
  runtest(FGYin11, FGYin11, FGYin11, "FGYin11", [5,10,20,40,80,160])
