type Concur[X,Y] = object
  x:X
  y:Y
  n:int
proc evolve*(S:Concur, t:float) =
  # Calculate forces at the same time
  S.x.evolve t
  S.y.evolve t

type Leapfrog[X,Y] = object
  x:X
  y:Y
  n:int
proc evolve*(S:LeapFrog, t:float) =
  let
    h = t / S.n.float
    hh = 0.5*h
  S.x.evolve hh
  S.y.evolve h
  for i in 1..<S.n:
    S.x.evolve h
    S.y.evolve h
  S.x.evolve hh

type SW92[X,Y] = object
  ## Sexton & Weingarten (1992)
  x:X
  y:Y
  n:int
proc evolve*(S:SW92, t:float) =
  ## Sexton & Weingarten (1992)
  let
    h = t / S.n.float
    h2 = 0.5*h
    h3 = h/3.0
    h23 = (2.0/3.0)*h
    h6 = h/6.0
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

type FGYin11[X,Y] = object
  ## Force Gradient Integrator, H. Yin (2011)
  x:X
  y:Y
  n:int
type FGupdate[X] = object
  s:X
proc evolve*(S:FGYin11, t:float) =
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

when isMainModule:
  import math, os

  let Verb = existsenv("VerbTest")
  # Simple 2D system
  type A = array[2,float]
  # Coulomb and Harmonic potential
  proc Vc(x:A):float =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
    -1 / r + r2 / 2.0
  proc Fc(x:A):A =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
      r3 = r*r2
    result[0] = -(x0 + x0 / r3)
    result[1] = -(x1 + x1 / r3)
  # Separating the log potential for this demo
  proc Vl0(x:A):float =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
    -ln(r/(r+1))
  proc Fl0(x:A):A =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
      d = 1/(r2*(1+r))
    result[0] = x0 * d
    result[1] = x1 * d
  proc Vl1(x:A):float =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
    -ln(r+1)
  proc Fl1(x:A):A =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
      d = 1/(r2+r)
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
    x:A
    p:A
  proc Vc(s:TestSys):float = Vc s.x
  proc Fc(s:TestSys):A = Fc s.x
  proc Vl0(s:TestSys):float = Vl0 s.x
  proc Fl0(s:TestSys):A = Fl0 s.x
  proc Vl1(s:TestSys):float = Vl1 s.x
  proc Fl1(s:TestSys):A = Fl1 s.x
  proc K(s:TestSys):float = K s.p
  proc dK(s:TestSys):A = dK s.p
  proc status(s:TestSys):string =
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
  proc evolve(s:MDVc|MDVl0|MDVl1, t:float) =
    let d = delta s
    s.s.p[0] += t*d[0]
    s.s.p[1] += t*d[1]

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
    let f = delta s
    s.s.p[0] += t*f[0]
    s.s.p[1] += t*f[1]
    s.s.x = x

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
      s.x = [1.0,0]
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
