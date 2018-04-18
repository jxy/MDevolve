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
  import math

  # Simple 2D system
  type A = array[2,float]
  proc V(x:A):float =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
    -1 / r
  proc F(x:A):A =
    let
      x0 = x[0]
      x1 = x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
      r3 = r*r2
    result[0] = - x0 / r3
    result[1] = - x1 / r3
  proc K(p:A):float =
    let
      p0 = p[0]
      p1 = p[1]
      k2 = p0*p0+p1*p1
    k2 / 2
  proc dK(p:A):A = p

  # Assemble it as Kepler
  type Kepler = object
    x:A
    p:A
  proc V(s:Kepler):float = V s.x
  proc F(s:Kepler):A = F s.x
  proc K(s:Kepler):float = K s.p
  proc dK(s:Kepler):A = dK s.p
  proc status(s:Kepler):string =
    result = "x: " & $s.x[0] & " " & $s.x[1]
    result &= "  p: " & $s.p[0] & " " & $s.p[1]
    let
      v = s.V
      k = s.K
      h = v + k
    result &= "  V: " & $v & "  K: " & $k & "  H: " & $h

  # Basic symplectic update
  type
    MDK = object
      s:ptr Kepler
    MDV = object
      s:ptr Kepler
  proc evolve(s:MDK, t:float) =
    let d = dK s.s[]
    s.s.x[0] += t*d[0]
    s.s.x[1] += t*d[1]
  proc evolve(s:MDV, t:float) =
    let f = F s.s[]
    s.s.p[0] += t*f[0]
    s.s.p[1] += t*f[1]

  # Force gradient update
  proc fgupdate(s:MDV):FGupdate[ptr Kepler] = FGupdate[ptr Kepler](s:s.s)
  proc evolve(s:FGupdate, t:float) =
    # t = 2/3 tau (Yin, 2011)
    let
      f = F s.s[]
      h = (3.0/32.0)*t*t  # 1r24=3r32**:2r3
    # Approximation is valid when h*f is small.
    var x = s.s.x
    x[0] += h*f[0]
    x[1] += h*f[1]
    let ff = F x
    s.s.p[0] += t*ff[0]
    s.s.p[1] += t*ff[1]

  var s:Kepler
  let
    t = 1.0
    n = 10
  template runtest(algo, label, ss:untyped) =
    echo "######## BEGIN ",label," ########"
    var dHs,dt,er = newseq[float]()
    for steps in ss:
      # echo "# steps: ",steps
      let H = algo[MDV,MDK](n:steps, x:MDV(s:s.addr), y:MDK(s:s.addr))
      s.x = [1.0,0]
      s.p = [0.0,0.75]
      let h0 = s.V + s.K
      for i in 0..<n:
        # echo (t*i.float)," ",s.status
        H.evolve t
      # echo (t*n.float)," ",s.status
      dt.add(t / steps.float)
      dHs.add(s.V + s.K - h0)
      for i in 0..<n:
        H.evolve(-t)
      # echo "# reversed to origin: ",s.status
      let h1 = s.V + s.K
      er.add((h1-h0)/h0)
      # echo "\n"
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
  runtest(Leapfrog, "Leapfrog", [50,100,200,400,800])
  runtest(SW92, "SW92", [50,100,200,400,800])
  runtest(FGYin11, "FGYin11", [25,50,100,200,400])
