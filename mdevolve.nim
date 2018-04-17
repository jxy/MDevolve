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

when isMainModule:
  import math
  type
    A = array[2,float]
    Kepler = object
      x:A
      p:A
  proc V(s:Kepler):float =
    let
      x0 = s.x[0]
      x1 = s.x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
    -1 / r
  proc F(s:Kepler):A =
    let
      x0 = s.x[0]
      x1 = s.x[1]
      r2 = x0*x0+x1*x1
      r = sqrt r2
      r3 = r*r2
    result[0] = - x0 / r3
    result[1] = - x1 / r3
  proc K(s:Kepler):float =
    let
      p0 = s.p[0]
      p1 = s.p[1]
      k2 = p0*p0+p1*p1
    k2 / 2
  proc dK(s:Kepler):A = s.p
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
  proc status(s:Kepler):string =
    result = "x: " & $s.x[0] & " " & $s.x[1]
    result &= "  p: " & $s.p[0] & " " & $s.p[1]
    let
      v = s.V
      k = s.K
      h = v + k
    result &= "  V: " & $v & "  K: " & $k & "  H: " & $h
  var s:Kepler
  let
    t = 1.0
    n = 10
  template runtest(algo, label, ss:untyped) =
    echo "######## BEGIN ",label," ########"
    var dHs,dt,er = newseq[float]()
    for steps in ss:
      # echo "# steps: ",steps
      var H = algo[MDK,MDV](n:steps, x:MDK(s:s.addr), y:MDV(s:s.addr))
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
