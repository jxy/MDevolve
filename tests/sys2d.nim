import math, tables

# Simple 2D system
type A = array[2,float]
type S = object
  x*:A
  rs:A  # Share some calculations
  old*:bool  # We update rs if old
# State save/load.
proc save*(x:var A, y:S) = x = y.x
proc load*(y:var S, x:A) =
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
#[ Potentials and their gradients (using maxima):
(%i1) Vh: 0.01 * (x^2+y^2) / 2;
                                        2    2
(%o1)                           0.005 (y  + x )
(%i13) Vl: -log(sqrt(x^2+y^2));
                                       2    2
                                  log(y  + x )
(%o13)                          - ------------
                                       2
(%i29) F(S):=ratsubst(r^2,x^2+y^2,[diff(S,x), diff(S,y)]);
                              2   2    2
(%o29)      F(S) := ratsubst(r , x  + y , [diff(S, x), diff(S, y)])
(%i30) F(Vh);

rat: replaced 0.01 by 1/100 = 0.01

rat: replaced 0.01 by 1/100 = 0.01
                                    x    y
(%o30)                            [---, ---]
                                   100  100
(%i31) F(Vl);
                                    x     y
(%o31)                           [- --, - --]
                                     2     2
                                    r     r
(%i32) FG(S):= ratsubst(r^2,x^2+y^2,2*matrix([diff(S,x,2), diff(S,x,1,y,1)], [diff(S,y,1,x,1), diff(S,y,2)]).F(S));
                          2   2    2
(%o32) FG(S) := ratsubst(r , x  + y , 2
                          [    diff(S, x, 2)     diff(S, x, 1, y, 1) ]
                         ([                                          ] . F(S)))
                          [ diff(S, y, 1, x, 1)     diff(S, y, 2)    ]
(%i33) FG(Vh);

rat: replaced 0.01 by 1/100 = 0.01

rat: replaced 0.01 by 1/100 = 0.01

rat: replaced 2.0e-4 by 1/5000 = 2.0e-4

rat: replaced 2.0e-4 by 1/5000 = 2.0e-4
                                   [  x   ]
                                   [ ---- ]
                                   [ 5000 ]
(%o33)                             [      ]
                                   [  y   ]
                                   [ ---- ]
                                   [ 5000 ]
(%i34) FG(Vl);
                                   [   2 x ]
                                   [ - --- ]
                                   [    4  ]
                                   [   r   ]
(%o34)                             [       ]
                                   [   2 y ]
                                   [ - --- ]
                                   [    4  ]
                                   [   r   ]
(%i35) Vl0:-log(sqrt(x^2+y^2)/(sqrt(x^2+y^2)+1));
                                         2    2
                                   sqrt(y  + x )
(%o35)                     - log(-----------------)
                                       2    2
                                 sqrt(y  + x ) + 1
(%i44) assume(r>0);
(%o44)                              [r > 0]
(%i45) facts();
(%o45)                              [r > 0]
(%i50) factor(F(Vl0));
                                x             y
(%o50)                   [- ----------, - ----------]
                             2             2
                            r  (r + 1)    r  (r + 1)
(%i51) factor(FG(Vl0));
                              [   2 (2 r + 1) x ]
                              [ - ------------- ]
                              [     4        3  ]
                              [    r  (r + 1)   ]
(%o51)                        [                 ]
                              [   2 (2 r + 1) y ]
                              [ - ------------- ]
                              [     4        3  ]
                              [    r  (r + 1)   ]
(%i52) Vl1:-log(sqrt(x^2+y^2)+1);
                                       2    2
(%o52)                     - log(sqrt(y  + x ) + 1)
(%i53) F(Vl1);
                                  x         y
(%o53)                       [- ------, - ------]
                                 2         2
                                r  + r    r  + r
(%i54) factor(F(Vl1));
                                 x            y
(%o54)                    [- ---------, - ---------]
                             r (r + 1)    r (r + 1)
(%i55) factor(FG(Vl1));
                               [      2 x     ]
                               [ - ---------- ]
                               [            3 ]
                               [   r (r + 1)  ]
(%o55)                         [              ]
                               [      2 y     ]
                               [ - ---------- ]
                               [            3 ]
                               [   r (r + 1)  ]
]#
# Harmonic potential, its force (-S'), and force gradient ([V,[T,V]] = 2*S'.S''.d/dp).
proc Vh(s:var S):float = 0.01 * s.r2 / 2.0
proc Fh(s:S):A =
  let
    x0 = s[0]
    x1 = s[1]
  result[0] = -0.01 * x0
  result[1] = -0.01 * x1
proc FGh(s:S):A =
  # 2 * S' . S''
  let
    x0 = s[0]
    x1 = s[1]
  # d2S is diagonal
  result[0] = 2E-4 * x0
  result[1] = 2E-4 * x1
# A log potential
proc Vl(s:var S):float = -ln(s.r)
proc Fl(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    d = 1/s.r2
  result[0] = x0 * d
  result[1] = x1 * d
proc FGl(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    d = 1/s.r2
    d2 = -2.0*d*d
  result[0] = x0 * d2
  result[1] = x1 * d2
# Separating the log potential for this demo
proc Vl0(s:var S):float = -ln(s.r/(s.r+1))
proc Fl0(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    d = 1/(s.r2*(1+s.r))
  result[0] = x0 * d
  result[1] = x1 * d
proc FGl0(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    r = s.r
    r2 = s.r2
    t = r+1.0
    d = -2.0*(r+t)/(r2*r2*t*t*t)
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
proc FGl1(s:var S):A =
  let
    x0 = s[0]
    x1 = s[1]
    r = s.r
    t = r+1.0
    d = -2.0/(r*t*t*t)
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
  x*:S
  p*:A
proc Vh*(s:var TestSys):float = Vh s.x
proc Fh(s:var TestSys):A = Fh s.x
proc FGh(s:var TestSys):A = FGh s.x
proc Vl*(s:var TestSys):float = Vl s.x
proc Fl(s:var TestSys):A = Fl s.x
proc FGl(s:var TestSys):A = FGl s.x
proc Vl0*(s:var TestSys):float = Vl0 s.x
proc Fl0(s:var TestSys):A = Fl0 s.x
proc FGl0(s:var TestSys):A = FGl0 s.x
proc Vl1*(s:var TestSys):float = Vl1 s.x
proc Fl1(s:var TestSys):A = Fl1 s.x
proc FGl1(s:var TestSys):A = FGl1 s.x
proc V*(s:var TestSys):float = s.Vh + s.Vl
proc F(s:var TestSys):A =
  let
    h = Fh s
    l = Fl s
  result[0] = h[0] + l[0]
  result[1] = h[1] + l[1]
proc FG(s:var TestSys):A =
  let
    h = FGh s
    l = FGl s
  result[0] = h[0] + l[0]
  result[1] = h[1] + l[1]
proc K*(s:TestSys):float = K s.p
proc dK(s:TestSys):A = dK s.p
proc status*(s:var TestSys):string =
  result = "x: " & $s.x[0] & " " & $s.x[1]
  result &= "  p: " & $s.p[0] & " " & $s.p[1]
  let
    vc = s.Vh
    vl = s.Vl
    vl0 = s.Vl0
    vl1 = s.Vl1
    k = s.K
    h = vc + vl0 + vl1 + k
  result &= "  Vh: " & $vc & " Vl: " & $vl & "  Vl0: " & $vl0 & "  Vl1: " & $vl1 & "  K: " & $k & "  H: " & $h

type Record = seq[tuple[name:string, t:float]]
type Summary = TableRef[string, tuple[t:float, n:int, extra:int]]
func `$`*(s:Summary):string =
  result = "Evolution summary:\n"
  var n = 0
  for name, x in s.pairs:
    result &= "    " & name & ":  " & $x.n & " iter over time " & $x.t & " with " & $x.extra & " extra\n"
    n += x.n
  result &= "    Total: " & $n
proc statusCheck*(tape:Record):bool =
  result = true
  const eps = 1e-12
  var r:Summary
  r.new()
  var short:Record
  var prev = ""
  for x in tape:
    if x.t > -eps and x.t < eps:
      short.add x
    if x.name notin r:
      r[x.name] = (x.t, 1, 0)
    else:
      let z = r[x.name]
      let extra = if prev == x.name: z.extra+1 else: z.extra
      r[x.name] = (z.t+x.t, z.n+1, extra)
    prev = x.name
  echo r
  var extra = 0
  for x in r.values:
    extra += x.extra
  if short.len > 0:
    echo "# Warning: ", short.len, " vanishing integration steps: ", short
    result = false
  if extra > 0:
    echo "# Error: ", $extra, " extra calls should be avoided."
    result = false
var record*:Record

# Declare the system and prepare for the evolution.
var s*:TestSys
# Basic symplectic update
proc updateK*(t:float) =
  record.add(("K", t))
  let d = dK s
  s.x[0] += t*d[0]
  s.x[1] += t*d[1]
  s.x.old = true
template mkUpdateV(force:untyped):untyped =
  proc `updateV force`*(t:float) =
    record.add(("V_" & astToStr(force), t))
    let d = force s
    s.p[0] += t*d[0]
    s.p[1] += t*d[1]
  proc `updateVfg force`*(t:float) =
    record.add(("Vfg_" & astToStr(force), t))
    let d = force s
    s.x[0] += t*d[0]
    s.x[1] += t*d[1]
    s.x.old = true
mkUpdateV F
mkUpdateV Fh
mkUpdateV Fl
mkUpdateV Fl0
mkUpdateV Fl1
template mkUpdateG(force:untyped):untyped =
  proc `updateG force`*(t:float) =
    record.add(("G_" & astToStr(force), t))
    let d = force s
    s.p[0] += t*d[0]
    s.p[1] += t*d[1]
mkUpdateG FG
mkUpdateG FGh
mkUpdateG FGl
mkUpdateG FGl0
mkUpdateG FGl1

proc updateVFG*(ts:openarray[float], ss:openarray[float]) =
  record.add(("VFG(combined)", ts.sum+ss.sum))
  let
    t = ts[0]
    s = ss[0]
  if t != 0: updateVF t
  if s != 0: updateGFG s

proc updateVFhl*(ts:openarray[float]) =
  record.add(("VFHl(combined)", ts.sum))
  let
    th = ts[0]
    tl = ts[1]
  if th != 0: updateVFh th
  if tl != 0: updateVFl tl

proc updateVFGhl*(ts:openarray[float], ss:openarray[float]) =
  record.add(("VFGHl(combined)", ts.sum+ss.sum))
  let
    th = ts[0]
    tl = ts[1]
    sh = ss[0]
    sl = ss[1]
  if th != 0: updateVFh th
  if sh != 0: updateGFGh sh
  if tl != 0: updateVFl tl
  if sl != 0: updateGFGl sl

proc updateVFhl01*(ts:openarray[float]) =
  record.add(("VFHl01(combined)", ts.sum))
  let
    th = ts[0]
    tl0 = ts[1]
    tl1 = ts[2]
  if th != 0: updateVFh th
  if tl0 != 0: updateVFl0 tl0
  if tl1 != 0: updateVFl1 tl1

var xsave*:A  # location to save x

import mdevolve

proc updateVFGa*(ts,gs:openarray[float]) =
  record.add(("VFGa(combined)", ts.sum+gs.sum))
  let
    t = ts[0]
    g = gs[0]
  if g != 0:
    if t != 0:
      # Approximate the force gradient update with a Taylor expansion.
      let (tf,tg) = approximateFGcoeff(t,g)
      xsave.save s.x
      updateVfgF tg
      updateVF tf
      s.x.load xsave
    else:
      quit("Force gradient without the force update.")
  elif t != 0:
    updateVF t
  else:
    quit("No updates required.")

proc updateVFGa2*(ts,gs:openarray[float]) =
  record.add(("VFGa(combined)", ts.sum+gs.sum))
  let
    t = ts[0]
    g = gs[0]
  if g != 0:
    if t != 0:
      # Approximate the force gradient update with two Taylor expansions.
      let (tf,tg) = approximateFGcoeff2(t,g)
      xsave.save s.x
      updateVfgF(tg[0])
      updateVF(tf[0])
      s.x.load xsave
      updateVfgF(tg[1])
      updateVF(tf[1])
      s.x.load xsave
    else:
      quit("Force gradient without the force update.")
  elif t != 0:
    updateVF t
  else:
    quit("No updates required.")
