#[
Copyright (c) 2018-2020 Xiao-Yong Jin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
]#

import mdevolve/seqset
import macros, math, strutils

const mdevolveVersion*:int = 1

type
  Integrator* = ref object
    id:string
    components:seqset[Updater]
    # interruption
    doInterrupt:bool
    interrupted:bool
    # actual integrators
    case isSolo:bool
    of false:
      md:seq[IntegratorDt]
      nstep:int
      order:int
    of true:
      need:seqset[Updater]
      step:proc(t:float)  ## update to coeff*dt for exp(tV) or coeff*dt^3 for exp(s[V,[T,V]]).
      # for Force Gradient, exp( t V + s [V,[T,V]] )
      isForceGradient:bool  ## Self is exp(s[V,[T,V]])
      case hasVTV:bool
      of false:
        discard
      of true:
        integrateVTV:Integrator ## This points to exp(s[V,[T,V]]), for self is V.
    # for integrator state
    t:float
    currentStep:int
    currentMD:int
  ParIntegrator* = ref object
    nonZeroStep*:float  ## Consider step size nonzero only if abs(step)>nonZeroStep.
    scale:seq[float]  ## absolute scale of the shared Integrator step for each in the list
    nstep:int
    list:seq[Integrator]
    shared:seqset[Integrator]
    components:seqset[Integrator]
    nonZeroStepWarn*:float  ## Warning of considered nonzero step as zero if abs(step)>nonZeroStepWarn.
  Integrators* = ref object
    ## An object that when indexed returns an Integrator
    id:string
    updater: Updater
    need:seqset[Updater]
    list:seq[Integrator]
  IntegratorDt = object
    I:Integrator
    f:float
  Updater = ref object
    case isSingle:bool
    of false:
      dts: seq[float]
      case hasVTV:bool
      of false:
        runs:proc(ts:openarray[float])
      of true:
        dtgs: seq[float]
        runfgs:proc(ts,ss:openarray[float])
    of true:
      dt: float
      run:proc(t:float)

func indent2(s: string, padding = "  "): string =
  ## Indent `count` spaces for each line excluding the first line.
  ## The output always has a newline character at the end.
  result = ""
  var s = s
  if s[^1] == '\n': s.setlen(s.len-1)
  var first = true
  for line in s.splitLines():
    if first: first = false
    else: result.add padding
    result.add line
    result.add "\n"

func `$`*(x:Integrator):string =
  if x.isSolo:
    result = x.id
  else:
    template mdlist(M:untyped):string =
      var r = ""
      for i in 0..<x.M.len:
        var s = ""
        for m in 0..<i:
          if x.M[i].I == x.M[m].I:
            s = x.M[m].I.id & "{" & $m & "}"
            break
        if s.len == 0: s = $x.M[i].I
        r &= $i & ": " & $x.M[i].f & " * " & indent2(s)
      #if r.len > 0: r.setlen(r.len-1)
      r
    let s = mdlist(md)
    result = ""
    if s.len > 0: result &= "(steps=" & $x.nstep & "){\n  " & indent2(s) & "}\n"
    result = x.id & ":dt^" & $x.order & "{\n  " & indent2(result) & "}"

func `$`*(x:ParIntegrator):string =
  result = "ParIntegrator(steps: " & $x.nstep
  result &= ", nonZeroStep: " & $x.nonZeroStep
  result &= ", nonZeroStepWarn: " & $x.nonZeroStepWarn
  result &= ", scale: " & $x.scale & ", list: @[\n"
  for c in x.list:
    result &= indent2("  " & $c & ",")
  result.setlen(result.len-2)
  result &= "],\n"
  result &= "  shared: " & $x.shared
  result &= ")"

proc steps*(x:Integrator):int =
  if x.isSolo:
    quit("mdevolve: steps: cannot get steps for the integrator: " & $x)
  else:
    return x.nstep
proc steps*(x:ParIntegrator):int =
  x.nstep

proc `steps=`*(x:Integrator, n:int) =
  if x.isSolo:
    quit("mdevolve: `steps=`: cannot set steps for the integrator: " & $x)
  else:
    x.nstep = n

proc `steps=`*(x:ParIntegrator, n:int) =
  x.nstep = n

proc forceGradient*(x:Integrator):Integrator =
  if x.isSolo and x.hasVTV:
    return x.integrateVTV
  else:
    quit("mdevolve: found no force gradient in integrator: " & $x)

proc `forceGradient=`*(x:Integrator, G:Integrator) =
  ## Designate `G` as a force gradient integrator, and associate it with x.
  if x.isSolo:
    G.isForceGradient = true
    let y = Integrator(id:x.id, components:x.components, isSolo:true,
      need:x.need, step:x.step, hasVTV:true, integrateVTV:G)
    x[] = y[]
  else:
    quit("mdevolve: cannot assign force gradient integrator to: " & $x)

proc build(xs:Integrators, n:int):Integrator =
  let u = xs.updater
  if u.dts.len <= n: u.dts.setlen(n+1)
  proc step(t:float) = u.dts[n] += t
  proc gstep(t:float) = u.dtgs[n] += t
  if u.hasVTV:
    if u.dtgs.len <= n: u.dtgs.setlen(n+1)
    let g = Integrator(id:xs.id & "[" & $n & "].VTV", components:newseqset(u),
      isSolo:true, need:xs.need, step:gstep, isForceGradient:true)
    Integrator(id:xs.id & "[" & $n & "]", components:newseqset(u),
      isSolo:true, need:xs.need, step:step, hasVTV:true, integrateVTV:g)
  else:
    Integrator(id:xs.id & "[" & $n & "]", components:newseqset(u),
      isSolo:true, need:xs.need, step:step)

proc `[]`*(xs:Integrators, n:int):Integrator =
  let lo = xs.list.len
  if lo <= n:
    for i in lo..n:
      xs.list.add xs.build n
  xs.list[n]

template wrap(f:proc(t:float)):(Updater, Integrator) =
  let fu = Updater(isSingle:true, dt:0, run:f)
  proc fs(t:float) {.gensym.} = fu.dt += t
  let fi = Integrator(id:asttostr(f), components:newseqset(fu), isSolo:true, step:fs)
  (fu, fi)

template wraps(f:proc(ts:openarray[float])):(Updater, Integrators) =
  let
    fu = Updater(isSingle:false, dts: @[], hasVTV:false, runs:f)
    fi = Integrators(id:asttostr(f), updater:fu)
  (fu, fi)

template wrapss(f:proc(ts,ss:openarray[float])):(Updater, Integrators) =
  let
    fu = Updater(isSingle:false, dts: @[], hasVTV:true, dtgs: @[], runfgs:f)
    fi = Integrators(id:asttostr(f), updater:fu)
  (fu, fi)

template newIntegratorPair*(f,g:proc(t:float)):(Integrator,Integrator) =
  ## Create a pair of integrators by wrapping update functions,
  ## `f` and `g`, inside `Integrator`s.
  ## `f` and `g` are assumed linear, and non-commuting, namely:
  ## - f(t0) f(t1) = f(t0+t1)
  ## - g(t0) g(t1) = g(t0+t1)
  ## - f(tf) g(tg) ≠ g(tg) f(tf)
  ## Such that the new integrators representing two non-commuting linear operators,
  ## f(t0) -> exp(t0 * F)
  ## g(t1) -> exp(t1 * G)
  let
    (fu, fi) = wrap(f)
    (gu, gi) = wrap(g)
  fi.need = newseqset(gu)
  gi.need = newseqset(fu)
  (fi, gi)

template newIntegratorPair*(f:proc(ts:openarray[float]),g:proc(t:float)):(Integrators,Integrator) =
  ## Similar to `newIntegratorPair(f,g:proc(t:float))`,
  ## except that `f` is a combined updater that has
  ## a type of `proc(ts:openarray[float])`.
  let
    (fu, fi) = wraps(f)
    (gu, gi) = wrap(g)
  fi.need = newseqset(gu)
  gi.need = newseqset(fu)
  (fi, gi)
template newIntegratorPair*(f:proc(t:float),g:proc(ts:openarray[float])):(Integrator,Integrators) =
  ## Similar to `newIntegratorPair*(f:proc(ts:openarray[float]),g:proc(t:float))`,
  ## with a reversed order of arguments.
  let (gi, fi) = newIntegratorPair(g, f)
  (fi, gi)

template newIntegratorPair*(f:proc(ts:openarray[float]),g:proc(ts:openarray[float])):(Integrators,Integrators) =
  ## Similar to `newIntegratorPair(f,g:proc(t:float))`,
  ## with both `f` and `g` combined updaters each has
  ## a type of `proc(ts:openarray[float])`.
  let
    (fu, fi) = wraps(f)
    (gu, gi) = wraps(g)
  fi.need = newseqset(gu)
  gi.need = newseqset(fu)
  (fi, gi)

template newIntegratorPair*(f:proc(ts,ss:openarray[float]),g:proc(t:float)):(Integrators,Integrator) =
  ## Similar to `newIntegratorPair*(f:proc(ts:openarray[float]),g:proc(t:float))`,
  ## except that `f` is a combined updater that has
  ## a type of `proc(ts,ss:openarray[float])`,
  ## accepting second time step for the force gradient term,
  ## with t to the first and s to the second argument evaluating exp(t*V + s*[V,[T,V]]).
  let
    (fu, fi) = wrapss(f)
    (gu, gi) = wrap(g)
  fi.need = newseqset(gu)
  gi.need = newseqset(fu)
  (fi, gi)
template newIntegratorPair*(f:proc(t:float),g:proc(ts,ss:openarray[float])):(Integrator,Integrators) =
  ## Similar to `newIntegratorPair*(f:proc(ts,ss:openarray[float]),g:proc(t:float))`,
  ## with a reversed order of arguments.
  let (gi, fi) = newIntegratorPair(g, f)
  (fi, gi)

proc createSingle(f:NimNode):(NimNode, NimNode, NimNode) =
  let
    fu = gensym(nskLet, $f & ":U")
    fi = gensym(nskLet, $f & ":I")
    body = quote do:
      let (`fu`, `fi`) = wrap(`f`)
  (body, fu, fi)

proc createMulti(fs:NimNode):(NimNode, NimNode, NimNode) =
  let
    body = newStmtList()
    fu = newPar()
    fi = newPar()
  for f in fs:
    let (b,u,i) = createSingle f
    body.add b
    fu.add u
    fi.add i
  (body, fu, fi)

proc createIntegrators(fs:NimNode):(NimNode, NimNode, NimNode) =
  let tfs = fs.gettype
  # echo tfs.treerepr
  proc isProcFloat(n:NimNode):bool =
    # echo n.treerepr
    n.kind == nnkBracketExpr and n.len == 3 and
      n[0].eqident"proc" and n[1].eqident"void" and n[2].eqident"float"
  if tfs.kind == nnkBracketExpr:
    if tfs.isProcFloat:
      return createSingle(fs)
    elif tfs[0].eqident"tuple":
      var allProcFloat = true
      for i in 1..<tfs.len:
        allProcFloat = allProcFloat and tfs[i].isProcFloat
      if allProcFloat:
        return createMulti(fs)
  error("newIntegratorPair: requires a single `proc(t:float)` or a tuple of those; got\n" &
    fs.repr & "\nof type:\n" & tfs.treerepr)

proc newIntegratorPairTuple(fs,gs:NimNode):NimNode =
  let
    (body0, fus, fis) = createIntegrators(fs)
    (body1, gus, gis) = createIntegrators(gs)
    fneed = newcall(bindsym"newseqset")
    gneed = newcall(bindsym"newseqset")
  if gus.kind == nnkPar:
    for gu in gus:
      fneed.add gu
  else:
    fneed.add gus
  if fus.kind == nnkPar:
    for fu in fus:
      gneed.add fu
  else:
    gneed.add fus
  result = newStmtList().add(body0, body1)
  if fis.kind == nnkPar:
    for fi in fis:
      result.add newAssignment(newDotExpr(fi,ident"need"), fneed)
  else:
    result.add newAssignment(newDotExpr(fis,ident"need"), fneed)
  if gis.kind == nnkPar:
    for gi in gis:
      result.add newAssignment(newDotExpr(gi,ident"need"), gneed)
  else:
    result.add newAssignment(newDotExpr(gis,ident"need"), gneed)
  result.add newPar(fis, gis)

macro newIntegratorPair*(fs:tuple, g:proc(t:float)):auto =
  ## Similar to templates,
  ## `newIntegratorPair(f:proc(ts:openarray[float]),g:proc(t:float))`,
  ## and relatives, with `f` being a tuple of a list of `proc(t:float)`
  ## instead of a pre-packaged `proc(ts:openarray[float])`.
  ## Note that for force gradient update, it must be associated with the
  ## force update explicitly with `forceGradient=`*(x:Integrator, G:Integrator).
  newIntegratorPairTuple(fs,g)

macro newIntegratorPair*(f:proc(t:float), gs:tuple):auto =
  ## Similar to the macro,
  ## `newIntegratorPair*(fs:tuple, g:proc(t:float))`,
  ## but with a reversed order of arguments.
  newIntegratorPairTuple(f,gs)

macro newIntegratorPair*(fs:tuple, gs:tuple):auto =
  ## Similar to the macro,
  ## `newIntegratorPair*(fs:tuple, g:proc(t:float))`,
  ## but both arguments are tuples of a list of `proc(t:float)`.
  newIntegratorPairTuple(fs,gs)

func integratorDt(I:Integrator, f:float):IntegratorDt =
  IntegratorDt(I: I, f: f)

proc newSerialEvolution*(id:string, mds:seq[tuple[i:Integrator,f:float]], steps:int, order:int):Integrator =
  var ms = newseq[IntegratorDt]()
  var c = newseqset[Updater]()
  for m in mds:
    ms.add integratorDt(m.i, m.f)
    c.add m.i.components
  Integrator(id: id, components: c, isSolo: false, md: ms, nstep: steps, order: order)

proc absScale(x:Integrator, s:Integrator):float =
  if x == s: return 1
  elif not x.isSolo:
    var f = 0.0
    for c in x.md:
      f += abs(c.f) * c.I.absScale(s)
    return f
  else: return 0

proc add*(r:ParIntegrator, xs:varargs[Integrator]) =
  ## Add a list of integrators to the ParIntegrator.
  ## There must be one and only one shared sub-integrator.
  var mds = newseqset[Integrator]()
  proc collect(x:Integrator) =
    mds.add x
    if not x.isSolo:
      for c in x.md: c.I.collect
  for x in xs:
    mds.empty
    x.collect
    r.shared.add r.components.addReturnIntersection mds
    r.list.add x
  if r.list.len > 1 and r.shared.len != 1:
    # We can add support for multiple shared integrators if the need arises.
    quit("Error: parallel evolution has " & $r.shared.len & " shared integrators:\n" & $r.shared)
  if r.shared.len == 1:
    let s = r.shared[0]
    s.doInterrupt = true
    for i in r.scale.len..<r.list.len:
      r.scale.add r.list[i].absScale s

proc newParallelEvolution*(xs:varargs[Integrator]):ParIntegrator =
  ## Create an integrator scheme from a list of integrators,
  ## where each integrator runs independently except for the one shared sub-integrator
  ## that gets synchronized to be at time steps specified from individual integrators.
  ## Currently only works with one shared sub-integrator.
  ## The scheme considers any normalized step smaller than `nonZeroStep` (default 1E-12) zero,
  ## and it prints a warning if such "zero" step is larger than `nonZeroStepWarn` (default 1E-15).
  result = ParIntegrator(nonZeroStep: 1E-12, scale: @[], nstep: 1,
    list: @[], shared: newseqset[Integrator](),
    components: newseqset[Integrator](),
    nonZeroStepWarn: 1E-15)
  result.add xs

proc update(x:Updater) =
  if x.isSingle:
    if x.dt != 0:
      x.run x.dt
      x.dt = 0
  else:
    var
      run = false
      rung = false
    for t in x.dts:
      if t != 0:  # any nonzero dt
        run = true
        break
    if x.hasVTV:
      for t in x.dtgs:
        if t != 0:  # any nonzero dt
          rung = true
          break
    if run or rung:
      if x.hasVTV:
        x.runfgs(x.dts, x.dtgs)
      else:
        x.runs(x.dts)
      if run:
        for m in x.dts.mitems: m = 0
      if rung:
        for m in x.dtgs.mitems: m = 0

proc evolve*(x:Integrator, t:float)

proc evolve(x:Integrator, dt:float, coeff:float) =
  if x.isSolo and x.isForceGradient:
    x.evolve(coeff*dt*dt*dt)
  else:
    x.evolve(coeff*dt)

proc evolve(x:Integrator) =
  if x.doInterrupt and not x.interrupted:
    x.interrupted = true
    return
  if x.interrupted:
    x.interrupted = false
  if x.t == 0:
    return
  if x.isSolo:
    # echo x.id,"  ",x.t
    for c in x.need: c.update
    x.step x.t
    return
  # On exit, save current state if interrupted.
  let dt = x.t / x.nstep.float
  var i = x.currentMD
  for s in x.currentStep..<x.nstep:
    while true:  # First loop start with i = x.currentMD, 0 otherwise.
      let md = x.md[i]
      md.I.evolve(dt, md.f)
      if md.I.interrupted:
        x.interrupted = true
        x.currentStep = s
        x.currentMD = i
        return
      inc i
      if i == x.md.len:
        i = 0
        break
  # reset to zero for normal exit
  x.currentStep = 0
  x.currentMD = 0

proc evolve*(x:Integrator, t:float) =
  if x.interrupted:
    # Ignore t, as we might have set x.t outside.
    x.evolve
  else:
    x.t = t
    x.currentStep = 0
    x.currentMD = 0
    x.evolve

proc evolve*(x:ParIntegrator, t:float) =
  # The shared integrator is referenced by all.  We need to take care of its state.
  let
    dt = t / x.nstep.float
    ls = x.list
    sc = x.scale
    si = x.shared[0]  # works only for one shared
  var maxZeroStep = 0.0
  for step in 0..<x.nstep:
    var
      first = true
      ss = newSeq[float](ls.len)  # signed step size
      af = newSeq[float](ls.len)  # absolute integration fraction
      running = newSeq[int](ls.len)
      next = newSeq[int](ls.len)
    for i in 0..<next.len:
      ss[i] = 0
      af[i] = 0
      running[i] = i
      next[i] = i
      ls[i].t = dt
    while true:
      for i in next:
        if first:
          si.interrupted = false
        else:
          si.interrupted = true
          si.t = ss[i]
        ls[i].evolve
        if ls[i].interrupted:
          ss[i] = si.t
          af[i] = abs(ss[i] / sc[i])
      first = false
      var fmin = 1.5 * abs(dt)  # some value with its abs larger than abs(dt)
      var imin = -1
      for k in countdown(running.len-1,0):
        let i = running[k]
        if not ls[i].interrupted:
          # ls[i] has finished.
          running.delete k  # changes the order, only works with countdown
          continue
        let f = af[i]
        if f < fmin - x.nonZeroStep:
          fmin = f
          imin = i
        elif f < fmin + x.nonZeroStep:  # tolerently equal
          # Always pick the smallest signed step.
          # This is IMPORTANT for reversibility.
          let
            ssi = ss[i]
            ssim = ss[imin]
          if ssi < ssim - x.nonZeroStep:
            fmin = f
            imin = i
          elif ssi < ssim + x.nonZeroStep:  # tolerently equal
            # Than we pick according to the order in the integrator list.
            if i < imin:
              fmin = f
              imin = i
      next.setlen 0
      if imin >= 0: next.add imin
      else: break
#[
      proc showss(p:string) =
        var s = p
        for i,f in ss:
          let ss = if i==imin: " #" elif ls[i].interrupted: " *" else: "  "
          s &= $f & ss & ", "
        s.setlen(s.len-2)
        s &= "    abs " & $af
        echo s
      showss("====")
]#
      for i in running:
        if i == imin: continue
        let
          f = ss[i] - ss[imin]
          a = abs(af[i] - fmin)
        if f > -x.nonZeroStep and f < x.nonZeroStep:
          if f != 0:
            let ff = abs f
            if ff > maxZeroStep: maxZeroStep = ff
          ss[i] = 0
        else: ss[i] = f
        if a < x.nonZeroStep:
          if a != 0:
            if a > maxZeroStep: maxZeroStep = a
          af[i] = 0
        else: af[i] = a
        if ss[i] == 0 and af[i] == 0: next.add i
      # showss(" << ")
  if maxZeroStep > x.nonZeroStepWarn:
    echo "Warning: Integrator(",si.id,") ignored small nonzero steps as large as ",maxZeroStep

proc finish*(x:Integrator) =
  ## Finish any dangling evolutions.  As `evolve` has lazy semantics,
  ## we call this procedure to make `evolve` finish their steps.
  for c in x.components: c.update

proc finish*(x:ParIntegrator) =
  for c in x.list: c.finish

proc approximateFGcoeff*(cv,cg:float):tuple[tf:float,tg:float] =
  ## Return the coefficients for approximate force gradient update using
  ## (cv*V + cg*[V,[T,V]])(x, p) -> (x+tf*S'(x + tg*S'(x)), p)
  ## Yin and Mawhinney, (2012)
  let tg = 2.0*cg/cv
  (cv, tg)

proc approximateFGcoeff2*(cv,cg:float, b = 0.2):tuple[tf:array[2,float],tg:array[2,float]] =
  ## Return the coefficients for approximate force gradient update using
  ## a second order approximation,
  ## (cv*V + cg*[V,[T,V]])(x, p) ->
  ##   (x + tf[0]*S'(x + tg[0]*S'(x)) + tf[1]*S'(x + tg[1]*S'(x)), p)
  ## Extending the approximation suggested in Yin and Mawhinney, (2012).
  ## We use two Taylor expansion to cancel the second order term.
  let
    h = 2.0*cg/cv
    a = (b*b+1.0)/(2.0*b)
    tf = [cv*(1.0+a)/2.0, cv*(1.0-a)/2.0]
    tg = [h*(1.0+b)/(1.0+a), h*(1.0-b)/(1.0-a)]
  (tf, tg)

proc mkLeapFrog*(T,V:Integrator, steps = 1):Integrator =
  newSerialEvolution("LeapFrog",
    @[(T, 0.5), (V, 1.0), (T, 0.5)],
    steps, 2)

proc mkOmelyan2MN*(T,V:Integrator, steps = 1,
    lambda = 0.1931833275037836):Integrator =
  ## Omelyan et. al. (2003), equation (31)
  newSerialEvolution("Omelyan2MN",
    @[(T, lambda), (V, 0.5), (T, 1.0-2*lambda), (V, 0.5), (T, lambda)],
    steps, 2)

proc mkSW92*(T,V:Integrator, steps = 1):Integrator =
  ## Sexton & Weingarten (1992), TVTVT
  result = mkOmelyan2MN(T,V,steps,1.0/6.0)
  result.id = "SW92"

proc mkOmelyan4MN4FP*(T,V:Integrator, steps = 1,
    rho = 0.1786178958448091,
    theta = -0.06626458266981843,
    lambda = 0.7123418310626056):Integrator =
  ## 4th order minimum norm, 4 force evaluation, position version.
  ## Omelyan et. al. (2003), equation (58) and (62)
  newSerialEvolution("Omelyan4MN4FP",
    @[(T, rho), (V, lambda), (T, theta), (V, 0.5-lambda), (T, 1.0-2*(theta+rho)),
      (V, 0.5-lambda), (T, theta), (V, lambda), (T, rho)],
    steps, 4)

proc mkOmelyan4MN5FV*(T,V:Integrator, steps = 1,
    rho = 0.2539785108410595,
    theta = -0.03230286765269967,
    vartheta = 0.08398315262876693,
    lambda = 0.6822365335719091):Integrator =
  ## 4th order minimum norm, 5 force evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (63) and (71)
  newSerialEvolution("Omelyan4MN5FV",
    @[(V, vartheta), (T, rho), (V, lambda), (T, theta), (V, 0.5-(lambda+vartheta)), (T, 1.0-2*(theta+rho)),
      (V, 0.5-(lambda+vartheta)), (T, theta), (V, lambda), (T, rho), (V, vartheta)],
    steps, 4)

proc mkOmelyan4MN5FP*(T,V:Integrator, steps = 1,
    rho = 0.2750081212332419,
    theta = -0.1347950099106792,
    vartheta= -0.08442961950707149,
    lambda = 0.3549000571574260):Integrator =
  ## 4th order minimum norm, 5 force evaluation, position version.
  ## Omelyan et. al. (2003), equation (72) and (80)
  newSerialEvolution("Omelyan4MN5FP",
    @[(T, rho), (V, vartheta), (T, theta), (V, lambda),
      (T, 0.5-(theta+rho)),
      (V, 1.0-2.0*(lambda+vartheta)),
      (T, 0.5-(theta+rho)),
      (V, lambda), (T, theta), (V, vartheta), (T, rho)],
    steps, 4)

proc mkOmelyan6MN7FV*(T,V:Integrator, steps = 1):Integrator =
  ## 6th order minimum norm, 7 force evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (83)
  let
    a = [ 0.2465881872786138,
          0.6047073875057809,
         -0.4009869039788007,
          0.09938265838881204]
    b = [ 0.08333333333333333,
          0.3977675859548440,
         -0.03933369314462574,
          0.05823277385644840]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..<b.len:
    mds.add (i:V, f:b[i])
    mds.add (i:T, f:a[i])
  mds.add (i:V, f:b[^1])
  for i in countdown(b.len-2,0):
    mds.add (i:T, f:a[i])
    mds.add (i:V, f:b[i])
  newSerialEvolution("Omelyan6MN7FV", mds, steps, 6)

proc mkOmelyan4MN2F2GV*(T,V:Integrator, steps = 1,
    xi = -17.0/18000.0):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 2 force evaluation, 2 gradient evaluation, velocity version.
  ## Omelyan et. al. (2002), equation (30) and (35).
  ## Omelyan et. al. (2003), equation (28).
  let
    G = V.forceGradient
    lambda = 1.0/6.0
    chi = 1.0/12.0 - 2.0*xi - 0.5*lambda*(1.0-lambda)
  newSerialEvolution("Omelyan4MN2F2GV",
    @[(V, lambda), (G, xi), (T, 0.5), (V, 1.0-2.0*lambda), (G, chi), (T, 0.5), (G, xi), (V, lambda)],
    steps, 4)

proc mkOmelyan4MN2F1GV*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 2 force evaluation, 1 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (29).
  ## A special case of Omelyan4MN2F2GV with only 1 gradient evaluation.
  let
    G = V.forceGradient
    lambda = 1.0/6.0
    chi = 1.0/12.0 - 0.5*lambda*(1.0-lambda)
  newSerialEvolution("Omelyan4MN2F1GV",
    @[(V, lambda), (T, 0.5), (V, 1.0-2.0*lambda), (G, chi), (T, 0.5), (V, lambda)],
    steps, 4)

proc mkOmelyan4MN3F3GP*(T,V:Integrator, steps = 1,
    lambda = 0.2825633404177051,
    chi = 0.003035236056708454):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 3 force evaluation, 3 gradient evaluation, position version.
  ## Omelyan et. al. (2003), equation (40), (41), and (42).
  let
    G = V.forceGradient
    theta = 0.5 - 1.0 / sqrt(24.0 * lambda)
    xi = (1.0 - 12.0*chi - sqrt(6.0 * lambda) * (1.0 - lambda)) / 24.0
  newSerialEvolution("Omelyan4MN3F3GP",
    @[(T, theta), (V, lambda), (G, xi), (T, 0.5-theta),
      (V, 1.0-2.0*lambda), (G, chi),
      (T, 0.5-theta), (G, xi), (V, lambda), (T, theta)],
    steps, 4)

proc mkOmelyan4MN3F2GP*(T,V:Integrator, steps = 1,
    lambda = 0.3152315246820299):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 3 force evaluation, 2 gradient evaluation, position version.
  ## Omelyan et. al. (2003), equation (43).
  let
    G = V.forceGradient
    theta = 0.5 - 1.0 / sqrt(24.0 * lambda)
    xi = (1.0 - sqrt(6.0 * lambda) * (1.0 - lambda)) / 24.0
  newSerialEvolution("Omelyan4MN3F2GP",
    @[(T, theta), (V, lambda), (G, xi), (T, 0.5-theta),
      (V, 1.0-2.0*lambda),
      (T, 0.5-theta), (G, xi), (V, lambda), (T, theta)],
    steps, 4)

proc mkOmelyan4MN3F1GP*(T,V:Integrator, steps = 1,
    lambda = 0.2470939580390842):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 3 force evaluation, 1 gradient evaluation, position version.
  ## Omelyan et. al. (2002), equation (38), (40), and (42).
  ## Omelyan et. al. (2003), equation (44)
  let
    G = V.forceGradient
    theta = 0.5 - 1.0 / sqrt(24.0 * lambda)
    chi = (1.0 - sqrt(6.0 * lambda) * (1.0 - lambda)) / 12.0
  newSerialEvolution("Omelyan4MN3F1GP",
    @[(T, theta), (V, lambda), (T, 0.5-theta),
      (V, 1.0-2.0*lambda), (G, chi),
      (T, 0.5-theta), (V, lambda), (T, theta)],
    steps, 4)

proc mkOmelyan4MN3F3GV*(T,V:Integrator, steps = 1,
    theta = 0.2728983001988755,
    chi = 0.002960781208329478):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 3 force evaluation, 3 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (34), (35), and (36).
  let
    G = V.forceGradient
    lambda = (6.0 + 1.0 / (theta*(theta-1.0))) / 12.0
    t1 = theta-1.0
    xi = - (6.0 + 288.0*chi - 1.0 / (theta*t1*t1)) / 288.0
  newSerialEvolution("Omelyan4MN3F3GV",
    @[(V, lambda), (G, xi), (T, theta), (V, 0.5-lambda), (G, chi),
      (T, 1-2.0*theta),
      (G, chi), (V, 0.5-lambda), (T, theta), (G, xi), (V, lambda)],
    steps, 4)

proc mkOmelyan4MN3F2GV*(T,V:Integrator, steps = 1,
    theta = 0.2813980611667719):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 3 force evaluation, 2 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (37)
  let
    G = V.forceGradient
    lambda = (6.0 + 1.0 / (theta*(theta-1.0))) / 12.0
    t1 = theta-1.0
    chi = (1.0 - 6.0*theta*(1.0-theta*(2.0-theta)))/(288.0*t1*t1*theta)
  newSerialEvolution("Omelyan4MN3F2GV",
    @[(V, lambda), (T, theta), (V, 0.5-lambda), (G, chi),
      (T, 1-2.0*theta),
      (G, chi), (V, 0.5-lambda), (T, theta), (V, lambda)],
    steps, 4)

proc mkOmelyan4MN3F1GV*(T,V:Integrator, steps = 1,
    theta = 0.2409202729169543):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 3 force evaluation, 1 gradient evaluation, velocity version.
  ## Omelyan et. al. (2002), equation (39) and (44).
  ## Omelyan et. al. (2003), equation (38)
  result = mkOmelyan4MN3F3GV(T,V, steps, theta = theta, chi = 0)
  result.id = "Omelyan4MN3F1GV"

proc mkOmelyan4MN4F2GVG*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 4 force evaluation, 1 gradient evaluation, velocity (gradient) version.
  ## Omelyan et. al. (2003), equation (48) and (53).
  let
    G = V.forceGradient
    theta = 0.1921125277429464
    vartheta = 0.05851872613455621
    lambda = 0.2852162240687091
    chi = 0.002427475259663050
    mu = 0.0004339598806816256
  newSerialEvolution("Omelyan4MN4F2GV",
    @[(V, vartheta), (G, mu), (T, theta), (V, lambda),
      (T, 0.5-theta),
      (V, 1.0-2.0*(lambda+vartheta)), (G, chi),
      (T, 0.5-theta),
      (V, lambda), (T, theta), (G, mu), (V, vartheta)],
    steps, 4)

proc mkOmelyan4MN4F2GV*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 4 force evaluation, 1 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (48) and (54).
  ## Cf. mkOmelyan4MN4F2GVG
  let
    G = V.forceGradient
    theta = 0.2189286596427438
    vartheta = 0.06840805970727767
    lambda = 0.3109406355938166
    xi = 0.001602503681334363
  newSerialEvolution("Omelyan4MN4F2GV",
    @[(V, vartheta), (T, theta), (V, lambda), (G, xi),
      (T, 0.5-theta),
      (V, 1.0-2.0*(lambda+vartheta)),
      (T, 0.5-theta),
      (G, xi), (V, lambda), (T, theta), (V, vartheta)],
    steps, 4)

proc mkOmelyan4MN5F2GV*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 5 force evaluation, 2 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (63) and (68)
  let
    G = V.forceGradient
    rho = 0.2029270564692829
    theta = 0.1926052063353027
    vartheta = 0.06668876199434440
    lambda = 0.2620356629687677
    xi = 0.001042387551227681
  newSerialEvolution("Omelyan4MN5F2GV",
    @[(V, vartheta), (T, rho), (V, lambda), (G, xi),
      (T, theta), (V, 0.5-(lambda+vartheta)),
      (T, 1.0-2*(theta+rho)),
      (V, 0.5-(lambda+vartheta)), (T, theta),
      (G, xi), (V, lambda), (T, rho), (V, vartheta)],
    steps, 4)

proc mkOmelyan4MN5F1GV*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 5 force evaluation, 1 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (63) and (70)
  let
    G = V.forceGradient
    rho = 0.2797644436188271
    theta = -0.001180329820696323
    vartheta = 0.08010998355755116
    lambda = -2.0220148671481104
    mu = 0.0003098750751031143
  newSerialEvolution("Omelyan4MN5F1GV",
    @[(V, vartheta), (G, mu), (T, rho), (V, lambda),
      (T, theta), (V, 0.5-(lambda+vartheta)),
      (T, 1.0-2*(theta+rho)),
      (V, 0.5-(lambda+vartheta)), (T, theta),
      (V, lambda), (T, rho), (G, mu), (V, vartheta)],
    steps, 4)

proc mkOmelyan4MN5F2GP*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 5 force evaluation, 2 gradient evaluation, position version.
  ## Omelyan et. al. (2003), equation (72) and (78)
  let
    G = V.forceGradient
    rho = 0.06419108866816235
    theta = 0.1919807940455741
    vartheta = 0.1518179640276466
    lambda = 0.2158369476787619
    xi = 0.0009628905212024874
  newSerialEvolution("Omelyan4MN5F2GV",
    @[(T, rho), (V, vartheta), (T, theta), (V, lambda), (G, xi),
      (T, 0.5-(theta+rho)),
      (V, 1.0-2.0*(lambda+vartheta)),
      (T, 0.5-(theta+rho)),
      (G, xi), (V, lambda), (T, theta), (V, vartheta), (T, rho)],
    steps, 4)

proc mkOmelyan4MN5F1GP*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 4.
  ## Minimum norm, 5 force evaluation, 2 gradient evaluation, position version.
  ## Omelyan et. al. (2003), equation (72) and (79)
  let
    G = V.forceGradient
    rho = 0.1255768596433302
    theta = -0.002407093745014925
    vartheta = -0.8938074259467744
    lambda = 1.1758501877269955
    chi = 0.002952744354631969
  newSerialEvolution("Omelyan4MN5F2GV",
    @[(T, rho), (V, vartheta), (T, theta), (V, lambda),
      (T, 0.5-(theta+rho)),
      (V, 1.0-2.0*(lambda+vartheta)), (G, chi),
      (T, 0.5-(theta+rho)),
      (V, lambda), (T, theta), (V, vartheta), (T, rho)],
    steps, 4)

proc mkOmelyan6MN5F5GV*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 6.
  ## Minimum norm, 5 force evaluation, 5 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (63) and (64)
  let
    G = V.forceGradient
    rho = 0.5309910490348568
    theta = -0.2573883543804353
    vartheta = 0.08281492492827128
    lambda = 0.008354543940755644
    xi = -0.0002401600937577623
    chi = 0.004267631995107088
    mu = -0.0001633190022736910
  newSerialEvolution("Omelyan4MN5F5GV",
    @[(V, vartheta), (G, mu), (T, rho), (V, lambda), (G, xi),
      (T, theta), (V, 0.5-(lambda+vartheta)), (G, chi),
      (T, 1.0-2*(theta+rho)),
      (G, chi), (V, 0.5-(lambda+vartheta)), (T, theta),
      (G, xi), (V, lambda), (T, rho), (G, mu), (V, vartheta)],
    steps, 6)

proc mkOmelyan6MN5F4GV*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 6.
  ## Minimum norm, 5 force evaluation, 4 gradient evaluation, velocity version.
  ## Omelyan et. al. (2003), equation (63) and (64)
  let
    G = V.forceGradient
    s = sqrt(5.0)
    t = sqrt(50.0+22.0*s)
    rho = (1.0 + 1.0/s)/2.0
    theta = -1.0/s
    vartheta = 1.0/12.0
    lambda = 5.0/12.0 - t/24.0
    xi = (15.0+5.0*s)/1152.0 - t*(1.0/2880.0 + s/1152.0)
    chi = -(11.0+5.0*s)/1152.0 + t*(1.0/2880.0 + s/1152.0)
  newSerialEvolution("Omelyan4MN5F4GV",
    @[(V, vartheta), (T, rho), (V, lambda), (G, xi),
      (T, theta), (V, 0.5-(lambda+vartheta)), (G, chi),
      (T, 1.0-2*(theta+rho)),
      (G, chi), (V, 0.5-(lambda+vartheta)), (T, theta),
      (G, xi), (V, lambda), (T, rho), (V, vartheta)],
    steps, 6)

proc mkOmelyan6MN5F5GP*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 6.
  ## Omelyan et. al. (2003), equation (72) and (73).
  let
    G = V.forceGradient
    rho = 0.1098059301577147
    theta = 0.4828099940251012
    vartheta = 0.2693816517677854
    lambda = 0.07611936345860829
    xi = -0.001803378129376054
    chi = 0.01083650107661986
    mu = 0.001011249349033012
  newSerialEvolution("Omelyan6MN5F3GV",
    @[(T, rho), (V, vartheta), (G, mu), (T, theta), (V, lambda), (G, xi),
      (T, 0.5-(theta+rho)),
      (V, 1.0-2.0*(lambda+vartheta)), (G, chi),
      (T, 0.5-(theta+rho)),
      (G, xi), (V, lambda), (T, theta), (G, mu), (V, vartheta), (T, rho)],
    steps, 6)

proc mkOmelyan6MN5F3GP*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 6.
  ## Omelyan et. al. (2002), equation (50) and (51).
  ## Omelyan et. al. (2003), equation (72) and (76).
  let
    G = V.forceGradient
    rho = 0.1097059723948682
    theta = 0.4140632267310831
    vartheta = 0.2693315848935301
    lambda = 1.1319803486515564
    chi = -0.01324638643416052
    mu = 0.0008642161339706166
  newSerialEvolution("Omelyan6MN5F3GV",
    @[(T, rho), (V, vartheta), (G, mu), (T, theta), (V, lambda), (T, 0.5-(theta+rho)),
      (V, 1.0-2.0*(lambda+vartheta)), (G, chi),
      (T, 0.5-(theta+rho)), (V, lambda), (T, theta), (G, mu), (V, vartheta), (T, rho)],
    steps, 6)

proc mkOmelyan8S11F11GP*(T,V:Integrator, steps = 1):Integrator =
  ## Force Gradient Integrator of order 8, Omelyan et. al. (2002), equation (53)
  let
    G = V.forceGradient
    a = [ 4.1009674738801111928784693005080E-1,
         -3.4123345756052780489101697378499E-1,
          2.5644714021068150492361761631743E-1,
          2.7765273975812438394100476242641E-1,
         -5.6926266869753773902939657321159E-1,
          4.6629949890124853576794423820194E-1]
    b = [ 4.8249309817414952912695842664785E-3,
          1.7492394861090375603419001374207E-1,
          2.9304366370957066164364546204288E-1,
          4.7448940168459770284238136482511E-2,
         -1.5299863411743974499219652320477E-3,
         -3.7422994259002571606842462603791E-2]
    c = [ 1.4743936907797528364717244760736E-4,
          2.3288450531932545357194967600155E-4,
          6.1648659635535962497705619884752E-3,
         -1.2307516860831240716732016960034E-2,
         -7.3296648559126385387017161643798E-5,
          1.5295860994523744731993293847001E-2]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..5:
    mds.add (i:T, f:a[i])
    mds.add (i:V, f:b[i])
    mds.add (i:G, f:c[i])
  mds.add (i:T, f:a[5])
  for i in countdown(4,0):
    mds.add (i:G, f:c[i])
    mds.add (i:V, f:b[i])
    mds.add (i:T, f:a[i])
  newSerialEvolution("Omelyan8S11F11GP", mds, steps, 8)

proc mkCreutzGocksch*(S:Integrator, steps = 1):Integrator =
  ## A Triplet concatenation that construct an integrator of order K + 2,
  ## from any self-adjoint integrotor of order K.
  ## Creutz and Gocksch, (1989)
  let
    k = S.order
    dk = 1.0 / (2.0 - pow(2.0, (1.0/float(k+1))))
  newSerialEvolution("CreutzGocksch" & $k, @[(S, dk), (S, 1.0-2.0*dk), (S, dk)], steps, k+2)

proc mkOmelyan8CK4P4*(S4:Integrator, steps = 1):Integrator =
  ## Omelyan 8th order composition scheme (Q = 8), using a 4th order integrator (K = 4), `S4`,
  ## with 7 evaluations of the integrator (2P-1, for P = 4).
  ## Omelyan et. al. (2002), equation (61)
  let
    d = [ 0.8461211474696757,
          0.1580128458008567,
         -1.090206660543938,
          1.172145334546811]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..<d.len:
    mds.add (i:S4, f:d[i])
  for i in countdown(d.len-2,0):
    mds.add (i:S4, f:d[i])
  newSerialEvolution("Omelyan8CK4P4", mds, steps, 8)

proc mkOmelyan10CK4P7*(S4:Integrator, steps = 1):Integrator =
  ## Omelyan 10th order composition scheme (Q = 10), using a 4th order integrator (K = 4), `S4`,
  ## with 13 evaluations of the integrator (2P-1, for P = 7).
  ## Omelyan et. al. (2002), equation (62)
  let
    d = [ 0.80523995769578082326628169802782,
         -0.49193105914623101022388138864143,
          0.35449258654398460535529269988483,
         -0.69573922271140223803036463461997,
          0.39959538030329256359349977087819,
          0.54979568601438452794128031563760]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..<d.len:
    mds.add (i:S4, f:d[i])
  mds.add (i:S4, f:1.0-2.0*d.sum)
  for i in countdown(d.len-1,0):
    mds.add (i:S4, f:d[i])
  newSerialEvolution("Omelyan10CK4P7", mds, steps, 10)

proc mkOmelyan12CK4P12*(S4:Integrator, steps = 1):Integrator =
  ## Omelyan 12th order composition scheme (Q = 12), using a 4th order integrator (K = 4), `S4`,
  ## with 23 evaluations of the integrator (2P-1, for P = 12).
  ## Omelyan et. al. (2002), equation (62)
  let
    d = [ 0.17385016093097855436061712858303,
          0.53377479890712207949282653990842,
          0.12130138614668307673802291966495,
          0.29650747033807195273440032505629,
         -0.59965999857335454018482312008233,
          0.09043581286204437145871130429094,
         -0.43979146257635806886778748138962,
         -0.30251552922346495057010240779104,
          0.59895872989247982114545906953712,
          0.31236416538275576151816280776696,
         -0.59081230769647833184090443445303]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..<d.len:
    mds.add (i:S4, f:d[i])
  mds.add (i:S4, f:1.0-2.0*d.sum)
  for i in countdown(d.len-1,0):
    mds.add (i:S4, f:d[i])
  newSerialEvolution("Omelyan12CK4P12", mds, steps, 12)

proc mkOmelyan10CK6P4*(S6:Integrator, steps = 1):Integrator =
  ## Omelyan 10th order composition scheme (Q = 10), using a 6th order integrator (K = 6), `S6`,
  ## with 7 evaluations of the integrator (2P-1, for P = 4).
  ## Omelyan et. al. (2002), equation (68)
  let
    d = [ 0.88480139304442862590773863625720,
          0.11922404430206648052593264029266,
         -1.0677277516805770678518370004925]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..<d.len:
    mds.add (i:S6, f:d[i])
  mds.add (i:S6, f:1.0-2.0*d.sum)
  for i in countdown(d.len-1,0):
    mds.add (i:S6, f:d[i])
  newSerialEvolution("Omelyan10CK6P4", mds, steps, 10)

proc mkOmelyan12CK6P7*(S6:Integrator, steps = 1):Integrator =
  ## Omelyan 12th order composition scheme (Q = 12), using a 6th order integrator (K = 6), `S6`,
  ## with 13 evaluations of the integrator (2P-1, for P = 7).
  ## Omelyan et. al. (2002), equation (70)
  let
    d = [ 0.64725339206305240605385248392083,
          0.44631941526959576960102601257986,
         -0.66447133641046221008529452937721,
         -0.58260619571844248816548809046510,
          0.64081619589013117205634311707157,
          0.31805596598883340430918587031701]
  var mds = newseq[tuple[i:Integrator,f:float]]()
  for i in 0..<d.len:
    mds.add (i:S6, f:d[i])
  mds.add (i:S6, f:1.0-2.0*d.sum)
  for i in countdown(d.len-1,0):
    mds.add (i:S6, f:d[i])
  newSerialEvolution("Omelyan12CK6P7", mds, steps, 12)
