#[
Copyright (c) 2018, 2019 Xiao-Yong Jin

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


import macros

const mdevolveVersion*:int = 0

template debug(xs:varargs[untyped]) =
  when false:
    echo xs

proc `@`(f,x:NimNode):NimNode =
  if f.kind in {nnkCall,nnkCommand}:
    result = f.copy.add x
  elif f.kind in {nnkIdent,nnkSym,nnkDotExpr}:
    result = f.newCall.add x
  elif f.kind == nnkStmtList:
    result = f.copyNimNode
    for i in 0..<f.len-1: result.add f[i].copy
    result.add f[^1]@x
  else:
    var m = "Calling:\n" & f.repr
    m &= "\nCannot call with arguments: " & x.repr
    error m
  #echo "@\n",result.repr
macro `@:`(f,x:untyped):untyped = f@x

type
  MDEvolveS1 = object
    cs:array[2,float]
    nstep:int
  MDEvolveS2 = object
    cs:array[3,float]
    nstep:int
  MDEvolveS3 = object
    cs:array[4,float]
    nstep:int
  MDEvolveS4 = object
    cs:array[5,float]
    nstep:int
  MDEvolveS5 = object
    cs:array[6,float]
    nstep:int

template mkDistinctProcs(T:untyped) =
  proc `steps=`(S:var T, n:int) {.used.} = S.nstep = n

macro runEvolution(cs,steps:untyped,fs:varargs[untyped]):untyped =
  let
    ncs = cs.len - 1
    nfs = fs.len - 1
  if nfs != 1 and nfs != ncs:
    error("runEvolution: requires number of functions be two or the same length of cs.")
  let
    first = fs[0]@cs[0]
    first2 = fs[0]@infix(newLit(2),"*",cs[0])
    rest = newstmtlist()
  for i in 1..<2*ncs:
    let j = nfs - abs(nfs - (i mod (2*nfs)))
    rest.add fs[j]@cs[ncs-abs(ncs-i)]
  result = quote do:
    `first`
    `rest`
    for i in 2..`steps`:
      `first2`
      `rest`
    `first`
  #echo "runEvolution:\n",result.repr

macro mkEvolve(T,nstep,coef:untyped,Fs:varargs[untyped]):untyped =
  let v = gensym(nskVar,"v" & $T)
  result = quote do:
    var `v`:`T`
    `v`.nstep = `nstep`
  for i in 0..<coef.len:
    let csi = coef[i]
    result.add quote do:
      `v`.cs[`i`] = `csi`
  let
    S = gensym(nskParam,"S")
    h = gensym(nskLet,"h")
    n = gensym(nskLet,"n")
  result.add quote do:
    proc evolve(`S`:`T`, t:float) {.used.} =
      let
        `n` = `S`.nstep
        `h` = t/`n`.float
  var cs = newNimNode nnkBracket
  for i in 0..<coef.len:
    let csi = gensym(nskLet,"c" & $i)
    cs.add csi
    result[^1].body.add quote do:
      let `csi` = `S`.cs[`i`]*`h`
  result[^1].body.add quote do:
    runEvolution(`cs`,`n`,`Fs`)
  result.add v
  #echo result.repr

const FinishedAbsT = 4.0  # Some number definitely larger than 1.

template searchMinFor(s:openarray, p:untyped):int =
  var n = -1
  var v:type(s[0])
  for i in 0..<s.len:
    let x {.inject.} = s[i]
    if p:
      n = i
      v = x
      break
  if n >= 0:
    for i in n+1 ..< s.len:
      let x {.inject.} = s[i]
      if p and x < v:
        n = i
        v = x
  n

macro minStepFrom(mds:varargs[untyped]):untyped =
  result = newStmtList()
  let
    n = newLit mds.len
    s = gensym(nskvar,"s")
  result.add quote do:
    var `s`:array[`n`,float]
  for i in 0..<mds.len:
    let md = mds[i]
    result.add quote do:
      `s`[`i`] = `md`.nextSharedStepAbsT
  let pred = quote do:
    x < FinishedAbsT
  result.add getast searchMinFor(s, pred)
  #echo "[info] minStepFrom"
  #echo result.repr

macro mkSharedEvolution*(MDs:varargs[untyped]):untyped =
  ## Combined integrators who share exact one updater.
  ## Restriction:
  ## (a) Each integrator must have the `shared` parameter set properly.
  ## (b) Each integrator must not be created by this `mkSharedEvolution`.
  ## We might use macros to lift restrictions later.
  ## Bugs:
  ## (a) Current implementation breaks if any integrator is a recursive one
  ## where the shared updater is embedded in a deeper layer.
  ## (b) When using multiple steps, intermediate force calculations are repeated
  ## at the end and the beginning of a step, instead of fusing together.
  let
    prep = newStmtList()
    cases = newNimNode(nnkCaseStmt).add quote do:
      minStepFrom `MDs`
    h = gensym(nskLet,"h")
    lastT = gensym(nskVar,"lastT")
    lastAbsT = gensym(nskVar,"lastAbsT")
  for i in 0..<MDs.len:
    let m = MDs[i]
    prep.add quote do:
      `m`.prepareShared
    let f = quote do:
      debug "# MD: ",`i`
      `m`.sharedStep(`lastT`,`lastAbsT`,`h`)
    cases.add newNimNode(nnkOfBranch).add(i.newLit,f)
  cases.add newNimNode(nnkElse).add quote do:
    break
  result = quote do:
    type MDShared = object
      nstep:int
    mkDistinctProcs MDShared
    proc evolve(S:MDShared, t:float) =
      let
        n = S.nstep
        `h` = t/n.float
      for i in 0..<n:
        debug "shared step ",i," of ",S.nstep
        `prep`
        var `lastT`, `lastAbsT` = 0.0
        while true:
          debug "last T, |T| = ",`lastT`," ",`lastAbsT`
          `cases`
    var vMDShared = MDShared(nstep:1)
    vMDShared
  #echo "[info] mkSharedEvolution"
  #echo result.repr

proc nextSharedStepIndex*(X:any, shared:static[int], s:int):tuple[a,b:int] =
  ## Return the absolute time
  const
    ns = X.cs.len
    n1 = ns-1
  if s == 0:
    when shared == 0: result = (a:0, b:1)
    else: result = (a:(-1), b:0)
  else:
    var i = s mod n1
    if shared == 0 and i == 0:  # End of one cycle
      result = (a:0, b:1)
    else:
      i = 2*i - shared  # works for shared == 0 or 1
      result = (a:n1-abs(n1-abs(i)), b:n1-abs(n1-abs(i+1)))

proc nextSharedStep*(X:any, s:int, ix:int):float =
  const
    ns = X.cs.len
    n1 = ns - 1
  let ix = if ix < ns: ix else: 2*n1-ix
  if ix < 0:  # The beginning.
    result = 0.0
  elif ix == 0 and not(s == 0 or s == X.nstep*n1):
    # The end of one cycle not the first or the last.
    result = 2.0*X.cs[ix]/X.nstep.float
  else:
    result = X.cs[ix]/X.nstep.float

proc `~`(x,y:float):bool =
  const CT = 1e-12
  if x == 0: result = abs(y) < CT
  elif y == 0: result = abs(x) < CT
  else:
    let
      ax = abs x
      ay = abs y
    if ax > ay: result = abs(x-y)/ay < CT
    else: result = abs(x-y)/ax < CT
template `!~`(x,y:float):bool = not(x~y)

macro mkSharedEvolve(T:untyped; shared:static[int]; Fs:varargs[untyped]):untyped =
  result = newStmtList()
  let nfs = Fs.len
  if shared < 0: return
  elif shared >= nfs: error("shared value, " & $shared & ", outside of Fs: " & $Fs)
  let
    Fs0 = Fs[0]
    FsS = Fs[shared]
    Fs1 = Fs[1]
  result.add quote do:
    # generate information for shared evolutions
    var
      absTnorm = 0.0
      evoStep = 0
      evoT = 0.0
      evoAbsT = 0.0
      evoAbsTnext = 0.0
    proc prepareShared(S:`T`) =
      evoStep = 0
      evoT = 0.0
      evoAbsT = 0.0
      evoAbsTnext = 0.0
      absTnorm = 0.0
      let nc = S.cs.len-1
      for i in 0..<S.cs.len-`shared`:
        absTnorm += S.cs[nc-abs(nc-(2*i+`shared`))].abs
      debug "sum absT: ",absTnorm
      absTnorm = 1.0/absTnorm
    proc nextSharedStepAbsT(S:`T`):float =
      if evoAbsTnext <= evoAbsT:
        evoAbsTnext = evoAbsT + S.nextSharedStep(evoStep, S.nextSharedStepIndex(`shared`,evoStep).a).abs*absTnorm
      evoAbsTnext
    proc sharedStep(S:`T`; lastT,lastAbsT:var float; t:float) =
      ## Evolve from the lastT to the minimal predetermined step.
      ## Update `lastT` and `lastAbstT`.
      debug "evoStep: ",evoStep," lastT: ",lastT," lastAbsT: ",lastAbsT
      debug "		evoT: ",evoT," evoAbsT: ",evoAbsT," evoAbsTnext: ",evoAbsTnext
      let
        (ix, ixp) = S.nextSharedStepIndex(`shared`,evoStep)
        maxStep = (S.cs.len-1)*S.nstep
      if ix < 0:
        let c1 = S.nextSharedStep(evoStep,ixp)
        debug "	Fs[0] ",c1
        `Fs0`@:(c1*t)
      else:
        let c = S.nextSharedStep(evoStep,ix)
        evoT += c
        evoAbsT = evoAbsTnext
        # FIXME if X is a combination, we need to call its sharedStep
        if evoT !~ lastT:
          debug "	Fs[",`shared`,"] ",(evoT-lastT)
          `FsS`@:((evoT-lastT)*t)
        if evoStep < maxStep or ix > 0:
          let c1 = S.nextSharedStep(evoStep,ixp)
          debug "	Fs[",ixp,"] ",c1
          let c1t = c1*t
          case ixp mod `nfs`:
          of 0: `Fs0`@:c1t
          of 1: `Fs1`@:c1t
        if evoT !~ lastT: lastT = evoT
        if evoAbsT !~ lastAbsT: lastAbsT = evoAbsT
      inc evoStep
      if evoStep > maxStep: evoAbsT = FinishedAbsT
  let cases = result[0][^1].body[^3][^1][0][4][0][^1][^1]  # the CaseStmt above
  for i in 2..<nfs:
    let ofb = cases[^1].copy  # One of the OfBranche
    ofb[0] = i.newLit
    ofb[1][0][1] = Fs[i]
    cases.add ofb.copy
  cases.add newNimNode(nnkElse).add newNimNode(nnkDiscardStmt).add newEmptyNode()
  #echo "[info] mkSharedEvolve"
  #echo result.repr

template mkLeapFrog*(T,V:untyped, steps = 1,
    shared:static[int] = -1):untyped =
  type LeapFrog {.borrow: `.`.} = distinct MDEvolveS1
  mkDistinctProcs LeapFrog
  mkSharedEvolve(LeapFrog, shared, T, V)
  mkEvolve(LeapFrog, steps, [0.5,1], T, V)

template mkOmelyan2MN*(T,V:untyped, steps = 1,
    lambda = 0.1931833275037836,
    shared:static[int] = -1):untyped =
  ## Omelyan et. al. (2002)
  type Omelyan2MN {.borrow: `.`.} = distinct MDEvolveS2
  mkDistinctProcs Omelyan2MN
  mkSharedEvolve(Omelyan2MN, shared, T, V)
  let l = lambda
  mkEvolve(Omelyan2MN, steps, [l,0.5,1.0-2*l], T, V)

template mkSW92*(T,V:untyped, steps = 1):untyped =
  ## Sexton & Weingarten (1992), TVTVT
  mkOmelyan2MN(T,V,steps,1.0/6.0)

template mkOmelyan4MN4FP*(T,V:untyped, steps = 1,
    rho = 0.1786178958448091,
    theta = -0.06626458266981843,
    lambda = 0.7123418310626056,
    shared:static[int] = -1):untyped =
  ## Omelyan et. al. (2003)
  type Omelyan4MN4FP {.borrow: `.`.} = distinct MDEvolveS4
  mkDistinctProcs Omelyan4MN4FP
  mkSharedEvolve(Omelyan4MN4FP, shared, T, V)
  let
    l = lambda
    r = rho
    t = theta
  mkEvolve(Omelyan4MN4FP, steps, [r, l, t, 0.5-l, 1.0-2*(t+r)], T, V)

template mkOmelyan4MN5FV*(V,T:untyped, steps = 1,
    theta = 0.08398315262876693,
    rho = 0.2539785108410595,
    lambda = 0.6822365335719091,
    mu = -0.03230286765269967,
    shared:static[int] = -1):untyped =
  ## Omelyan et. al. (2003)
  type Omelyan4MN5FV {.borrow: `.`.} = distinct MDEvolveS5
  mkDistinctProcs Omelyan4MN5FV
  mkSharedEvolve(Omelyan4MN5FV, shared, V, T)
  let
    l = lambda
    r = rho
    t = theta
    m = mu
  mkEvolve(Omelyan4MN5FV, steps, [t, r, l, m, 0.5-(l+t), 1.0-2*(m+r)], V, T)

template mkFGYin11*(V,T,Vfg,save,load:untyped, steps = 1,
    shared:static[int] = -1):untyped =
  ## Force Gradient Integrator, H. Yin (2011)
  ## It is 4th order only when T = p^2/2.
  ## In addition to the usual updater, `V` and `T`, we need:
  ## `Vfg` that updates `x` using the force directly;
  ## `save` and `load` that save and load `x`.
  proc fVfg(t:float) {.gensym.} =
    let h = (3.0/32.0)*t*t  # 1r24=3r32**:2r3
    save
    Vfg@:h
    V@:t
    load
  type FGYin11 {.borrow: `.`.} = distinct MDEvolveS2
  mkDistinctProcs FGYin11
  mkSharedEvolve(FGYin11, shared, V, T, fVfg)
  mkEvolve(FGYin11, steps, [1.0/6.0, 0.5, 2.0/3.0], V, T, fVfg)
