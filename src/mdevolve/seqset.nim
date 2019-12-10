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

type seqset*[T] = object
  s:seq[T]

func `$`*[T](x:seqset[T]):string = $x.s

func newseqset*[T](xs:varargs[T]):seqset[T] =
  result.s.newseq(xs.len)
  for i in 0..<xs.len:
    result.s[i] = xs[i]

iterator items*[T](x:seqset[T]):T =
  for c in x.s: yield c

func len*(x:seqset):int = x.s.len

proc empty*[T](x:var seqset[T]) = x.s.setlen 0

func contains*[T](x:seqset[T], e:T):bool =
  for c in x:
    if e == c:
      return true
  return false

func `[]`*[T](x:seqset[T], i:int):T = x.s[i]

proc add*[T](x:var seqset[T], e:T) =
  if e notin x:
    x.s.add e
proc add*[T](x:var seqset[T], y:seqset[T]) =
  for c in y: x.add c

proc addReturnIntersection*[T](x:var seqset[T], y:seqset[T]):seqset[T] =
  for c in y:
    if c in x:
      result.s.add c
    else:
      x.s.add c

func intersection*[T](x,y:seqset[T]):seqset[T] =
  for c in y:
    if c in x:
      result.s.add c
