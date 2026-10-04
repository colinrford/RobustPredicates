//
//  ExactArithmetic.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

// MARK: - Error-free transformations

/// `a + b` as `hi + lo` exactly, where `hi` is the rounded sum (Knuth, TAOCP Vol. 2, §4.2.2).
@inlinable
func twoSum(_ a: Double, _ b: Double) -> (hi: Double, lo: Double) {
  let s = a + b
  let bv = s - a
  let av = s - bv
  return (s, (a - av) + (b - bv))
}

/// `a - b` as `hi + lo` exactly, where `hi` is the rounded difference (Knuth, TAOCP Vol. 2, §4.2.2).
@inlinable
func twoDiff(_ a: Double, _ b: Double) -> (hi: Double, lo: Double) {
  let s = a - b
  let bv = a - s
  let av = s + bv
  return (s, (a - av) + (bv - b))
}

/// `a − b` as an expansion: `[lo, hi]`, or `[hi]` when the difference is exact.
@inlinable
func twoDiffE(_ a: Double, _ b: Double) -> [Double] {
  let (hi, lo) = twoDiff(a, b)
  return lo == 0 ? [hi] : [lo, hi]
}

/// placeholder; rounding error of `x = a - b`
@inlinable
func twoDiffTail(_ a: Double, _ b: Double, _ x: Double) -> Double {
  let bvirt = (a - x)
  let avirt = x + bvirt
  return (a - avirt) + (bvirt - b)
}

/// `a · b` as `hi + lo` exactly; `lo` comes from one fused multiply-add
/// (TwoProductFMA: Ogita, Rump & Oishi 2005, Algorithm 3.5).
@inlinable
func twoProd(_ a: Double, _ b: Double) -> (hi: Double, lo: Double) {
  let p = a * b
  return (p, (-p).addingProduct(a, b))
}

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

/// placeholder
@inlinable
func twoOneSum(_ a1: Double, _ a0: Double, _ b: Double) -> (Double, Double, Double) {
  let (i, x0) = twoSum(a0, b)
  let (x2, x1) = twoSum(a1, i)
  return (x2, x1, x0)
}

/// placeholder;
@inlinable
func twoOneDiff(_ a1: Double, _ a0: Double, _ b: Double) -> (Double, Double, Double) {
  let (i, x0) = twoDiff(a0, b)
  let (x2, x1) = twoSum(a1, i)
  return (x2, x1, x0)
}

/// placeholder
@inlinable
func twoTwoSum(_ a1: Double, _ a0: Double, _ b1: Double, _ b0: Double) -> (Double, Double, Double, Double) {
  let (j, z, x0) = twoOneSum(a1, a0, b0)
  let (x3, x2, x1) = twoOneSum(j, z, b1)
  return (x3, x2, x1, x0)
}

/// placeholder
@inlinable
func twoTwoDiff(_ a1: Double, _ a0: Double, _ b1: Double, _ b0: Double) -> (Double, Double, Double, Double) {
  let (j, z, x0) = twoOneDiff(a1, a0, b0)
  let (x3, x2, x1) = twoOneDiff(j, z, b1)
  return (x3, x2, x1, x0)
}

/*
@inlinable
Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b , _j, x1, x0); \
  Two_One_Sum(a3, a2, _j, x4, x3, x2)

@inlinable
Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0); \
  Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)

@inlinable
Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, \
                      x1, x0) \
  Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0); \
  Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)

@inlinable
Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, \
                      x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b , _j, x3, x2, x1, x0); \
  Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)

@inlinable
Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, \
                      x6, x5, x4, x3, x2, x1, x0) \
  Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, \
                _1, _0, x0); \
  Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, \
                x3, x2, x1)

@inlinable
Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, \
                       x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0) \
  Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, \
                _2, _1, _0, x1, x0); \
  Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, \
                x7, x6, x5, x4, x3, x2) */

// MARK: - Expansions
// An expansion is an array of Doubles whose exact sum is the value it represents:
// nonoverlapping, in increasing magnitude, with zeros removed (Shewchuk 1997, §2).

/// Sum of two expansions (Shewchuk's fast-expansion-sum with zero elimination).
func expansionSum(_ e: [Double], _ f: [Double]) -> [Double] {
  if e.isEmpty { return f }
  if f.isEmpty { return e }
  // Merge by increasing magnitude.
  var g: [Double] = []
  g.reserveCapacity(e.count + f.count)
  var (i, j) = (0, 0)
  while i < e.count && j < f.count {
    if (f[j] > e[i]) == (f[j] > -e[i]) {
      g.append(e[i]); i += 1
    } else {
      g.append(f[j]); j += 1
    }
  }
  g.append(contentsOf: e[i...])
  g.append(contentsOf: f[j...])
  
  var h: [Double] = []
  h.reserveCapacity(g.count)
  let (q, small) = twoSum(g[1], g[0])
  if small != 0 { h.append(small) }
  var acc = q
  for k in 2..<g.count {
    let (sum, err) = twoSum(acc, g[k])
    acc = sum
    if err != 0 { h.append(err) }
  }
  if acc != 0 || h.isEmpty { h.append(acc) }
  return h
}

func expansionNegate(_ e: [Double]) -> [Double] { e.map { -$0 } }

/// Expansion × scalar (Shewchuk's scale-expansion with zero elimination).
func expansionScale(_ e: [Double], _ b: Double) -> [Double] {
  guard !e.isEmpty else { return [] }
  var h: [Double] = []
  h.reserveCapacity(e.count * 2)
  var (q, small) = twoProd(e[0], b)
  if small != 0 { h.append(small) }
  for i in 1..<e.count {
    let (p, err1) = twoProd(e[i], b)
    let (sum, err2) = twoSum(q, err1)
    if err2 != 0 { h.append(err2) }
    let (newQ, err3) = twoSum(p, sum)
    if err3 != 0 { h.append(err3) }
    q = newQ
  }
  if q != 0 || h.isEmpty { h.append(q) }
  return h
}

/// Exact product of two expansions (distribute + accumulate).
func expansionProduct(_ e: [Double], _ f: [Double]) -> [Double] {
  var acc: [Double] = []
  for b in f {
    acc = expansionSum(acc, expansionScale(e, b))
  }
  return acc
}

/// The most significant nonzero component, whose sign is the expansion's (Shewchuk 1997, §2.8);
/// 0 if the expansion is zero.
func mostSignificantComponent(_ e: [Double]) -> Double {
  for c in e.reversed() where c != 0 { return c }
  return 0
}

/// Produce a one-word estimate of an expansion's value.
func estimate(_ e: [Double]) -> Double {
  e.dropFirst().reduce(e[0], +)
}
