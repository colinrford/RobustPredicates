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

/// `a · b` as `hi + lo` exactly; `lo` comes from one fused multiply-add
/// (TwoProductFMA: Ogita, Rump & Oishi 2005, Algorithm 3.5).
@inlinable
func twoProd(_ a: Double, _ b: Double) -> (hi: Double, lo: Double) {
  let p = a * b
  return (p, (-p).addingProduct(a, b))
}

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
    if abs(e[i]) < abs(f[j]) { g.append(e[i]); i += 1 } else { g.append(f[j]); j += 1 }
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

/// The largest nonzero component, which carries the expansion's sign; 0 if the expansion is zero.
func expansionSign(_ e: [Double]) -> Double {
  for c in e.reversed() where c != 0 { return c }
  return 0
}
