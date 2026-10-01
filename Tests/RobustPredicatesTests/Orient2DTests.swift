//
//  Orient2DTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

@Suite("orient2d")
struct Orient2DTests {
  
  // Kettner, Mehlhorn, Pion, Schirra & Yap 2008, §4: p = (0.5 + i·u, 0.5 + j·u)
  // with u = 2⁻⁵³, against q = (12, 12) and r = (24, 24).
  @Test func kettnerGridMatchesOracle() {
    let q = SIMD2(12.0, 12.0), r = SIMD2(24.0, 24.0)
    var wrong: [SIMD2<Double>] = []
    var naiveWrong = 0
    for p in ulpGrid(from: SIMD2(0.5, 0.5), size: 256) {
      let want = orient2dOracle(p, q, r)
      if orient2d(p, q, r) != want || Orientation(sign: orient2dExact(p, q, r)) != want {
        wrong.append(p)
      }
      if naiveOrient2d(p, q, r) != want { naiveWrong += 1 }
    }
    #expect(wrong.isEmpty, "first failures: \(wrong.prefix(3))")
    #expect(naiveWrong > 0, "the grid must defeat the unfiltered determinant")
  }
  
  @Test(arguments: [101 as UInt64, 102, 103])
  func nearCollinearGridsMatchOracle(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var wrong = 0
    var naiveWrong = 0
    for _ in 0..<20 {
      let (a, b, c0) = nearCollinearTriple(&rng)
      for c in ulpGrid(from: c0, size: 16) {
        let want = orient2dOracle(a, b, c)
        if orient2d(a, b, c) != want || Orientation(sign: orient2dExact(a, b, c)) != want {
          wrong += 1
        }
        if naiveOrient2d(a, b, c) != want { naiveWrong += 1 }
      }
    }
    #expect(wrong == 0)
    #expect(naiveWrong > 0, "the grids must defeat the unfiltered determinant")
  }
  
  // Each transform is exact, so an exact predicate's answer changes predictably.
  @Test(arguments: [111 as UInt64, 112])
  func exactSymmetriesHoldNearDegeneracy(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var seen: Set<Orientation> = []
    for _ in 0..<10 {
      let (a, b, c0) = nearCollinearTriple(&rng)
      for c in ulpGrid(from: c0, size: 8) {
        let o = orient2d(a, b, c)
        seen.insert(o)
        #expect(orient2d(b, c, a) == o)
        #expect(orient2d(c, a, b) == o)
        #expect(orient2d(b, a, c) == o.flipped)
        #expect(orient2d(-a, -b, -c) == o)
        #expect(orient2d(swappedXY(a), swappedXY(b), swappedXY(c)) == o.flipped)
        for k in [-64, 37, 64] {
          let s = powerOfTwo(k)
          #expect(orient2d(a * s, b * s, c * s) == o)
        }
      }
    }
    #expect(seen.isSuperset(of: [.ccw, .cw]), "a constant answer would satisfy every identity")
  }
  
  @Test(arguments: [21 as UInt64, 22])
  func collinearIntegerTriplesAreCollinear(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<500 {
      let a = randomIntegerPoint(&rng, in: -1000...1000)
      let v = randomIntegerPoint(&rng, in: -50...50)
      if v == .zero { continue }
      let t = Double(Int.random(in: 1...20, using: &rng))
      let s = Double(Int.random(in: 21...40, using: &rng))
      #expect(orient2d(a, a + t * v, a + s * v) == .collinear)
      // (t·v) × (s·v + v⊥) = t·|v|² > 0.
      #expect(orient2d(a, a + t * v, a + s * v + SIMD2(-v.y, v.x)) == .ccw)
    }
  }
  
  // At this magnitude the products are exact in Double, so this only checks
  // the fast path's arithmetic; the grids above exercise the exact fallback.
  @Test(arguments: [7 as UInt64, 8, 9])
  func randomIntegerInputsMatchOracle(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<2000 {
      let a = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      let b = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      let c = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      #expect(orient2d(a, b, c) == orient2dOracle(a, b, c))
    }
  }
}
