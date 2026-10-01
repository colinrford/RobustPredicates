//
//  InCircleTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

@Suite("inCircle")
struct InCircleTests {

  @Test func squareCornerGridMatchesOracle() {
    let a = SIMD2(1.0, 1.0), b = SIMD2(3.0, 1.0), c = SIMD2(3.0, 3.0)
    var wrong: [SIMD2<Double>] = []
    var naiveWrong = 0
    for d in ulpGrid(from: SIMD2(1.0, 3.0), size: 64) {
      let want = inCircleOracle(a, b, c, d)
      if inCircle(a, b, c, d) != want || CirclePosition(sign: inCircleExact(a, b, c, d)) != want {
        wrong.append(d)
      }
      if naiveInCircle(a, b, c, d) != want { naiveWrong += 1 }
    }
    #expect(wrong.isEmpty, "first failures: \(wrong.prefix(3))")
    #expect(naiveWrong > 0, "the grid must defeat the unfiltered determinant")
  }

  @Test(arguments: [121 as UInt64, 122, 123])
  func nearCocircularGridsMatchOracle(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var wrong = 0
    var naiveWrong = 0
    for _ in 0..<20 {
      let (a, b, c, d0) = nearCocircularQuad(&rng)
      for d in ulpGrid(from: d0, size: 16) {
        let want = inCircleOracle(a, b, c, d)
        if inCircle(a, b, c, d) != want || CirclePosition(sign: inCircleExact(a, b, c, d)) != want {
          wrong += 1
        }
        if naiveInCircle(a, b, c, d) != want { naiveWrong += 1 }
      }
    }
    #expect(wrong == 0)
    #expect(naiveWrong > 0, "the grids must defeat the unfiltered determinant")
  }

  // Each transform is exact, so an exact predicate's answer changes predictably.
  @Test(arguments: [131 as UInt64, 132])
  func exactSymmetriesHoldNearDegeneracy(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var seen: Set<CirclePosition> = []
    for _ in 0..<10 {
      let (a, b, c, d0) = nearCocircularQuad(&rng)
      for d in ulpGrid(from: d0, size: 8) {
        let p = inCircle(a, b, c, d)
        seen.insert(p)
        #expect(inCircle(b, c, a, d) == p)
        #expect(inCircle(b, a, c, d) == p.flipped)
        #expect(inCircle(d, b, c, a) == p.flipped)
        #expect(inCircle(-a, -b, -c, -d) == p)
        #expect(inCircle(swappedXY(a), swappedXY(b), swappedXY(c), swappedXY(d)) == p.flipped)
        for k in [-64, 37, 64] {
          let s = powerOfTwo(k)
          #expect(inCircle(a * s, b * s, c * s, d * s) == p)
        }
      }
    }
    #expect(seen.isSuperset(of: [.inside, .outside]), "a constant answer would satisfy every identity")
  }

  @Test func rectangleCornersAreCocircular() {
    for w in 1...6 {
      for h in 1...6 {
        let a = SIMD2(0.0, 0.0)
        let b = SIMD2(Double(w), 0.0)
        let c = SIMD2(Double(w), Double(h))
        let d = SIMD2(0.0, Double(h))
        #expect(inCircle(a, b, c, d) == .on, "rect \(w)×\(h)")
        #expect(inCircle(a, b, c, SIMD2(Double(w) / 2, Double(h) / 2)) == .inside)
        #expect(inCircle(a, b, c, SIMD2(-Double(w + h), -Double(w + h))) == .outside)
      }
    }
  }

  @Test func latticePointsOnARadiusFiveCircleAreCocircular() {
    let offsets: [SIMD2<Double>] = [
      SIMD2(5, 0), SIMD2(4, 3), SIMD2(3, 4), SIMD2(0, 5), SIMD2(-3, 4), SIMD2(-4, 3),
      SIMD2(-5, 0), SIMD2(-4, -3), SIMD2(-3, -4), SIMD2(0, -5), SIMD2(3, -4), SIMD2(4, -3),
    ]
    let center = SIMD2(7.0, -11.0)
    let points = offsets.map { $0 + center }
    for i in points.indices {
      for j in points.indices where j != i {
        for k in points.indices where k != i && k != j {
          let ccw = orient2d(points[i], points[j], points[k]) == .ccw
          #expect(inCircle(points[i], points[j], points[k], center) == (ccw ? .inside : .outside))
          for l in points.indices where l != i && l != j && l != k {
            #expect(inCircle(points[i], points[j], points[k], points[l]) == .on)
          }
        }
      }
    }
  }

  @Test func oneUlpResolvesCocircularity() {
    let a = SIMD2(0.0, 0.0)
    let b = SIMD2(1.0, 0.0)
    let c = SIMD2(1.0, 1.0)
    #expect(inCircle(a, b, c, SIMD2(0.0, 1.0)) == .on)
    #expect(inCircle(a, b, c, SIMD2(0.0, 1.0 + 1.0.ulp)) == .outside)
    #expect(inCircle(a, b, c, SIMD2(0.0, 1.0 - 0.5.ulp)) == .inside)
  }

  @Test(arguments: [11 as UInt64, 12, 13])
  func randomIntegerInputsMatchOracle(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<1000 {
      let a = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      let b = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      let c = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      let d = randomIntegerPoint(&rng, in: -1_000_000...1_000_000)
      #expect(inCircle(a, b, c, d) == inCircleOracle(a, b, c, d))
    }
  }
}
