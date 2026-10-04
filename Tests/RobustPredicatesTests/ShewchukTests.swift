//
//  ShewchukTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

private struct Agreement<Outcome: Hashable> {
  var outcomes: Set<Outcome> = []
  var disagreements = 0
  var naiveWrong = 0

  mutating func record(_ got: Outcome, shewchuk: Outcome, naive: Outcome) {
    outcomes.insert(shewchuk)
    if got != shewchuk { disagreements += 1 }
    if naive != shewchuk { naiveWrong += 1 }
  }
}

// The same near-degenerate and degenerate inputs as the oracle tests, checked
// against Shewchuk's predicates.c.
@Suite("predicates.c")
struct ShewchukTests {

  @Test func orient2dAgrees() {
    var t = Agreement<Orientation>()
    let q = SIMD2(12.0, 12.0), r = SIMD2(24.0, 24.0)
    for p in ulpGrid(from: SIMD2(0.5, 0.5), size: 256) {
      t.record(orient2d(p, q, r), shewchuk: shewchukOrient2d(p, q, r), naive: naiveOrient2d(p, q, r))
    }
    for seed: UInt64 in 101...103 {
      var rng = SplitMix64(seed: seed)
      for _ in 0..<20 {
        let (a, b, c0) = nearCollinearTriple(&rng)
        for c in ulpGrid(from: c0, size: 16) {
          t.record(orient2d(a, b, c), shewchuk: shewchukOrient2d(a, b, c), naive: naiveOrient2d(a, b, c))
        }
      }
    }
    #expect(t.disagreements == 0)
    #expect(t.outcomes == [.ccw, .collinear, .cw])
    #expect(t.naiveWrong > 0)
  }

  @Test func inCircleAgrees() {
    var t = Agreement<CirclePosition>()
    let a = SIMD2(1.0, 1.0), b = SIMD2(3.0, 1.0), c = SIMD2(3.0, 3.0)
    for d in ulpGrid(from: SIMD2(1.0, 3.0), size: 64) {
      t.record(inCircle(a, b, c, d), shewchuk: shewchukInCircle(a, b, c, d), naive: naiveInCircle(a, b, c, d))
    }
    for seed: UInt64 in 121...123 {
      var rng = SplitMix64(seed: seed)
      for _ in 0..<20 {
        let (a, b, c, d0) = nearCocircularQuad(&rng)
        for d in ulpGrid(from: d0, size: 16) {
          t.record(inCircle(a, b, c, d), shewchuk: shewchukInCircle(a, b, c, d), naive: naiveInCircle(a, b, c, d))
        }
      }
    }
    for w in 1...6 {
      for h in 1...6 {
        let a = SIMD2(0.0, 0.0), b = SIMD2(Double(w), 0.0)
        let c = SIMD2(Double(w), Double(h)), d = SIMD2(0.0, Double(h))
        t.record(inCircle(a, b, c, d), shewchuk: shewchukInCircle(a, b, c, d), naive: naiveInCircle(a, b, c, d))
      }
    }
    #expect(t.disagreements == 0)
    #expect(t.outcomes == [.inside, .on, .outside])
    #expect(t.naiveWrong > 0)
  }

  @Test func orient3dAgrees() {
    var t = Agreement<PlaneSide>()
    for seed: UInt64 in 141...143 {
      var rng = SplitMix64(seed: seed)
      for _ in 0..<20 {
        let (a, b, c, d0) = nearCoplanarQuad(&rng)
        for xz in ulpGrid(from: SIMD2(d0.x, d0.z), size: 16) {
          let d = SIMD3(xz.x, d0.y, xz.y)
          t.record(orient3d(a, b, c, d), shewchuk: shewchukOrient3d(a, b, c, d), naive: naiveOrient3d(a, b, c, d))
        }
      }
    }
    var rng = SplitMix64(seed: 61)
    for _ in 0..<500 {
      let z = Double(Int64.random(in: -1000...1000, using: &rng))
      func onPlane() -> SIMD3<Double> {
        let p = randomIntegerPoint(&rng, in: -1000...1000)
        return SIMD3(p.x, p.y, z)
      }
      let (a, b, c, d) = (onPlane(), onPlane(), onPlane(), onPlane())
      t.record(orient3d(a, b, c, d), shewchuk: shewchukOrient3d(a, b, c, d), naive: naiveOrient3d(a, b, c, d))
    }
    #expect(t.disagreements == 0)
    #expect(t.outcomes == [.above, .on, .below])
    #expect(t.naiveWrong > 0)
  }
}
