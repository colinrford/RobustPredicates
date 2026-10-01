//
//  Orient3DTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

/// Lifting d by +z puts it on the side the normal (b−a)×(c−a) points to exactly
/// when a, b, c run counterclockwise seen from +z, i.e. in their xy projection.
private func offPlaneSideIsConsistent(
  _ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>
) -> Bool {
  let xy = { (p: SIMD3<Double>) in SIMD2(p.x, p.y) }
  let side = orient3d(a, b, c, d + SIMD3(0, 0, 1))
  switch orient2d(xy(a), xy(b), xy(c)) {
  case .ccw: return side == .above
  case .cw: return side == .below
  case .collinear: return true
  }
}

@Suite("orient3d")
struct Orient3DTests {
  
  @Test(arguments: [141 as UInt64, 142, 143])
  func nearCoplanarGridsMatchOracle(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var wrong = 0
    var naiveWrong = 0
    for _ in 0..<20 {
      let (a, b, c, d0) = nearCoplanarQuad(&rng)
      for xz in ulpGrid(from: SIMD2(d0.x, d0.z), size: 16) {
        let d = SIMD3(xz.x, d0.y, xz.y)
        let want = orient3dOracle(a, b, c, d)
        if orient3d(a, b, c, d) != want || PlaneSide(sign: orient3dExact(a, b, c, d)) != want {
          wrong += 1
        }
        if naiveOrient3d(a, b, c, d) != want { naiveWrong += 1 }
      }
    }
    #expect(wrong == 0)
    #expect(naiveWrong > 0, "the grids must defeat the unfiltered determinant")
  }
  
  // Each transform is exact, so an exact predicate's answer changes predictably.
  @Test(arguments: [151 as UInt64, 152])
  func exactSymmetriesHoldNearDegeneracy(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var seen: Set<PlaneSide> = []
    for _ in 0..<10 {
      let (a, b, c, d0) = nearCoplanarQuad(&rng)
      for xz in ulpGrid(from: SIMD2(d0.x, d0.z), size: 8) {
        let d = SIMD3(xz.x, d0.y, xz.y)
        let o = orient3d(a, b, c, d)
        seen.insert(o)
        #expect(orient3d(b, c, a, d) == o)
        #expect(orient3d(b, a, c, d) == o.flipped)
        #expect(orient3d(d, b, c, a) == o.flipped)
        #expect(orient3d(-a, -b, -c, -d) == o.flipped)
        #expect(orient3d(swappedXY(a), swappedXY(b), swappedXY(c), swappedXY(d)) == o.flipped)
        for k in [-64, 37, 64] {
          let s = powerOfTwo(k)
          #expect(orient3d(a * s, b * s, c * s, d * s) == o)
        }
      }
    }
    #expect(seen.isSuperset(of: [.above, .below]), "a constant answer would satisfy every identity")
  }
  
  @Test(arguments: [61 as UInt64, 62])
  func constantZPointsAreCoplanar(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<500 {
      let z = Double(Int64.random(in: -1000...1000, using: &rng))
      func onPlane() -> SIMD3<Double> {
        let p = randomIntegerPoint(&rng, in: -1000...1000)
        return SIMD3(p.x, p.y, z)
      }
      let (a, b, c, d) = (onPlane(), onPlane(), onPlane(), onPlane())
      #expect(orient3d(a, b, c, d) == .on)
      #expect(offPlaneSideIsConsistent(a, b, c, d))
    }
  }
  
  @Test(arguments: [63 as UInt64, 64])
  func tiltedIntegerPlanePointsAreCoplanar(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<500 {
      let k = Double(Int64.random(in: -1000...1000, using: &rng))
      func onPlane() -> SIMD3<Double> {
        let p = randomIntegerPoint(&rng, in: -1000...1000)
        return SIMD3(p.x, p.y, k - p.x - p.y)
      }
      let (a, b, c, d) = (onPlane(), onPlane(), onPlane(), onPlane())
      #expect(orient3d(a, b, c, d) == .on)
      #expect(offPlaneSideIsConsistent(a, b, c, d))
    }
  }
  
  @Test(arguments: [51 as UInt64, 52, 53])
  func randomIntegerInputsMatchOracle(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<1000 {
      let a = randomIntegerPoint3(&rng, in: -1_000_000...1_000_000)
      let b = randomIntegerPoint3(&rng, in: -1_000_000...1_000_000)
      let c = randomIntegerPoint3(&rng, in: -1_000_000...1_000_000)
      let d = randomIntegerPoint3(&rng, in: -1_000_000...1_000_000)
      #expect(orient3d(a, b, c, d) == orient3dOracle(a, b, c, d))
    }
  }
}
