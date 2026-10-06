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

/// The Kettner grid and the near-collinear grids.
private func orient2dInputs() -> [(SIMD2<Double>, SIMD2<Double>, SIMD2<Double>)] {
  let q = SIMD2(12.0, 12.0), r = SIMD2(24.0, 24.0)
  var inputs = ulpGrid(from: SIMD2(0.5, 0.5), size: 256).map { ($0, q, r) }
  for seed: UInt64 in 101...103 {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<20 {
      let (a, b, c0) = nearCollinearTriple(&rng)
      inputs += ulpGrid(from: c0, size: 16).map { (a, b, $0) }
    }
  }
  return inputs
}

/// The square-corner grid, the near-cocircular grids, the rectangle corners, and
/// grids around the fourth corner of rectangles whose sides mix magnitudes, so
/// coordinate differences round and stage D of the adaptive inCircle runs.
private func inCircleInputs() -> [(SIMD2<Double>, SIMD2<Double>, SIMD2<Double>, SIMD2<Double>)] {
  let a = SIMD2(1.0, 1.0), b = SIMD2(3.0, 1.0), c = SIMD2(3.0, 3.0)
  var inputs = ulpGrid(from: SIMD2(1.0, 3.0), size: 64).map { (a, b, c, $0) }
  for seed: UInt64 in 121...123 {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<20 {
      let (a, b, c, d0) = nearCocircularQuad(&rng)
      inputs += ulpGrid(from: d0, size: 16).map { (a, b, c, $0) }
    }
  }
  for w in 1...6 {
    for h in 1...6 {
      inputs.append((SIMD2(0.0, 0.0), SIMD2(Double(w), 0.0), SIMD2(Double(w), Double(h)), SIMD2(0.0, Double(h))))
    }
  }
  for seed: UInt64 in 124...126 {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<20 {
      let x0 = randomDouble(&rng, exponents: -20...0), x1 = randomDouble(&rng, exponents: 0...20)
      let y0 = randomDouble(&rng, exponents: -20...0), y1 = randomDouble(&rng, exponents: 0...20)
      let a = SIMD2(x0, y0), b = SIMD2(x1, y0), c = SIMD2(x1, y1)
      // Rotated so each tail is nonzero somewhere.
      for d in ulpGrid(from: SIMD2(x0, y1), size: 4) {
        inputs += [(a, b, c, d), (b, c, a, d), (c, a, b, d)]
      }
    }
  }
  return inputs
}

/// The near-coplanar grids, integer points on horizontal planes, and grids around
/// a point on the plane x = y (and its axis permutations) with coordinates of mixed
/// magnitude, so coordinate differences round and stage D of the adaptive orient3d runs.
private func orient3dInputs() -> [(SIMD3<Double>, SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)] {
  var inputs: [(SIMD3<Double>, SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)] = []
  for seed: UInt64 in 141...143 {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<20 {
      let (a, b, c, d0) = nearCoplanarQuad(&rng)
      inputs += ulpGrid(from: SIMD2(d0.x, d0.z), size: 16).map { (a, b, c, SIMD3($0.x, d0.y, $0.y)) }
    }
  }
  var rng = SplitMix64(seed: 61)
  for _ in 0..<500 {
    let z = Double(Int64.random(in: -1000...1000, using: &rng))
    func onPlane() -> SIMD3<Double> {
      let p = randomIntegerPoint(&rng, in: -1000...1000)
      return SIMD3(p.x, p.y, z)
    }
    inputs.append((onPlane(), onPlane(), onPlane(), onPlane()))
  }
  for seed: UInt64 in 144...146 {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<20 {
      func onDiagonal() -> (u: Double, z: Double) {
        (randomDouble(&rng, exponents: -20...20), randomDouble(&rng, exponents: -20...20))
      }
      let (a, b, c, d) = (onDiagonal(), onDiagonal(), onDiagonal(), onDiagonal())
      for g in ulpGrid(from: SIMD2(d.u, d.z), size: 4) {
        for point in [{ (u: Double, v: Double, z: Double) in SIMD3(u, v, z) },
                      { (u: Double, v: Double, z: Double) in SIMD3(u, z, v) },
                      { (u: Double, v: Double, z: Double) in SIMD3(z, u, v) }] {
          inputs.append((point(a.u, a.u, a.z), point(b.u, b.u, b.z), point(c.u, c.u, c.z), point(g.x, d.u, g.y)))
        }
      }
    }
  }
  return inputs
}

// The same near-degenerate and degenerate inputs as the oracle tests, checked
// against Shewchuk's predicates.c.
@Suite("predicates.c")
struct ShewchukTests {

  @Test func orient2dAgrees() {
    var t = Agreement<Orientation>()
    for (a, b, c) in orient2dInputs() {
      t.record(orient2d(a, b, c), shewchuk: shewchukOrient2d(a, b, c), naive: naiveOrient2d(a, b, c))
    }
    #expect(t.disagreements == 0)
    #expect(t.outcomes == [.ccw, .collinear, .cw])
    #expect(t.naiveWrong > 0)
  }

  // Each stage returns its own approximation, so a wrong bound or stage
  // changes the value even when the sign survives.
  @Test func orient2dAdaptMatchesBitForBit() {
    var mismatches = 0
    for (a, b, c) in orient2dInputs() {
      let detsum = abs((a.x - c.x) * (b.y - c.y)) + abs((a.y - c.y) * (b.x - c.x))
      let ours = orient2dAdapt(a, b, c, detsum: detsum)
      let theirs = shewchukOrient2dAdapt(a, b, c, detsum: detsum)
      if ours.bitPattern != theirs.bitPattern { mismatches += 1 }
    }
    #expect(mismatches == 0)
  }

  @Test func inCircleAgrees() {
    var t = Agreement<CirclePosition>()
    for (a, b, c, d) in inCircleInputs() {
      t.record(inCircle(a, b, c, d), shewchuk: shewchukInCircle(a, b, c, d), naive: naiveInCircle(a, b, c, d))
    }
    #expect(t.disagreements == 0)
    #expect(t.outcomes == [.inside, .on, .outside])
    #expect(t.naiveWrong > 0)
  }

  @Test func inCircleAdaptMatchesBitForBit() {
    var mismatches = 0
    for (a, b, c, d) in inCircleInputs() {
      let adx = a.x - d.x, ady = a.y - d.y
      let bdx = b.x - d.x, bdy = b.y - d.y
      let cdx = c.x - d.x, cdy = c.y - d.y
      let permanent = (abs(bdx * cdy) + abs(cdx * bdy)) * (adx * adx + ady * ady)
        + (abs(cdx * ady) + abs(adx * cdy)) * (bdx * bdx + bdy * bdy)
        + (abs(adx * bdy) + abs(bdx * ady)) * (cdx * cdx + cdy * cdy)
      let ours = inCircleAdapt(a, b, c, d, permanent: permanent)
      let theirs = shewchukInCircleAdapt(a, b, c, d, permanent: permanent)
      if ours.bitPattern != theirs.bitPattern { mismatches += 1 }
    }
    #expect(mismatches == 0)
  }

  @Test func orient3dAgrees() {
    var t = Agreement<PlaneSide>()
    for (a, b, c, d) in orient3dInputs() {
      t.record(orient3d(a, b, c, d), shewchuk: shewchukOrient3d(a, b, c, d), naive: naiveOrient3d(a, b, c, d))
    }
    #expect(t.disagreements == 0)
    #expect(t.outcomes == [.above, .on, .below])
    #expect(t.naiveWrong > 0)
  }

  @Test func orient3dAdaptMatchesBitForBit() {
    var mismatches = 0
    for (a, b, c, d) in orient3dInputs() {
      let adx = a.x - d.x, ady = a.y - d.y, adz = a.z - d.z
      let bdx = b.x - d.x, bdy = b.y - d.y, bdz = b.z - d.z
      let cdx = c.x - d.x, cdy = c.y - d.y, cdz = c.z - d.z
      let permanent = (abs(bdx * cdy) + abs(cdx * bdy)) * abs(adz)
        + (abs(cdx * ady) + abs(adx * cdy)) * abs(bdz)
        + (abs(adx * bdy) + abs(bdx * ady)) * abs(cdz)
      let ours = orient3dAdapt(a, b, c, d, permanent: permanent)
      let theirs = shewchukOrient3dAdapt(a, b, c, d, permanent: permanent)
      if ours.bitPattern != theirs.bitPattern { mismatches += 1 }
    }
    #expect(mismatches == 0)
  }
}
