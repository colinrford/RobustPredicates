//
//  InputRangeTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

/// Coordinates that are 0 or have magnitude in [2⁻ˡ, 2ᵘ].
private struct Domain {
  let u, l: Int

  func contains(_ xs: [Double]) -> Bool {
    xs.allSatisfy { $0 == 0 || (powerOfTwo(-l)...powerOfTwo(u)).contains(abs($0)) }
  }

  func reachesUpper(_ xs: [Double]) -> Bool { xs.contains { abs($0) >= powerOfTwo(u - 3) } }
  func reachesLower(_ xs: [Double]) -> Bool { xs.contains { $0 != 0 && abs($0) <= powerOfTwo(3 - l) } }
}

private struct Tally<Outcome: Hashable> {
  var coordinates: [Double] = []
  var outcomes: Set<Outcome> = []
  var wrong = 0
  var naiveWrong = 0

  mutating func record(_ got: Outcome, naive: Outcome, want: Outcome) {
    outcomes.insert(want)
    if got != want { wrong += 1 }
    if naive != want { naiveWrong += 1 }
  }
}

private func orient2dTally(scale s: Double) -> Tally<Orientation> {
  let q = SIMD2(-24.0, -24.0) * s, r = SIMD2(24.0, 24.0) * s
  var t = Tally<Orientation>(coordinates: [q.x, q.y, r.x, r.y])
  for p in ulpGrid(from: SIMD2(0.5, 0.5) * s, size: 64) {
    t.coordinates += [p.x, p.y]
    t.record(orient2d(p, q, r), naive: naiveOrient2d(p, q, r), want: orient2dOracle(p, q, r))
  }
  return t
}

private func inCircleTally(seed: UInt64, centers: ClosedRange<Double>, scale s: Double) -> Tally<CirclePosition> {
  var rng = SplitMix64(seed: seed)
  var t = Tally<CirclePosition>()
  for _ in 0..<64 {
    let (a0, b0, c0, d0) = nearCocircularQuad(&rng, centers: centers)
    let (a, b, c) = (a0 * s, b0 * s, c0 * s)
    t.coordinates += [a.x, a.y, b.x, b.y, c.x, c.y]
    for d in ulpGrid(from: d0 * s, size: 4) {
      t.coordinates += [d.x, d.y]
      t.record(inCircle(a, b, c, d), naive: naiveInCircle(a, b, c, d), want: inCircleOracle(a, b, c, d))
    }
  }
  return t
}

private func orient3dTally(seed: UInt64, range: ClosedRange<Double>, scale s: Double) -> Tally<PlaneSide> {
  var rng = SplitMix64(seed: seed)
  var t = Tally<PlaneSide>()
  for _ in 0..<64 {
    let (a0, b0, c0, d0) = nearCoplanarQuad(&rng, in: range)
    let (a, b, c) = (a0 * s, b0 * s, c0 * s)
    t.coordinates += [a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z]
    for z in ulpSteps(from: d0.z * s, count: 16) {
      let d = SIMD3(d0.x * s, d0.y * s, z)
      t.coordinates += [d.x, d.y, d.z]
      t.record(orient3d(a, b, c, d), naive: naiveOrient3d(a, b, c, d), want: orient3dOracle(a, b, c, d))
    }
  }
  return t
}

// Inside the domain, near-degenerate inputs at each end must match the oracle;
// well outside it, some must not.
@Suite("Input range")
struct InputRangeTests {

  @Test func orient2dIsExactAtTheEdgesOfItsDomain() {
    let domain = Domain(u: 510, l: 459)
    let upper = orient2dTally(scale: powerOfTwo(domain.u - 5))
    let lower = orient2dTally(scale: powerOfTwo(1 - domain.l))
    #expect(domain.contains(upper.coordinates) && domain.reachesUpper(upper.coordinates))
    #expect(domain.contains(lower.coordinates) && domain.reachesLower(lower.coordinates))
    for t in [upper, lower] {
      #expect(t.wrong == 0)
      #expect(t.outcomes == [.ccw, .collinear, .cw])
      #expect(t.naiveWrong > 0)
    }
    #expect(orient2dTally(scale: powerOfTwo(domain.u + 8)).wrong > 0)
    #expect(orient2dTally(scale: powerOfTwo(-domain.l - 200)).wrong > 0)
  }

  @Test func inCircleIsExactAtTheEdgesOfItsDomain() {
    let domain = Domain(u: 254, l: 203)
    let upper = inCircleTally(seed: 191, centers: -1...1, scale: powerOfTwo(domain.u - 2))
    let lower = inCircleTally(seed: 192, centers: 4...5, scale: powerOfTwo(-domain.l - 1))
    #expect(domain.contains(upper.coordinates) && domain.reachesUpper(upper.coordinates))
    #expect(domain.contains(lower.coordinates) && domain.reachesLower(lower.coordinates))
    for t in [upper, lower] {
      #expect(t.wrong == 0)
      #expect(t.outcomes.isSuperset(of: [.inside, .outside]))
      #expect(t.naiveWrong > 0)
    }
    #expect(inCircleTally(seed: 191, centers: -1...1, scale: powerOfTwo(domain.u + 8)).wrong > 0)
    #expect(inCircleTally(seed: 192, centers: 4...5, scale: powerOfTwo(-domain.l - 100)).wrong > 0)
  }

  @Test func orient3dIsExactAtTheEdgesOfItsDomain() {
    let domain = Domain(u: 339, l: 288)
    let upper = orient3dTally(seed: 193, range: -1...1, scale: powerOfTwo(domain.u - 2))
    let lower = orient3dTally(seed: 194, range: 4...5, scale: powerOfTwo(-domain.l - 1))
    #expect(domain.contains(upper.coordinates) && domain.reachesUpper(upper.coordinates))
    #expect(domain.contains(lower.coordinates) && domain.reachesLower(lower.coordinates))
    for t in [upper, lower] {
      #expect(t.wrong == 0)
      #expect(t.outcomes.isSuperset(of: [.above, .below]))
      #expect(t.naiveWrong > 0)
    }
    #expect(orient3dTally(seed: 193, range: -1...1, scale: powerOfTwo(domain.u + 8)).wrong > 0)
    #expect(orient3dTally(seed: 194, range: 4...5, scale: powerOfTwo(-domain.l - 100)).wrong > 0)
  }
}
