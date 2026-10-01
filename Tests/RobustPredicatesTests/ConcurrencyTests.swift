//
//  ConcurrencyTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

private func requireSendable<T: Sendable>(_: T.Type) {}

private actor OrientationTally {
  private(set) var counts: [Orientation: Int] = [:]
  
  func classify(_ points: [SIMD2<Double>], against q: SIMD2<Double>, _ r: SIMD2<Double>) {
    for p in points { counts[orient2d(p, q, r), default: 0] += 1 }
  }
}

// These compile only if the public API is safe to share across isolation
// domains. The grids are near-degenerate, so the exact fallbacks run
// concurrently too.
@Suite("Concurrency")
struct ConcurrencyTests {
  
  @Test func publicAPIIsSendable() {
    requireSendable(Orientation.self)
    requireSendable(CirclePosition.self)
    requireSendable(PlaneSide.self)
    
    let o2: @Sendable (SIMD2<Double>, SIMD2<Double>, SIMD2<Double>) -> Orientation = orient2d
    let ic: @Sendable (SIMD2<Double>, SIMD2<Double>, SIMD2<Double>, SIMD2<Double>) -> CirclePosition = inCircle
    let o3: @Sendable (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) -> PlaneSide = orient3d
    #expect(o2(SIMD2(0, 0), SIMD2(1, 0), SIMD2(0, 1)) == .ccw)
    #expect(ic(SIMD2(0, 0), SIMD2(1, 0), SIMD2(1, 1), SIMD2(0, 1)) == .on)
    #expect(o3(SIMD3(0, 0, 0), SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)) == .above)
  }
  
  @Test func taskGroupMatchesOracle() async {
    let q = SIMD2(12.0, 12.0), r = SIMD2(24.0, 24.0)
    let grid = ulpGrid(from: SIMD2(0.5, 0.5), size: 128)
    let expected = grid.map { orient2dOracle($0, q, r) }
    
    let chunk = 1024
    let parallel = await withTaskGroup(of: (Int, [Orientation]).self) { group in
      for start in stride(from: 0, to: grid.count, by: chunk) {
        let slice = Array(grid[start..<min(start + chunk, grid.count)])
        group.addTask { (start, slice.map { orient2d($0, q, r) }) }
      }
      var results = [Orientation](repeating: .collinear, count: grid.count)
      for await (start, part) in group {
        results.replaceSubrange(start..<start + part.count, with: part)
      }
      return results
    }
    #expect(parallel == expected)
  }
  
  @Test(arguments: [161 as UInt64])
  func concurrentExactFallbacksAgree(seed: UInt64) async {
    var rng = SplitMix64(seed: seed)
    let quads = (0..<32).map { _ in nearCocircularQuad(&rng) }
    let planes = (0..<32).map { _ in nearCoplanarQuad(&rng) }
    
    let mismatches = await withTaskGroup(of: Int.self) { group in
      for (a, b, c, d) in quads {
        group.addTask {
          ulpGrid(from: d, size: 8).count { inCircle(a, b, c, $0) != inCircleOracle(a, b, c, $0) }
        }
      }
      for (a, b, c, d) in planes {
        group.addTask {
          ulpSteps(from: d.z, count: 64).count {
            let p = SIMD3(d.x, d.y, $0)
            return orient3d(a, b, c, p) != orient3dOracle(a, b, c, p)
          }
        }
      }
      return await group.reduce(0, +)
    }
    #expect(mismatches == 0)
  }
  
  @Test func actorAndMainActorCallersAgree() async {
    let q = SIMD2(12.0, 12.0), r = SIMD2(24.0, 24.0)
    let grid = ulpGrid(from: SIMD2(0.5, 0.5), size: 64)
    
    let tally = OrientationTally()
    await tally.classify(grid, against: q, r)
    let fromActor = await tally.counts
    
    let fromMainActor = await MainActor.run {
      Dictionary(grouping: grid) { orient2d($0, q, r) }.mapValues(\.count)
    }
    let expected = Dictionary(grouping: grid) { orient2dOracle($0, q, r) }.mapValues(\.count)
    #expect(fromActor == expected)
    #expect(fromMainActor == expected)
    #expect(Set(expected.keys) == [.ccw, .collinear, .cw])
  }
}
