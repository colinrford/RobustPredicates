//
//  WorkloadTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

// MARK: - Convex hull

/// The corners of a random convex quadrilateral, `perSide` rounded points on
/// each side, and `interior` points well inside. The rounded side points land
/// within an ulp or so of the boundary, on either side of it.
private func nearDegenerateHullInput(
  _ rng: inout SplitMix64, perSide: Int, interior: Int
) -> [SIMD2<Double>] {
  let jitter = { (rng: inout SplitMix64) in randomPoint(&rng, in: -0.25...0.25) }
  let corners = [SIMD2(-1.0, -1.0), SIMD2(1.0, -1.0), SIMD2(1.0, 1.0), SIMD2(-1.0, 1.0)]
    .map { $0 + jitter(&rng) }
  var points = corners
  for (a, b) in zip(corners, corners.dropFirst() + [corners[0]]) {
    for _ in 0..<perSide {
      points.append(a + (b - a) * Double.random(in: 0...1, using: &rng))
    }
  }
  for _ in 0..<interior { points.append(randomPoint(&rng, in: -0.5...0.5)) }
  return points
}

/// The hull of each vertical strip in parallel, then the hull of their corners.
/// This is the same hull, since every corner of the whole hull is a corner of
/// its strip's hull.
private func parallelConvexHull(
  _ points: [SIMD2<Double>], strips: Int, probe: OverlapProbe
) async -> [SIMD2<Double>] {
  let sorted = points.sorted { ($0.x, $0.y) < ($1.x, $1.y) }
  let size = (sorted.count + strips - 1) / strips
  let corners = await withTaskGroup(of: [SIMD2<Double>].self) { group in
    for start in stride(from: 0, to: sorted.count, by: size) {
      let strip = Array(sorted[start..<min(start + size, sorted.count)])
      group.addTask { probe.measure { convexHull(strip) } }
    }
    return await group.reduce(into: []) { $0 += $1 }
  }
  return convexHull(corners)
}

/// Whether `hull` is a strictly convex counterclockwise polygon with corners
/// from `points` that contains them all, decided by the exact oracle.
private func isExactHull(_ hull: [SIMD2<Double>], of points: [SIMD2<Double>]) -> Bool {
  guard hull.count >= 3, Set(hull).isSubset(of: points) else { return false }
  let next = { (i: Int) in hull[(i + 1) % hull.count] }
  let convex = hull.indices.allSatisfy { orient2dOracle(hull[$0], next($0), next($0 + 1)) == .ccw }
  return convex && points.allSatisfy { p in
    hull.indices.allSatisfy { orient2dOracle(hull[$0], next($0), p) != .cw }
  }
}

// MARK: - Delaunay check

/// A `(2k + 1)²`-point integer grid rotated by the angle with cosine `c / 5` and
/// sine `s / 5`, each square split along one diagonal. Exactly rotated, every
/// square is cocircular; unless the rotation is by a multiple of 90°, rounding
/// the coordinates leaves them only nearly so.
private func rotatedGrid(k: Int, cos c: Int, sin s: Int) -> Triangulation {
  let n = 2 * k + 1
  var points: [SIMD2<Double>] = []
  for j in -k...k {
    for i in -k...k {
      points.append(SIMD2(Double(c * i - s * j) / 5, Double(s * i + c * j) / 5))
    }
  }
  var triangles: [SIMD3<Int>] = []
  for j in 0..<n - 1 {
    for i in 0..<n - 1 {
      let p = j * n + i
      triangles.append(SIMD3(p, p + 1, p + n + 1))
      triangles.append(SIMD3(p, p + n + 1, p + n))
    }
  }
  return Triangulation(points: points, triangles: triangles)
}

/// The indices of edges whose far vertex `x` lies strictly inside the circle
/// through `u`, `v`, `w`, checked in parallel chunks.
///
/// A triangulation is Delaunay exactly when no edge fails this local test:
/// B. Delaunay, "Sur la sphère vide", Bull. Acad. Sci. URSS, Classe Sci. Mat.
/// Nat. 6:793–800, 1934. The edge test is Lawson's: C. L. Lawson, "Software
/// for C¹ surface interpolation", Mathematical Software III, pp. 161–194, 1977.
private func illegalEdges(
  _ t: Triangulation, _ edges: [Triangulation.InteriorEdge], chunks: Int, probe: OverlapProbe
) async -> Set<Int> {
  let size = (edges.count + chunks - 1) / chunks
  let points = t.points
  return await withTaskGroup(of: [Int].self) { group in
    for start in stride(from: 0, to: edges.count, by: size) {
      let range = start..<min(start + size, edges.count)
      group.addTask {
        probe.measure {
          range.filter { i in
            let e = edges[i]
            return inCircle(points[e.u], points[e.v], points[e.w], points[e.x]) == .inside
          }
        }
      }
    }
    return await group.reduce(into: []) { $0.formUnion($1) }
  }
}

@Suite("Parallel workloads")
struct WorkloadTests {
  
  @Test(arguments: [171 as UInt64, 172])
  func parallelHullIsExact(seed: UInt64) async {
    var rng = SplitMix64(seed: seed)
    let points = nearDegenerateHullInput(&rng, perSide: 250, interior: 100)
    let probe = OverlapProbe()
    let hull = await parallelConvexHull(points, strips: 8, probe: probe)
    #expect(hull == convexHull(points))
    #expect(isExactHull(hull, of: points))
    #expect(hull.count > 4, "some rounded side points must be corners for this to test anything")
    #expect(!isExactHull(convexHull(points, orient: naiveOrient2d), of: points))
    #expect(probe.peak >= 2, "the strips must have run at the same time")
  }
  
  @Test func parallelDelaunayCheckIsExact() async {
    let t = rotatedGrid(k: 16, cos: 4, sin: 3)
    let edges = t.interiorEdges
    let probe = OverlapProbe()
    let illegal = await illegalEdges(t, edges, chunks: 8, probe: probe)
    
    let p = t.points
    let oracle = edges.map { inCircleOracle(p[$0.u], p[$0.v], p[$0.w], p[$0.x]) }
    #expect(illegal == Set(oracle.indices.filter { oracle[$0] == .inside }))
    #expect(Set(oracle) == [.inside, .on, .outside])
    #expect(edges.map { naiveInCircle(p[$0.u], p[$0.v], p[$0.w], p[$0.x]) } != oracle)
    #expect(probe.peak >= 2, "the chunks must have run at the same time")
  }
  
  // Unrotated, every diagonal's four points are exactly cocircular, so the
  // triangulation is Delaunay, just not uniquely.
  @Test func unrotatedGridIsDelaunay() async {
    let k = 16
    let t = rotatedGrid(k: k, cos: 5, sin: 0)
    let edges = t.interiorEdges
    let illegal = await illegalEdges(t, edges, chunks: 8, probe: OverlapProbe())
    #expect(illegal.isEmpty)
    
    let p = t.points
    let onCircle = edges.count { inCircle(p[$0.u], p[$0.v], p[$0.w], p[$0.x]) == .on }
    #expect(onCircle == (2 * k) * (2 * k))
  }
}
