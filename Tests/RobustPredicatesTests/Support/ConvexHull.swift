//
//  ConvexHull.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

@testable import RobustPredicates

/// The convex hull's corners in counterclockwise order, starting from the
/// lexicographically smallest point; boundary points that aren't corners are
/// dropped. `orient` is a parameter so the naive predicate can be swapped in.
///
/// Andrew's monotone chain: A. M. Andrew, "Another efficient algorithm for
/// convex hulls in two dimensions", Information Processing Letters 9(5):216–219, 1979.
func convexHull(
  _ points: [SIMD2<Double>],
  orient: (SIMD2<Double>, SIMD2<Double>, SIMD2<Double>) -> Orientation = orient2d
) -> [SIMD2<Double>] {
  let sorted = points.sorted { ($0.x, $0.y) < ($1.x, $1.y) }
  guard sorted.count >= 3 else { return sorted }

  func chain(_ points: some Sequence<SIMD2<Double>>) -> [SIMD2<Double>] {
    var chain: [SIMD2<Double>] = []
    for p in points {
      while chain.count >= 2 && orient(chain[chain.count - 2], chain[chain.count - 1], p) != .ccw {
        chain.removeLast()
      }
      chain.append(p)
    }
    return chain
  }
  // Each chain ends where the other starts.
  return chain(sorted).dropLast() + chain(sorted.reversed()).dropLast()
}
