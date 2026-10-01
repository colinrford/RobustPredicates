//
//  Triangulation.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

/// Points and the counterclockwise index triples of the triangles over them.
struct Triangulation {
  var points: [SIMD2<Double>]
  var triangles: [SIMD3<Int>]

  /// An edge `u`–`v` shared by two triangles: `(u, v, w)` on its left and
  /// `(v, u, x)` on its right.
  struct InteriorEdge: Sendable {
    let u, v, w, x: Int
  }

  /// Every edge shared by two triangles, once each.
  var interiorEdges: [InteriorEdge] {
    var apex: [SIMD2<Int>: Int] = [:]
    for t in triangles {
      apex[SIMD2(t[0], t[1])] = t[2]
      apex[SIMD2(t[1], t[2])] = t[0]
      apex[SIMD2(t[2], t[0])] = t[1]
    }
    return apex.compactMap { edge, w in
      guard edge[0] < edge[1], let x = apex[SIMD2(edge[1], edge[0])] else { return nil }
      return InteriorEdge(u: edge[0], v: edge[1], w: w, x: x)
    }
    .sorted { ($0.u, $0.v) < ($1.u, $1.v) }
  }
}
