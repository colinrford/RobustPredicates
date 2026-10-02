//
//  Orient2D.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

extension Orientation {
  /// The orientation for an orient2d determinant of sign `s`; positive means counterclockwise.
  @inlinable
  init(sign s: Double) {
    self = s > 0 ? .ccw : (s < 0 ? .cw : .collinear)
  }
}

// MARK: - Exact fallbacks (fully exact from the original coordinates)

/// The orient2d determinant, computed exactly; only the sign of the result is meaningful.
@usableFromInline
func orient2dExact(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
  let acx = twoDiffE(a.x, c.x), acy = twoDiffE(a.y, c.y)
  let bcx = twoDiffE(b.x, c.x), bcy = twoDiffE(b.y, c.y)
  let left = expansionProduct(acx, bcy)
  let right = expansionProduct(acy, bcx)
  return expansionSign(expansionSum(left, expansionNegate(right)))
}

/// Shewchuk 1997, §4.3, Table 1.
@inlinable
var ccwErrBoundA: Double { (3 + 16 * shewchukEpsilon) * shewchukEpsilon }

/// Returns the orientation of `a`, `b`, `c`, exactly.
///
/// The result is exact when every coordinate is finite and either 0 or of
/// magnitude in [2⁻⁴⁵⁹, 2⁵¹⁰].
///
/// Orientations assume the usual y-up axes.
///
/// - Parameters:
///   - a: The first point.
///   - b: The second point.
///   - c: The point tested against the directed line from `a` to `b`.
/// - Returns: ``Orientation/ccw`` if `c` is left of the line from `a` to `b`,
///   ``Orientation/cw`` if right, ``Orientation/collinear`` if on it.
@inlinable
public func orient2d(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Orientation {
  let detleft = (a.x - c.x) * (b.y - c.y)
  let detright = (a.y - c.y) * (b.x - c.x)
  let det = detleft - detright
  
  // If the products differ in sign or one is zero, the subtraction can't
  // cancel, so the rounded det already has the right sign.
  var detsum: Double
  if detleft > 0 {
    if detright <= 0 { return Orientation(sign: det) }
    detsum = detleft + detright
  } else if detleft < 0 {
    if detright >= 0 { return Orientation(sign: det) }
    detsum = -detleft - detright
  } else {
    return Orientation(sign: det)
  }
  
  // Error bound: Shewchuk 1997, §4.3, Table 1, where detsum is |x₅| + |x₆|.
  let errbound = ccwErrBoundA * detsum
  if det >= errbound || -det >= errbound { return Orientation(sign: det) }
  return Orientation(sign: orient2dExact(a, b, c))
}
