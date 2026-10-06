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

// MARK: - Exact

/// The orient2d determinant, computed exactly; only the sign of the result is meaningful.
@usableFromInline
func orient2dExact(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Double {
  let acx = twoDiffE(a.x, c.x)
  let acy = twoDiffE(a.y, c.y)
  let bcx = twoDiffE(b.x, c.x)
  let bcy = twoDiffE(b.y, c.y)
  let left = expansionProduct(acx, bcy)
  let right = expansionProduct(acy, bcx)
  return mostSignificantComponent(expansionSum(left, expansionNegate(right)))
}

// MARK: - Error bounds

/// Shewchuk 1997, §4.3, Table 1.
@inlinable
var ccwErrBoundA: Double { (3 + 16 * shewchukEpsilon) * shewchukEpsilon }

@inlinable
var ccwErrBoundB: Double { (2 + 12 * shewchukEpsilon) * shewchukEpsilon }

@inlinable
var ccwErrBoundC: Double { (9 + 64 * shewchukEpsilon) * shewchukEpsilon * shewchukEpsilon }

// MARK: - Adaptive

/// placeholder
@usableFromInline
func orient2dAdapt(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, detsum: Double) -> Double {

  let acx = a.x - c.x
  let bcx = b.x - c.x
  let acy = a.y - c.y
  let bcy = b.y - c.y

  let (detleft, detlefttail) = twoProd(acx, bcy)
  let (detright, detrighttail) = twoProd(acy, bcx)

  let (b3, b2, b1, b0) = twoTwoDiff(detleft, detlefttail, detright, detrighttail)
  let bExp = [b0, b1, b2, b3]

  var det = estimate(bExp)
  var errBound = ccwErrBoundB * detsum
  if ((det >= errBound) || (-det >= errBound)) {
    return det
  }

  let acxtail = twoDiffTail(a.x, c.x, acx)
  let bcxtail = twoDiffTail(b.x, c.x, bcx)
  let acytail = twoDiffTail(a.y, c.y, acy)
  let bcytail = twoDiffTail(b.y, c.y, bcy)

  if ((acxtail == 0.0) && (acytail == 0.0) && (bcxtail == 0.0) && (bcytail == 0.0)) {
    return det
  }

  errBound = ccwErrBoundC * detsum + resultErrBound * abs(det)
  det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
  if ((det >= errBound) || (-det >= errBound)) {
    return det
  }

  let (s1, s0) = twoProd(acxtail, bcy)
  let (t1, t0) = twoProd(acytail, bcx)
  let (u3, u2, u1, u0) = twoTwoDiff(s1, s0, t1, t0)
  let c1 = expansionSum(bExp, [u0, u1, u2, u3])

  let (s1b, s0b) = twoProd(acx, bcytail)
  let (t1b, t0b) = twoProd(acy, bcxtail)
  let (v3, v2, v1, v0) = twoTwoDiff(s1b, s0b, t1b, t0b)
  let c2 = expansionSum(c1, [v0, v1, v2, v3])

  let (s1c, s0c) = twoProd(acxtail, bcytail)
  let (t1c, t0c) = twoProd(acytail, bcxtail)
  let (w3, w2, w1, w0) = twoTwoDiff(s1c, s0c, t1c, t0c)
  let d = expansionSum(c2, [w0, w1, w2, w3])

  return d[d.count - 1]
}

// MARK: Orient2D Predicate

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
  return Orientation(sign: orient2dAdapt(a, b, c, detsum: detsum))
}
