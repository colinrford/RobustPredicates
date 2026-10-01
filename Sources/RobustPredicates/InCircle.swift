//
//  InCircle.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

extension CirclePosition {
  /// The position for an incircle determinant of sign `s`; positive means inside.
  @inlinable init(sign s: Double) {
    self = s > 0 ? .inside : (s < 0 ? .outside : .on)
  }
}

// MARK: - Exact fallbacks (fully exact from the original coordinates)

/// The incircle determinant, computed exactly; only the sign of the result is meaningful.
@usableFromInline
func inCircleExact(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>
) -> Double {
  let adx = twoDiffE(a.x, d.x), ady = twoDiffE(a.y, d.y)
  let bdx = twoDiffE(b.x, d.x), bdy = twoDiffE(b.y, d.y)
  let cdx = twoDiffE(c.x, d.x), cdy = twoDiffE(c.y, d.y)

  func cross(_ ux: [Double], _ uy: [Double], _ vx: [Double], _ vy: [Double]) -> [Double] {
    expansionSum(expansionProduct(ux, vy), expansionNegate(expansionProduct(uy, vx)))
  }
  func lift(_ x: [Double], _ y: [Double]) -> [Double] {
    expansionSum(expansionProduct(x, x), expansionProduct(y, y))
  }
  
  let aTerm = expansionProduct(lift(adx, ady), cross(bdx, bdy, cdx, cdy))
  let bTerm = expansionProduct(lift(bdx, bdy), cross(cdx, cdy, adx, ady))
  let cTerm = expansionProduct(lift(cdx, cdy), cross(adx, ady, bdx, bdy))
  return expansionSign(expansionSum(expansionSum(aTerm, bTerm), cTerm))
}

/// Shewchuk 1997, §4.4, Table 5.
@inlinable
var iccErrBoundA: Double { (10 + 96 * shewchukEpsilon) * shewchukEpsilon }

/// Returns where a point lies relative to the circle through three others, exactly.
///
/// The circle passes through `a`, `b`, `c`, which must not be collinear. List them
/// counterclockwise; for a clockwise triangle, ``CirclePosition/inside`` and
/// ``CirclePosition/outside`` swap. The result is exact for all finite inputs,
/// provided no intermediate value overflows or underflows.
///
/// - Parameters:
///   - a: A point on the circle.
///   - b: A point on the circle.
///   - c: A point on the circle.
///   - d: The point to classify.
/// - Returns: ``CirclePosition/inside``, ``CirclePosition/on``, or ``CirclePosition/outside``.
@inlinable
public func inCircle(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>
) -> CirclePosition {
  let adx = a.x - d.x, ady = a.y - d.y
  let bdx = b.x - d.x, bdy = b.y - d.y
  let cdx = c.x - d.x, cdy = c.y - d.y
  
  let bdxcdy = bdx * cdy, cdxbdy = cdx * bdy
  let alift = adx * adx + ady * ady
  let cdxady = cdx * ady, adxcdy = adx * cdy
  let blift = bdx * bdx + bdy * bdy
  let adxbdy = adx * bdy, bdxady = bdx * ady
  let clift = cdx * cdx + cdy * cdy
  
  let det = alift * (bdxcdy - cdxbdy)
  + blift * (cdxady - adxcdy)
  + clift * (adxbdy - bdxady)
  
  // Error bound: Shewchuk 1997, §4.4, Table 5; the permanent is α_a + α_b + α_c.
  let permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift
  + (abs(cdxady) + abs(adxcdy)) * blift
  + (abs(adxbdy) + abs(bdxady)) * clift
  let errbound = iccErrBoundA * permanent
  if det > errbound || -det > errbound {
    return CirclePosition(sign: det)
  }
  return CirclePosition(sign: inCircleExact(a, b, c, d))
}
