//
//  Orient3D.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

extension PlaneSide {
  /// The side for an orient3d determinant of sign `s`; positive means below.
  @inlinable init(sign s: Double) {
    self = s > 0 ? .below : (s < 0 ? .above : .on)
  }
}

// MARK: - Exact fallbacks (fully exact from the original coordinates)

/// The orient3d determinant, computed exactly; only the sign of the result is meaningful.
@usableFromInline
func orient3dExact(
  _ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>
) -> Double {
  let adx = twoDiffE(a.x, d.x)
  let ady = twoDiffE(a.y, d.y)
  let adz = twoDiffE(a.z, d.z)
  let bdx = twoDiffE(b.x, d.x)
  let bdy = twoDiffE(b.y, d.y)
  let bdz = twoDiffE(b.z, d.z)
  let cdx = twoDiffE(c.x, d.x)
  let cdy = twoDiffE(c.y, d.y)
  let cdz = twoDiffE(c.z, d.z)
  
  func cross(_ ux: [Double], _ uy: [Double], _ vx: [Double], _ vy: [Double]) -> [Double] {
    expansionSum(expansionProduct(ux, vy), expansionNegate(expansionProduct(uy, vx)))
  }
  
  let aTerm = expansionProduct(adz, cross(bdx, bdy, cdx, cdy))
  let bTerm = expansionProduct(bdz, cross(cdx, cdy, adx, ady))
  let cTerm = expansionProduct(cdz, cross(adx, ady, bdx, bdy))
  return expansionSign(expansionSum(expansionSum(aTerm, bTerm), cTerm))
}

/// Returns which side of the plane through `a`, `b`, `c` the point `d` lies on, exactly.
///
/// The plane is oriented by the normal (b − a) × (c − a); ``PlaneSide/above`` is the
/// side it points to. If `a`, `b`, `c` are collinear, the plane is undefined and every
/// `d` is reported ``PlaneSide/on``. The result is exact for all finite inputs,
/// provided no intermediate value overflows or underflows.
///
/// - Parameters:
///   - a: A point on the plane.
///   - b: A point on the plane.
///   - c: A point on the plane.
///   - d: The point to classify.
/// - Returns: ``PlaneSide/above``, ``PlaneSide/on``, or ``PlaneSide/below``.
@inlinable
public func orient3d(
  _ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>
) -> PlaneSide {
  // The determinant with rows a − d, b − d, c − d. It is positive when d lies
  // opposite the normal, which is Shewchuk's convention.
  let adx = a.x - d.x
  let ady = a.y - d.y
  let adz = a.z - d.z
  let bdx = b.x - d.x
  let bdy = b.y - d.y
  let bdz = b.z - d.z
  let cdx = c.x - d.x
  let cdy = c.y - d.y
  let cdz = c.z - d.z
  
  let bdxcdy = bdx * cdy
  let cdxbdy = cdx * bdy
  let cdxady = cdx * ady
  let adxcdy = adx * cdy
  let adxbdy = adx * bdy
  let bdxady = bdx * ady
  
  let det = adz * (bdxcdy - cdxbdy)
  + bdz * (cdxady - adxcdy)
  + cdz * (adxbdy - bdxady)
  
  // Error bound: Shewchuk 1997, §4.4, Table 3; the permanent is α_a + α_b + α_c.
  let permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz)
  + (abs(cdxady) + abs(adxcdy)) * abs(bdz)
  + (abs(adxbdy) + abs(bdxady)) * abs(cdz)
  let errbound = o3dErrBoundA * permanent
  if det > errbound || -det > errbound {
    return PlaneSide(sign: det)
  }
  return PlaneSide(sign: orient3dExact(a, b, c, d))
}
