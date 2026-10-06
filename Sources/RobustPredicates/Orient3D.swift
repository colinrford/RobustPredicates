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

// MARK: - Exact fallbacks

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
  return mostSignificantComponent(expansionSum(expansionSum(aTerm, bTerm), cTerm))
}

// MARK: - Error bounds
/// Shewchuk 1997, §4.4, Table 3.
@inlinable
var o3dErrBoundA: Double { (7 + 56 * shewchukEpsilon) * shewchukEpsilon }

@inlinable
var o3dErrBoundB: Double { (3 + 28 * shewchukEpsilon) * shewchukEpsilon }

@inlinable
var o3dErrBoundC: Double { (26 + 288 * shewchukEpsilon) * shewchukEpsilon * shewchukEpsilon }

// MARK: - Adaptive
@usableFromInline
func orient3dAdapt(_ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>, permanent: Double) -> Double {

  let adx = a.x - d.x
  let bdx = b.x - d.x
  let cdx = c.x - d.x
  let ady = a.y - d.y
  let bdy = b.y - d.y
  let cdy = c.y - d.y
  let adz = a.z - d.z
  let bdz = b.z - d.z
  let cdz = c.z - d.z

  let (bdxcdy1, bdxcdy0) = twoProd(bdx, cdy)
  let (cdxbdy1, cdxbdy0) = twoProd(cdx, bdy)
  let (bc3, bc2, bc1, bc0) = twoTwoDiff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0)
  let bc = [bc0, bc1, bc2, bc3]
  let adet = expansionScale(bc, adz)

  let (cdxady1, cdxady0) = twoProd(cdx, ady)
  let (adxcdy1, adxcdy0) = twoProd(adx, cdy)
  let (ca3, ca2, ca1, ca0) = twoTwoDiff(cdxady1, cdxady0, adxcdy1, adxcdy0)
  let ca = [ca0, ca1, ca2, ca3]
  let bdet = expansionScale(ca, bdz)

  let (adxbdy1, adxbdy0) = twoProd(adx, bdy)
  let (bdxady1, bdxady0) = twoProd(bdx, ady)
  let (ab3, ab2, ab1, ab0) = twoTwoDiff(adxbdy1, adxbdy0, bdxady1, bdxady0)
  let ab = [ab0, ab1, ab2, ab3]
  let cdet = expansionScale(ab, cdz)

  let abdet = expansionSum(adet, bdet)
  var fin1 = expansionSum(abdet, cdet)

  var det = estimate(fin1)
  var errBound = o3dErrBoundB * permanent
  if ((det >= errBound) || (-det >= errBound)) {
    return det
  }

  let adxtail = twoDiffTail(a.x, d.x, adx)
  let bdxtail = twoDiffTail(b.x, d.x, bdx)
  let cdxtail = twoDiffTail(c.x, d.x, cdx)
  let adytail = twoDiffTail(a.y, d.y, ady)
  let bdytail = twoDiffTail(b.y, d.y, bdy)
  let cdytail = twoDiffTail(c.y, d.y, cdy)
  let adztail = twoDiffTail(a.z, d.z, adz)
  let bdztail = twoDiffTail(b.z, d.z, bdz)
  let cdztail = twoDiffTail(c.z, d.z, cdz)

  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
      && (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
    return det
  }

  errBound = o3dErrBoundC * permanent + resultErrBound * abs(det)
  det += (adz * ((bdx * cdytail + cdy * bdxtail)
                 - (bdy * cdxtail + cdx * bdytail))
          + adztail * (bdx * cdy - bdy * cdx))
       + (bdz * ((cdx * adytail + ady * cdxtail)
                 - (cdy * adxtail + adx * cdytail))
          + bdztail * (cdx * ady - cdy * adx))
       + (cdz * ((adx * bdytail + bdy * adxtail)
                 - (ady * bdxtail + bdx * adytail))
          + cdztail * (adx * bdy - ady * bdx))
  
  if ((det >= errBound) || (-det >= errBound)) {
    return det
  }

  let at_b: [Double]
  let at_c: [Double]
  if (adxtail == 0.0) {
    if (adytail == 0.0) {
      at_b = [0.0]
      at_c = [0.0]
    } else {
      let (at_blarge, at_b0) = twoProd(-adytail, bdx)
      at_b = [at_b0, at_blarge]
      let (at_clarge, at_c0) = twoProd(adytail, cdx)
      at_c = [at_c0, at_clarge]
    }
  } else {
    if (adytail == 0.0) {
      let (at_blarge, at_b0) = twoProd(adxtail, bdy)
      at_b = [at_b0, at_blarge]
      let (at_clarge, at_c0) = twoProd(-adxtail, cdy)
      at_c = [at_c0, at_clarge]
    } else {
      let (adxt_bdy1, adxt_bdy0) = twoProd(adxtail, bdy)
      let (adyt_bdx1, adyt_bdx0) = twoProd(adytail, bdx)
      let (at_blarge, at_b2, at_b1, at_b0) = twoTwoDiff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0)
      at_b = [at_b0, at_b1, at_b2, at_blarge]
      
      let (adyt_cdx1, adyt_cdx0) = twoProd(adytail, cdx)
      let (adxt_cdy1, adxt_cdy0) = twoProd(adxtail, cdy)
      let (at_clarge, at_c2, at_c1, at_c0) = twoTwoDiff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0)
      at_c = [at_c0, at_c1, at_c2, at_clarge]
    }
  }
  
  let bt_c: [Double]
  let bt_a: [Double]
  if (bdxtail == 0.0) {
    if (bdytail == 0.0) {
      bt_c = [0.0]
      bt_a = [0.0]
    } else {
      let (bt_clarge, bt_c0) = twoProd(-bdytail, cdx)
      bt_c = [bt_c0, bt_clarge]
      let (bt_alarge, bt_a0) = twoProd(bdytail, adx)
      bt_a = [bt_a0, bt_alarge]
    }
  } else {
    if (bdytail == 0.0) {
      let (bt_clarge, bt_c0) = twoProd(bdxtail, cdy)
      bt_c = [bt_c0, bt_clarge]
      
      let (bt_alarge, bt_a0) = twoProd(-bdxtail, ady)
      bt_a = [bt_a0, bt_alarge]
    } else {
      let (bdxt_cdy1, bdxt_cdy0) = twoProd(bdxtail, cdy)
      let (bdyt_cdx1, bdyt_cdx0) = twoProd(bdytail, cdx)
      let (bt_clarge, bt_c2, bt_c1, bt_c0) = twoTwoDiff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0)
      bt_c = [bt_c0, bt_c1, bt_c2, bt_clarge]
      
      let (bdyt_adx1, bdyt_adx0) = twoProd(bdytail, adx)
      let (bdxt_ady1, bdxt_ady0) = twoProd(bdxtail, ady)
      let (bt_alarge, bt_a2, bt_a1, bt_a0) = twoTwoDiff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0)
      bt_a = [bt_a0, bt_a1, bt_a2, bt_alarge]
    }
  }
  
  let ct_a: [Double]
  let ct_b: [Double]
  if (cdxtail == 0.0) {
    if (cdytail == 0.0) {
      ct_a = [0.0]
      ct_b = [0.0]
    } else {
      let (ct_alarge, ct_a0) = twoProd(-cdytail, adx)
      ct_a = [ct_a0, ct_alarge]
      let (ct_blarge, ct_b0) = twoProd(cdytail, bdx)
      ct_b = [ct_b0, ct_blarge]
    }
  } else {
    if (cdytail == 0.0) {
      let (ct_alarge, ct_a0) = twoProd(cdxtail, ady)
      ct_a = [ct_a0, ct_alarge]
      let (ct_blarge, ct_b0) = twoProd(-cdxtail, bdy)
      ct_b = [ct_b0, ct_blarge]
    } else {
      let (cdxt_ady1, cdxt_ady0) = twoProd(cdxtail, ady)
      let (cdyt_adx1, cdyt_adx0) = twoProd(cdytail, adx)
      let (ct_alarge, ct_a2, ct_a1, ct_a0) = twoTwoDiff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0)
      ct_a = [ct_a0, ct_a1, ct_a2, ct_alarge]
      
      let (cdyt_bdx1, cdyt_bdx0) = twoProd(cdytail, bdx)
      let (cdxt_bdy1, cdxt_bdy0) = twoProd(cdxtail, bdy)
      let (ct_blarge, ct_b2, ct_b1, ct_b0) = twoTwoDiff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0)
      ct_b = [ct_b0, ct_b1, ct_b2, ct_blarge]
    }
  }

  let bct = expansionSum(bt_c, ct_b)
  let wa = expansionScale(bct, adz)
  fin1 = expansionSum(fin1, wa)

  let cat = expansionSum(ct_a, at_c)
  let wb = expansionScale(cat, bdz)
  fin1 = expansionSum(fin1, wb)

  let abt = expansionSum(at_b, bt_a)
  let wc = expansionScale(abt, cdz)
  fin1 = expansionSum(fin1, wc)

  if (adztail != 0.0) {
    let v = expansionScale(bc, adztail)
    fin1 = expansionSum(fin1, v)
  }
  if (bdztail != 0.0) {
    let v = expansionScale(ca, bdztail)
    fin1 = expansionSum(fin1, v)
  }
  if (cdztail != 0.0) {
    let v = expansionScale(ab, cdztail)
    fin1 = expansionSum(fin1, v)
  }

  if (adxtail != 0.0) {
    if (bdytail != 0.0) {
      let (adxt_bdyt1, adxt_bdyt0) = twoProd(adxtail, bdytail)
      let (u3, u2, u1, u0) = twoOneProd(adxt_bdyt1, adxt_bdyt0, cdz)
      let u = [u0, u1, u2, u3]
      fin1 = expansionSum(fin1, u)
      
      if (cdztail != 0.0) {
        let (u3, u2, u1, u0) = twoOneProd(adxt_bdyt1, adxt_bdyt0, cdztail)
        let u = [u0, u1, u2, u3]
        fin1 = expansionSum(fin1, u)
      }
    }
    if (cdytail != 0.0) {
      let (adxt_cdyt1, adxt_cdyt0) = twoProd(-adxtail, cdytail)
      let (u3, u2, u1, u0) = twoOneProd(adxt_cdyt1, adxt_cdyt0, bdz)
      let u = [u0, u1, u2, u3]
      fin1 = expansionSum(fin1, u)
      
      if (bdztail != 0.0) {
        let (u3, u2, u1, u0) = twoOneProd(adxt_cdyt1, adxt_cdyt0, bdztail)
        let u = [u0, u1, u2, u3]
        fin1 = expansionSum(fin1, u)
      }
    }
  }
  if (bdxtail != 0.0) {
    if (cdytail != 0.0) {
      let (bdxt_cdyt1, bdxt_cdyt0) = twoProd(bdxtail, cdytail)
      let (u3, u2, u1, u0) = twoOneProd(bdxt_cdyt1, bdxt_cdyt0, adz)
      let u  = [u0, u1, u2, u3]
      fin1 = expansionSum(fin1, u)
      
      if (adztail != 0.0) {
        let (u3, u2, u1, u0) = twoOneProd(bdxt_cdyt1, bdxt_cdyt0, adztail)
        let u = [u0, u1, u2, u3]
        fin1 = expansionSum(fin1, u)
      }
    }
    if (adytail != 0.0) {
      let (bdxt_adyt1, bdxt_adyt0) = twoProd(-bdxtail, adytail)
      let (u3, u2, u1, u0) = twoOneProd(bdxt_adyt1, bdxt_adyt0, cdz)
      let u = [u0, u1, u2, u3]
      fin1 = expansionSum(fin1, u)
      
      if (cdztail != 0.0) {
        let (u3, u2, u1, u0) = twoOneProd(bdxt_adyt1, bdxt_adyt0, cdztail)
        let u = [u0, u1, u2, u3]
        fin1 = expansionSum(fin1, u)
      }
    }
  }
  if (cdxtail != 0.0) {
    if (adytail != 0.0) {
      let (cdxt_adyt1, cdxt_adyt0) = twoProd(cdxtail, adytail)
      let (u3, u2, u1, u0) = twoOneProd(cdxt_adyt1, cdxt_adyt0, bdz)
      let u  = [u0, u1, u2, u3]
      fin1 = expansionSum(fin1, u)
      
      if (bdztail != 0.0) {
        let (u3, u2, u1, u0) = twoOneProd(cdxt_adyt1, cdxt_adyt0, bdztail)
        let u = [u0, u1, u2, u3]
        fin1 = expansionSum(fin1, u)
      }
    }
    if (bdytail != 0.0) {
      let (cdxt_bdyt1, cdxt_bdyt0) = twoProd(-cdxtail, bdytail)
      let (u3, u2, u1, u0) = twoOneProd(cdxt_bdyt1, cdxt_bdyt0, adz)
      let u = [u0, u1, u2, u3]
      fin1 = expansionSum(fin1, u)
      
      if (adztail != 0.0) {
        let (u3, u2, u1, u0) = twoOneProd(cdxt_bdyt1, cdxt_bdyt0, adztail)
        let u = [u0, u1, u2, u3]
        fin1 = expansionSum(fin1, u)
      }
    }
  }

  if (adztail != 0.0) {
    let w = expansionScale(bct, adztail)
    fin1 = expansionSum(fin1, w)
  }
  
  if (bdztail != 0.0) {
    let w = expansionScale(cat, bdztail)
    fin1 = expansionSum(fin1, w)
  }
  
  if (cdztail != 0.0) {
    let w = expansionScale(abt, cdztail)
    fin1 = expansionSum(fin1, w)
  }

  return fin1[fin1.count - 1]
}

// MARK: Orient3D Predicate

/// Returns which side of the plane through `a`, `b`, `c` the point `d` lies on, exactly.
///
/// The plane is oriented by the normal (b − a) × (c − a); ``PlaneSide/above`` is the
/// side it points to. If `a`, `b`, `c` are collinear, the plane is undefined and every
/// `d` is reported ``PlaneSide/on``.
///
/// The result is exact when every coordinate is finite and either 0 or of
/// magnitude in [2⁻²⁸⁸, 2³³⁹].
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
  return PlaneSide(sign: orient3dAdapt(a, b, c, d, permanent: permanent))
}
