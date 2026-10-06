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

// MARK: - Exact

/// The incircle determinant, computed exactly; only the sign of the result is meaningful.
@usableFromInline
func inCircleExact(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>
) -> Double {
  let adx = twoDiffE(a.x, d.x)
  let ady = twoDiffE(a.y, d.y)
  let bdx = twoDiffE(b.x, d.x)
  let bdy = twoDiffE(b.y, d.y)
  let cdx = twoDiffE(c.x, d.x)
  let cdy = twoDiffE(c.y, d.y)
  
  func cross(_ ux: [Double], _ uy: [Double], _ vx: [Double], _ vy: [Double]) -> [Double] {
    expansionSum(expansionProduct(ux, vy), expansionNegate(expansionProduct(uy, vx)))
  }
  func lift(_ x: [Double], _ y: [Double]) -> [Double] {
    expansionSum(expansionProduct(x, x), expansionProduct(y, y))
  }
  
  let aTerm = expansionProduct(lift(adx, ady), cross(bdx, bdy, cdx, cdy))
  let bTerm = expansionProduct(lift(bdx, bdy), cross(cdx, cdy, adx, ady))
  let cTerm = expansionProduct(lift(cdx, cdy), cross(adx, ady, bdx, bdy))
  return mostSignificantComponent(expansionSum(expansionSum(aTerm, bTerm), cTerm))
}

// MARK: - Error Bounds

/// Shewchuk 1997, §4.4, Table 5.
@inlinable
var iccErrBoundA: Double { (10 + 96 * shewchukEpsilon) * shewchukEpsilon }

@inlinable
var iccErrBoundB: Double { (4 + 48 * shewchukEpsilon) * shewchukEpsilon }

@inlinable
var iccErrBoundC: Double { (44 + 576 * shewchukEpsilon) * shewchukEpsilon * shewchukEpsilon }

// MARK: - Adaptive
@usableFromInline
func inCircleAdapt(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>, permanent: Double) -> Double {
  
  let adx = a.x - d.x
  let bdx = b.x - d.x
  let cdx = c.x - d.x
  let ady = a.y - d.y
  let bdy = b.y - d.y
  let cdy = c.y - d.y
  
  let (bdxcdy1, bdxcdy0) = twoProd(bdx, cdy)
  let (cdxbdy1, cdxbdy0) = twoProd(cdx, bdy)
  let (bc3, bc2, bc1, bc0) = twoTwoDiff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0)
  let bc = [bc0, bc1, bc2, bc3]
  let axbc = expansionScale(bc, adx)
  let axxbc = expansionScale(axbc, adx)
  let aybc = expansionScale(bc, ady)
  let ayybc = expansionScale(aybc, ady)
  let adet = expansionSum(axxbc, ayybc)
  
  let (cdxady1, cdxady0) = twoProd(cdx, ady)
  let (adxcdy1, adxcdy0) = twoProd(adx, cdy)
  let (ca3, ca2, ca1, ca0) = twoTwoDiff(cdxady1, cdxady0, adxcdy1, adxcdy0)
  let ca = [ca0, ca1, ca2, ca3]
  let bxca = expansionScale(ca, bdx)
  let bxxca = expansionScale(bxca, bdx)
  let byca = expansionScale(ca, bdy)
  let byyca = expansionScale(byca, bdy)
  let bdet = expansionSum(bxxca, byyca)
  
  let (adxbdy1, adxbdy0) = twoProd(adx, bdy)
  let (bdxady1, bdxady0) = twoProd(bdx, ady)
  let (ab3, ab2, ab1, ab0) = twoTwoDiff(adxbdy1, adxbdy0, bdxady1, bdxady0)
  let ab = [ab0, ab1, ab2, ab3]
  let cxab = expansionScale(ab, cdx)
  let cxxab = expansionScale(cxab, cdx)
  let cyab = expansionScale(ab, cdy)
  let cyyab = expansionScale(cyab, cdy)
  let cdet = expansionSum(cxxab, cyyab)
  
  let abdet = expansionSum(adet, bdet)
  var fin1 = expansionSum(abdet, cdet)
  
  var det = estimate(fin1)
  var errBound = iccErrBoundB * permanent
  if ((det >= errBound) || (-det >= errBound)) {
    return det
  }
  
  let adxtail = twoDiffTail(a.x, d.x, adx)
  let adytail = twoDiffTail(a.y, d.y, ady)
  let bdxtail = twoDiffTail(b.x, d.x, bdx)
  let bdytail = twoDiffTail(b.y, d.y, bdy)
  let cdxtail = twoDiffTail(c.x, d.x, cdx)
  let cdytail = twoDiffTail(c.y, d.y, cdy)
  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) {
    return det
  }
  
  errBound = iccErrBoundC * permanent + resultErrBound * abs(det)
  det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
                                     - (bdy * cdxtail + cdx * bdytail))
          + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
  + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
                                - (cdy * adxtail + adx * cdytail))
     + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
  + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
                                - (ady * bdxtail + bdx * adytail))
     + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx))
  
  if ((det >= errBound) || (-det >= errBound)) {
    return det
  }
  
  var aa: [Double] = []
  if ((bdxtail != 0.0) || (bdytail != 0.0)
      || (cdxtail != 0.0) || (cdytail != 0.0)) {
    let (adxadx1, adxadx0) = twoProd(adx, adx)
    let (adyady1, adyady0) = twoProd(ady, ady)
    let (aa3, aa2, aa1, aa0) = twoTwoSum(adxadx1, adxadx0, adyady1, adyady0)
    aa = [aa0, aa1, aa2, aa3]
  }
  
  var bb: [Double] = []
  if ((cdxtail != 0.0) || (cdytail != 0.0)
      || (adxtail != 0.0) || (adytail != 0.0)) {
    let (bdxbdx1, bdxbdx0) = twoProd(bdx, bdx)
    let (bdybdy1, bdybdy0) = twoProd(bdy, bdy)
    let (bb3, bb2, bb1, bb0) = twoTwoSum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0)
    bb = [bb0, bb1, bb2, bb3]
  }
  
  var cc: [Double] = []
  if ((adxtail != 0.0) || (adytail != 0.0)
      || (bdxtail != 0.0) || (bdytail != 0.0)) {
    let (cdxcdx1, cdxcdx0) = twoProd(cdx, cdx)
    let (cdycdy1, cdycdy0) = twoProd(cdy, cdy)
    let (cc3, cc2, cc1, cc0) = twoTwoSum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0)
    cc = [cc0, cc1, cc2, cc3]
  }
  
  var axtbc: [Double] = []
  var aytbc: [Double] = []
  var bxtca: [Double] = []
  var bytca: [Double] = []
  var cxtab: [Double] = []
  var cytab: [Double] = []
  
  if (adxtail != 0.0) {
    axtbc = expansionScale(bc, adxtail)
    let temp16a = expansionScale(axtbc, 2.0 * adx)
    
    let axtcc = expansionScale(cc, adxtail)
    let temp16b = expansionScale(axtcc, bdy)
    
    let axtbb = expansionScale(bb, adxtail)
    let temp16c = expansionScale(axtbb, -cdy)
    
    let temp32a = expansionSum(temp16a, temp16b)
    let temp48 = expansionSum(temp16c, temp32a)
    fin1 = expansionSum(fin1, temp48)
  }
  if (adytail != 0.0) {
    aytbc = expansionScale(bc, adytail)
    let temp16a = expansionScale(aytbc, 2.0 * ady)
    
    let aytbb = expansionScale(bb, adytail)
    let temp16b = expansionScale(aytbb, cdx)
    
    let aytcc = expansionScale(cc, adytail)
    let temp16c = expansionScale(aytcc, -bdx)
    
    let temp32a = expansionSum(temp16a, temp16b)
    let temp48 = expansionSum(temp16c, temp32a)
    fin1 = expansionSum(fin1, temp48)
  }
  if (bdxtail != 0.0) {
    bxtca = expansionScale(ca, bdxtail)
    let temp16a = expansionScale(bxtca, 2.0 * bdx)
    
    let bxtaa = expansionScale(aa, bdxtail)
    let temp16b = expansionScale(bxtaa, cdy)
    
    let bxtcc = expansionScale(cc, bdxtail)
    let temp16c = expansionScale(bxtcc, -ady)
    
    let temp32a = expansionSum(temp16a, temp16b)
    let temp48 = expansionSum(temp16c, temp32a)
    fin1 = expansionSum(fin1, temp48)
  }
  if (bdytail != 0.0) {
    bytca = expansionScale(ca, bdytail)
    let temp16a = expansionScale(bytca, 2.0 * bdy)
    
    let bytcc = expansionScale(cc, bdytail)
    let temp16b = expansionScale(bytcc, adx)
    
    let bytaa = expansionScale(aa, bdytail)
    let temp16c = expansionScale(bytaa, -cdx)
    
    let temp32a = expansionSum(temp16a, temp16b)
    let temp48 = expansionSum(temp16c, temp32a)
    fin1 = expansionSum(fin1, temp48)
  }
  if (cdxtail != 0.0) {
    cxtab = expansionScale(ab, cdxtail)
    let temp16a = expansionScale(cxtab, 2.0 * cdx)
    
    let cxtbb = expansionScale(bb, cdxtail)
    let temp16b = expansionScale(cxtbb, ady)
    
    let cxtaa = expansionScale(aa, cdxtail)
    let temp16c = expansionScale(cxtaa, -bdy)
    
    let temp32a = expansionSum(temp16a, temp16b)
    let temp48 = expansionSum(temp16c, temp32a)
    fin1 = expansionSum(fin1, temp48)
  }
  if (cdytail != 0.0) {
    cytab = expansionScale(ab, cdytail)
    let temp16a = expansionScale(cytab, 2.0 * cdy)
    
    let cytaa = expansionScale(aa, cdytail)
    let temp16b = expansionScale(cytaa, bdx)
    
    let cytbb = expansionScale(bb, cdytail)
    let temp16c = expansionScale(cytbb, -adx)
    
    let temp32a = expansionSum(temp16a, temp16b)
    let temp48 = expansionSum(temp16c, temp32a)
    fin1 = expansionSum(fin1, temp48)
  }
  
  if ((adxtail != 0.0) || (adytail != 0.0)) {
    let bct: [Double]
    let bctt: [Double]
    if ((bdxtail != 0.0) || (bdytail != 0.0)
        || (cdxtail != 0.0) || (cdytail != 0.0)) {
      let (ti1, ti0) = twoProd(bdxtail, cdy)
      let (tj1, tj0) = twoProd(bdx, cdytail)
      let (u3, u2, u1, u0) = twoTwoSum(ti1, ti0, tj1, tj0)
      let u = [u0, u1, u2, u3]
      let (ti1b, ti0b) = twoProd(cdxtail, -bdy)
      let (tj1b, tj0b) = twoProd(cdx, -bdytail)
      let (v3, v2, v1, v0) = twoTwoSum(ti1b, ti0b, tj1b, tj0b)
      let v = [v0, v1, v2, v3]
      bct = expansionSum(u, v)
      
      let (ti1c, ti0c) = twoProd(bdxtail, cdytail)
      let (tj1c, tj0c) = twoProd(cdxtail, bdytail)
      let (bctt3, bctt2, bctt1, bctt0) = twoTwoDiff(ti1c, ti0c, tj1c, tj0c)
      bctt = [bctt0, bctt1, bctt2, bctt3]
    } else {
      bct = [0.0]
      bctt = [0.0]
    }
    
    if (adxtail != 0.0) {
      var temp16a = expansionScale(axtbc, adxtail)
      let axtbct = expansionScale(bct, adxtail)
      var temp32a = expansionScale(axtbct, 2.0 * adx)
      let temp48 = expansionSum(temp16a, temp32a)
      fin1 = expansionSum(fin1, temp48)
      
      if (bdytail != 0.0) {
        let temp8 = expansionScale(cc, adxtail)
        temp16a = expansionScale(temp8, bdytail)
        fin1 = expansionSum(fin1, temp16a)
      }
      
      if (cdytail != 0.0) {
        let temp8 = expansionScale(bb, -adxtail)
        temp16a = expansionScale(temp8, cdytail)
        fin1 = expansionSum(fin1, temp16a)
      }
      
      temp32a = expansionScale(axtbct, adxtail)
      let axtbctt = expansionScale(bctt, adxtail)
      temp16a = expansionScale(axtbctt, 2.0 * adx)
      let temp16b = expansionScale(axtbctt, adxtail)
      let temp32b = expansionSum(temp16a, temp16b)
      let temp64 = expansionSum(temp32a, temp32b)
      fin1 = expansionSum(fin1, temp64)
    }
    if (adytail != 0.0) {
      var temp16a = expansionScale(aytbc, adytail)
      let aytbct = expansionScale(bct, adytail)
      var temp32a = expansionScale(aytbct, 2.0 * ady)
      let temp48 = expansionSum(temp16a, temp32a)
      fin1 = expansionSum(fin1, temp48)
      
      temp32a = expansionScale(aytbct, adytail)
      let aytbctt = expansionScale(bctt, adytail)
      temp16a = expansionScale(aytbctt, 2.0 * ady)
      let temp16b = expansionScale(aytbctt, adytail)
      let temp32b = expansionSum(temp16a, temp16b)
      let temp64 = expansionSum(temp32a, temp32b)
      fin1 = expansionSum(fin1, temp64)
    }
  }
  if ((bdxtail != 0.0) || (bdytail != 0.0)) {
    let cat: [Double]
    let catt: [Double]
    if ((cdxtail != 0.0) || (cdytail != 0.0)
        || (adxtail != 0.0) || (adytail != 0.0)) {
      let (ti1, ti0) = twoProd(cdxtail, ady)
      let (tj1, tj0) = twoProd(cdx, adytail)
      let (u3, u2, u1, u0) = twoTwoSum(ti1, ti0, tj1, tj0)
      let u = [u0, u1, u2, u3]
      let (ti1b, ti0b) = twoProd(adxtail, -cdy)
      let (tj1b, tj0b) = twoProd(adx, -cdytail)
      let (v3, v2, v1, v0) = twoTwoSum(ti1b, ti0b, tj1b, tj0b)
      let v = [v0, v1, v2, v3]
      cat = expansionSum(u, v)
      
      let (ti1c, ti0c) = twoProd(cdxtail, adytail)
      let (tj1c, tj0c) = twoProd(adxtail, cdytail)
      let (catt3, catt2, catt1, catt0) = twoTwoDiff(ti1c, ti0c, tj1c, tj0c)
      catt = [catt0, catt1, catt2, catt3]
    } else {
      cat = [0.0]
      catt = [0.0]
    }
    
    if (bdxtail != 0.0) {
      var temp16a = expansionScale(bxtca, bdxtail)
      let bxtcat = expansionScale(cat, bdxtail)
      var temp32a = expansionScale(bxtcat, 2.0 * bdx)
      let temp48 = expansionSum(temp16a, temp32a)
      fin1 = expansionSum(fin1, temp48)
      
      if (cdytail != 0.0) {
        let temp8 = expansionScale(aa, bdxtail)
        temp16a = expansionScale(temp8, cdytail)
        fin1 = expansionSum(fin1, temp16a)
      }
      if (adytail != 0.0) {
        let temp8 = expansionScale(cc, -bdxtail)
        temp16a = expansionScale(temp8, adytail)
        fin1 = expansionSum(fin1, temp16a)
      }
      
      temp32a = expansionScale(bxtcat, bdxtail)
      let bxtcatt = expansionScale(catt, bdxtail)
      temp16a = expansionScale(bxtcatt, 2.0 * bdx)
      let temp16b = expansionScale(bxtcatt, bdxtail)
      let temp32b = expansionSum(temp16a, temp16b)
      let temp64 = expansionSum(temp32a, temp32b)
      fin1 = expansionSum(fin1, temp64)
    }
    if (bdytail != 0.0) {
      var temp16a = expansionScale(bytca, bdytail)
      let bytcat = expansionScale(cat, bdytail)
      var temp32a = expansionScale(bytcat, 2.0 * bdy)
      let temp48 = expansionSum(temp16a, temp32a)
      fin1 = expansionSum(fin1, temp48)
      
      temp32a = expansionScale(bytcat, bdytail)
      let bytcatt = expansionScale(catt, bdytail)
      temp16a = expansionScale(bytcatt, 2.0 * bdy)
      let temp16b = expansionScale(bytcatt, bdytail)
      let temp32b = expansionSum(temp16a, temp16b)
      let temp64 = expansionSum(temp32a, temp32b)
      fin1 = expansionSum(fin1, temp64)
    }
  }
  if ((cdxtail != 0.0) || (cdytail != 0.0)) {
    let abt: [Double]
    let abtt: [Double]
    if ((adxtail != 0.0) || (adytail != 0.0)
        || (bdxtail != 0.0) || (bdytail != 0.0)) {
      let (ti1, ti0) = twoProd(adxtail, bdy)
      let (tj1, tj0) = twoProd(adx, bdytail)
      let (u3, u2, u1, u0) = twoTwoSum(ti1, ti0, tj1, tj0)
      let u = [u0, u1, u2, u3]
      let (ti1b, ti0b) = twoProd(bdxtail, -ady)
      let (tj1b, tj0b) = twoProd(bdx, -adytail)
      let (v3, v2, v1, v0) = twoTwoSum(ti1b, ti0b, tj1b, tj0b)
      let v = [v0, v1, v2, v3]
      abt = expansionSum(u, v)
      
      let (ti1c, ti0c) = twoProd(adxtail, bdytail)
      let (tj1c, tj0c) = twoProd(bdxtail, adytail)
      let (abtt3, abtt2, abtt1, abtt0) = twoTwoDiff(ti1c, ti0c, tj1c, tj0c)
      abtt = [abtt0, abtt1, abtt2, abtt3]
    } else {
      abt = [0.0]
      abtt = [0.0]
    }
    
    if (cdxtail != 0.0) {
      var temp16a = expansionScale(cxtab, cdxtail)
      let cxtabt = expansionScale(abt, cdxtail)
      var temp32a = expansionScale(cxtabt, 2.0 * cdx)
      let temp48 = expansionSum(temp16a, temp32a)
      fin1 = expansionSum(fin1, temp48)
      
      if (adytail != 0.0) {
        let temp8 = expansionScale(bb, cdxtail)
        temp16a = expansionScale(temp8, adytail)
        fin1 = expansionSum(fin1, temp16a)
      }
      
      if (bdytail != 0.0) {
        let temp8 = expansionScale(aa, -cdxtail)
        temp16a = expansionScale(temp8, bdytail)
        fin1 = expansionSum(fin1, temp16a)
      }
      
      temp32a = expansionScale(cxtabt, cdxtail)
      let cxtabtt = expansionScale(abtt, cdxtail)
      temp16a = expansionScale(cxtabtt, 2.0 * cdx)
      let temp16b = expansionScale(cxtabtt, cdxtail)
      let temp32b = expansionSum(temp16a, temp16b)
      let temp64 = expansionSum(temp32a, temp32b)
      fin1 = expansionSum(fin1, temp64)
    }
    
    if (cdytail != 0.0) {
      var temp16a = expansionScale(cytab, cdytail)
      let cytabt = expansionScale(abt, cdytail)
      var temp32a = expansionScale(cytabt, 2.0 * cdy)
      let temp48 = expansionSum(temp16a, temp32a)
      fin1 = expansionSum(fin1, temp48)
      
      temp32a = expansionScale(cytabt, cdytail)
      let cytabtt = expansionScale(abtt, cdytail)
      temp16a = expansionScale(cytabtt, 2.0 * cdy)
      let temp16b = expansionScale(cytabtt, cdytail)
      let temp32b = expansionSum(temp16a, temp16b)
      let temp64 = expansionSum(temp32a, temp32b)
      fin1 = expansionSum(fin1, temp64)
    }
  }
  
  return fin1[fin1.count - 1]
}

// MARK: inCircle Predicate

/// Returns where a point lies relative to the circle through three others, exactly.
///
/// The circle passes through `a`, `b`, `c`, which must not be collinear. List them
/// counterclockwise; for a clockwise triangle, ``CirclePosition/inside`` and
/// ``CirclePosition/outside`` swap.
///
/// The result is exact when every coordinate is finite and either 0 or of
/// magnitude in [2⁻²⁰³, 2²⁵⁴].
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
  let adx = a.x - d.x
  let ady = a.y - d.y
  let bdx = b.x - d.x
  let bdy = b.y - d.y
  let cdx = c.x - d.x
  let cdy = c.y - d.y
  
  let bdxcdy = bdx * cdy
  let cdxbdy = cdx * bdy
  let alift = adx * adx + ady * ady
  let cdxady = cdx * ady
  let adxcdy = adx * cdy
  let blift = bdx * bdx + bdy * bdy
  let adxbdy = adx * bdy
  let bdxady = bdx * ady
  let clift = cdx * cdx + cdy * cdy
  
  let det = alift * (bdxcdy - cdxbdy)
  + blift * (cdxady - adxcdy)
  + clift * (adxbdy - bdxady)
  
  // Error bound: Shewchuk 1997, §4.4, Table 5; the permanent is α_a + α_b + α_c.
  let permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift + (abs(cdxady) + abs(adxcdy)) * blift + (abs(adxbdy) + abs(bdxady)) * clift
  let errbound = iccErrBoundA * permanent
  if det > errbound || -det > errbound {
    return CirclePosition(sign: det)
  }
  return CirclePosition(sign: inCircleAdapt(a, b, c, d, permanent: permanent))
}
