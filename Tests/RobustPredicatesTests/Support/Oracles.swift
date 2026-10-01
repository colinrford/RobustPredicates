//
//  Oracles.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

@testable import RobustPredicates

// Exact determinant signs over `Dyadic`. Each convention is written out here
// rather than borrowed from the library's `init(sign:)`.

func orient2dOracle(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Orientation {
  let acx = Dyadic(a.x) - Dyadic(c.x), acy = Dyadic(a.y) - Dyadic(c.y)
  let bcx = Dyadic(b.x) - Dyadic(c.x), bcy = Dyadic(b.y) - Dyadic(c.y)
  switch (acx * bcy - acy * bcx).signum {
  case 1: return .ccw
  case -1: return .cw
  default: return .collinear
  }
}

func inCircleOracle(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>
) -> CirclePosition {
  let adx = Dyadic(a.x) - Dyadic(d.x), ady = Dyadic(a.y) - Dyadic(d.y)
  let bdx = Dyadic(b.x) - Dyadic(d.x), bdy = Dyadic(b.y) - Dyadic(d.y)
  let cdx = Dyadic(c.x) - Dyadic(d.x), cdy = Dyadic(c.y) - Dyadic(d.y)
  let det = (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
    + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy)
    + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady)
  switch det.signum {
  case 1: return .inside
  case -1: return .outside
  default: return .on
  }
}

func orient3dOracle(
  _ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>
) -> PlaneSide {
  let adx = Dyadic(a.x) - Dyadic(d.x), ady = Dyadic(a.y) - Dyadic(d.y), adz = Dyadic(a.z) - Dyadic(d.z)
  let bdx = Dyadic(b.x) - Dyadic(d.x), bdy = Dyadic(b.y) - Dyadic(d.y), bdz = Dyadic(b.z) - Dyadic(d.z)
  let cdx = Dyadic(c.x) - Dyadic(d.x), cdy = Dyadic(c.y) - Dyadic(d.y), cdz = Dyadic(c.z) - Dyadic(d.z)
  let det = adz * (bdx * cdy - cdx * bdy)
    + bdz * (cdx * ady - adx * cdy)
    + cdz * (adx * bdy - bdx * ady)
  switch det.signum {
  case 1: return .below
  case -1: return .above
  default: return .on
  }
}

// The fast-path determinants with no error bound: what the filter guards against.

func naiveOrient2d(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Orientation {
  Orientation(sign: (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x))
}

func naiveInCircle(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>
) -> CirclePosition {
  let adx = a.x - d.x, ady = a.y - d.y
  let bdx = b.x - d.x, bdy = b.y - d.y
  let cdx = c.x - d.x, cdy = c.y - d.y
  return CirclePosition(sign: (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
    + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy)
    + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady))
}

func naiveOrient3d(
  _ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>
) -> PlaneSide {
  let adx = a.x - d.x, ady = a.y - d.y, adz = a.z - d.z
  let bdx = b.x - d.x, bdy = b.y - d.y, bdz = b.z - d.z
  let cdx = c.x - d.x, cdy = c.y - d.y, cdz = c.z - d.z
  return PlaneSide(sign: adz * (bdx * cdy - cdx * bdy)
    + bdz * (cdx * ady - adx * cdy)
    + cdz * (adx * bdy - bdx * ady))
}

// MARK: - Symmetry helpers

extension Orientation {
  var flipped: Orientation {
    switch self {
    case .ccw: .cw
    case .cw: .ccw
    case .collinear: .collinear
    }
  }
}

extension CirclePosition {
  var flipped: CirclePosition {
    switch self {
    case .inside: .outside
    case .outside: .inside
    case .on: .on
    }
  }
}

extension PlaneSide {
  var flipped: PlaneSide {
    switch self {
    case .above: .below
    case .below: .above
    case .on: .on
    }
  }
}

/// 2^k as a Double; multiplying by it is exact absent overflow and underflow.
func powerOfTwo(_ k: Int) -> Double {
  Double(sign: .plus, exponent: k, significand: 1)
}

func swappedXY(_ p: SIMD2<Double>) -> SIMD2<Double> { SIMD2(p.y, p.x) }
func swappedXY(_ p: SIMD3<Double>) -> SIMD3<Double> { SIMD3(p.y, p.x, p.z) }
