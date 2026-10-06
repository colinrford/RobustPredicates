//
//  Shewchuk.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

@testable import RobustPredicates
import ShewchukPredicates

// exactinit sets predicates.c's globals, so it runs once, before any predicate.
private let initialized: Void = exactinitC()

func shewchukOrient2d(_ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>) -> Orientation {
  _ = initialized
  var pa = [a.x, a.y], pb = [b.x, b.y], pc = [c.x, c.y]
  return Orientation(sign: unsafe orient2dC(&pa, &pb, &pc))
}

func shewchukOrient2dAdapt(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, detsum: Double
) -> Double {
  _ = initialized
  var pa = [a.x, a.y], pb = [b.x, b.y], pc = [c.x, c.y]
  return unsafe orient2dAdaptC(&pa, &pb, &pc, detsum: detsum)
}

func shewchukInCircle(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>
) -> CirclePosition {
  _ = initialized
  var pa = [a.x, a.y], pb = [b.x, b.y], pc = [c.x, c.y], pd = [d.x, d.y]
  return CirclePosition(sign: unsafe inCircleC(&pa, &pb, &pc, &pd))
}

func shewchukInCircleAdapt(
  _ a: SIMD2<Double>, _ b: SIMD2<Double>, _ c: SIMD2<Double>, _ d: SIMD2<Double>, permanent: Double
) -> Double {
  _ = initialized
  var pa = [a.x, a.y], pb = [b.x, b.y], pc = [c.x, c.y], pd = [d.x, d.y]
  return unsafe inCircleAdaptC(&pa, &pb, &pc, &pd, permanent: permanent)
}

func shewchukOrient3d(
  _ a: SIMD3<Double>, _ b: SIMD3<Double>, _ c: SIMD3<Double>, _ d: SIMD3<Double>
) -> PlaneSide {
  _ = initialized
  var pa = [a.x, a.y, a.z], pb = [b.x, b.y, b.z], pc = [c.x, c.y, c.z], pd = [d.x, d.y, d.z]
  return PlaneSide(sign: unsafe orient3dC(&pa, &pb, &pc, &pd))
}
