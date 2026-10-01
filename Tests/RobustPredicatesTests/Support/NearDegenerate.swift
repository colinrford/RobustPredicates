//
//  NearDegenerate.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

/// `count` consecutive Doubles starting at `x`.
func ulpSteps(from x: Double, count: Int) -> [Double] {
  var steps = [x]
  steps.reserveCapacity(count)
  while steps.count < count { steps.append(steps[steps.count - 1].nextUp) }
  return steps
}

/// The `size × size` grid of consecutive Doubles starting at `p` on each axis.
func ulpGrid(from p: SIMD2<Double>, size: Int) -> [SIMD2<Double>] {
  let xs = ulpSteps(from: p.x, count: size)
  let ys = ulpSteps(from: p.y, count: size)
  return ys.flatMap { y in xs.map { x in SIMD2(x, y) } }
}

/// A rounded point on segment a–b, so (a, b, c) is collinear up to rounding.
func nearCollinearTriple(
  _ rng: inout SplitMix64
) -> (a: SIMD2<Double>, b: SIMD2<Double>, c: SIMD2<Double>) {
  let a = randomPoint(&rng, in: -1...1)
  let b = randomPoint(&rng, in: -1...1)
  let t = Double.random(in: 0...1, using: &rng)
  return (a, b, a + (b - a) * t)
}

/// Four rounded points on one circle.
func nearCocircularQuad(
  _ rng: inout SplitMix64
) -> (a: SIMD2<Double>, b: SIMD2<Double>, c: SIMD2<Double>, d: SIMD2<Double>) {
  let center = randomPoint(&rng, in: -1...1)
  let radius = Double.random(in: 0.5...2, using: &rng)
  func onCircle() -> SIMD2<Double> {
    let t = Double.random(in: -4...4, using: &rng)
    return center + radius * SIMD2((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))
  }
  return (onCircle(), onCircle(), onCircle(), onCircle())
}

/// A rounded point d in the plane through a, b, c.
func nearCoplanarQuad(
  _ rng: inout SplitMix64
) -> (a: SIMD3<Double>, b: SIMD3<Double>, c: SIMD3<Double>, d: SIMD3<Double>) {
  let a = randomPoint3(&rng, in: -1...1)
  let b = randomPoint3(&rng, in: -1...1)
  let c = randomPoint3(&rng, in: -1...1)
  let s = Double.random(in: 0...1, using: &rng)
  let t = Double.random(in: 0...1, using: &rng)
  return (a, b, c, a + (b - a) * s + (c - a) * t)
}
