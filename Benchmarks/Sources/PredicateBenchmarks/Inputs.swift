//
//  Inputs.swift
//  PredicateBenchmarks
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

/// The points of each call, flattened: `stride` coordinates per call, in the
/// order predicates.c takes them.
struct InputSet {
  let name: String
  let stride: Int
  var coordinates: [Double] = []

  var count: Int { coordinates.count / stride }
}

/// SplitMix64, as in the tests.
struct SplitMix64: RandomNumberGenerator {
  var state: UInt64
  init(seed: UInt64) { state = seed }
  mutating func next() -> UInt64 {
    state &+= 0x9E3779B97F4A7C15
    var z = state
    z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
    z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
    return z ^ (z >> 31)
  }
}

/// `count` consecutive Doubles from `x`.
func ulpSteps(from x: Double, count: Int) -> [Double] {
  Array(sequence(first: x, next: \.nextUp).prefix(count))
}

private func random(_ rng: inout SplitMix64, _ n: Int) -> [Double] {
  (0..<n).map { _ in Double.random(in: -1...1, using: &rng) }
}

private let randomCalls = 100_000
private let grids = 400  // near-degenerate configurations, each with a 16 × 16 grid

func orient2dInputs() -> [InputSet] {
  var rng = SplitMix64(seed: 1)
  var easy = InputSet(name: "random", stride: 6)
  for _ in 0..<randomCalls { easy.coordinates += random(&rng, 6) }

  // Kettner, Mehlhorn, Pion, Schirra & Yap 2008.
  var kettner = InputSet(name: "Kettner grid", stride: 6)
  let steps = ulpSteps(from: 0.5, count: 256)
  for y in steps { for x in steps { kettner.coordinates += [x, y, 12, 12, 24, 24] } }

  var near = InputSet(name: "near-collinear", stride: 6)
  for _ in 0..<grids {
    let a = random(&rng, 2), b = random(&rng, 2)
    let t = Double.random(in: 0...1, using: &rng)
    let c = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]
    for y in ulpSteps(from: c[1], count: 16) {
      for x in ulpSteps(from: c[0], count: 16) { near.coordinates += a + b + [x, y] }
    }
  }
  return [easy, kettner, near]
}

func inCircleInputs() -> [InputSet] {
  var rng = SplitMix64(seed: 2)
  var easy = InputSet(name: "random", stride: 8)
  for _ in 0..<randomCalls { easy.coordinates += random(&rng, 8) }

  var near = InputSet(name: "near-cocircular", stride: 8)
  for _ in 0..<grids {
    let center = random(&rng, 2)
    let radius = Double.random(in: 0.5...2, using: &rng)
    func onCircle() -> [Double] {
      let t = Double.random(in: -4...4, using: &rng)
      return [center[0] + radius * (1 - t * t) / (1 + t * t), center[1] + radius * 2 * t / (1 + t * t)]
    }
    let abc = onCircle() + onCircle() + onCircle(), d = onCircle()
    for y in ulpSteps(from: d[1], count: 16) {
      for x in ulpSteps(from: d[0], count: 16) { near.coordinates += abc + [x, y] }
    }
  }
  return [easy, near]
}

func orient3dInputs() -> [InputSet] {
  var rng = SplitMix64(seed: 3)
  var easy = InputSet(name: "random", stride: 12)
  for _ in 0..<randomCalls { easy.coordinates += random(&rng, 12) }

  var near = InputSet(name: "near-coplanar", stride: 12)
  for _ in 0..<grids {
    let a = random(&rng, 3), b = random(&rng, 3), c = random(&rng, 3)
    let s = Double.random(in: 0...1, using: &rng), t = Double.random(in: 0...1, using: &rng)
    let d = (0..<3).map { a[$0] + (b[$0] - a[$0]) * s + (c[$0] - a[$0]) * t }
    for z in ulpSteps(from: d[2], count: 16) {
      for x in ulpSteps(from: d[0], count: 16) { near.coordinates += a + b + c + [x, d[1], z] }
    }
  }
  return [easy, near]
}
