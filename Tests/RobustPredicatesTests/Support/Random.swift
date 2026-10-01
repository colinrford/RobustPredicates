//
//  Random.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

/// Seeded RNG for reproducible property tests (SplitMix64).
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

/// A Double with a uniformly random sign, 52-bit fraction, and binary exponent.
func randomDouble(_ rng: inout SplitMix64, exponents: ClosedRange<Int>) -> Double {
  let significand = 1 + Double(rng.next() >> 12) * 0x1p-52
  let sign: FloatingPointSign = rng.next() & 1 == 0 ? .plus : .minus
  return Double(sign: sign, exponent: Int.random(in: exponents, using: &rng), significand: significand)
}

func randomPoint(_ rng: inout SplitMix64, in range: ClosedRange<Double>) -> SIMD2<Double> {
  SIMD2(Double.random(in: range, using: &rng), Double.random(in: range, using: &rng))
}

func randomPoint3(_ rng: inout SplitMix64, in range: ClosedRange<Double>) -> SIMD3<Double> {
  SIMD3(Double.random(in: range, using: &rng),
        Double.random(in: range, using: &rng),
        Double.random(in: range, using: &rng))
}

func randomIntegerPoint(_ rng: inout SplitMix64, in range: ClosedRange<Int64>) -> SIMD2<Double> {
  SIMD2(Double(Int64.random(in: range, using: &rng)), Double(Int64.random(in: range, using: &rng)))
}

func randomIntegerPoint3(_ rng: inout SplitMix64, in range: ClosedRange<Int64>) -> SIMD3<Double> {
  SIMD3(Double(Int64.random(in: range, using: &rng)),
        Double(Int64.random(in: range, using: &rng)),
        Double(Int64.random(in: range, using: &rng)))
}
