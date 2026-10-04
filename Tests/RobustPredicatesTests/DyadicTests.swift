//
//  DyadicTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing

/// The oracle is checked against Int128 on inputs small enough for it to hold.
@Suite("Dyadic oracle")
struct DyadicTests {

  @Test(arguments: [1 as UInt64, 2])
  func integerDoublesConvertExactly(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<2000 {
      let n = Int64.random(in: -(1 << 53)...(1 << 53), using: &rng)
      #expect(Dyadic(Double(n)) == Dyadic(Int128(n)))
    }
  }

  @Test(arguments: [3 as UInt64, 4])
  func arithmeticMatchesInt128(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<2000 {
      let (a, b, c, d) = (
        Int64.random(in: -(1 << 40)...(1 << 40), using: &rng),
        Int64.random(in: -(1 << 40)...(1 << 40), using: &rng),
        Int64.random(in: -(1 << 40)...(1 << 40), using: &rng),
        Int64.random(in: -(1 << 40)...(1 << 40), using: &rng))
      let want = Int128(a) * Int128(b) - Int128(c) * Int128(d) + Int128(a)
      let got = Dyadic(Double(a)) * Dyadic(Double(b)) - Dyadic(Double(c)) * Dyadic(Double(d))
        + Dyadic(Double(a))
      #expect(got == Dyadic(want))
      #expect(got.signum == want.signum())
    }
  }

  @Test(arguments: [5 as UInt64])
  func powerOfTwoScalingIsExact(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<2000 {
      let x = randomDouble(&rng, exponents: -200...200)
      let k = Int.random(in: -300...300, using: &rng)
      #expect(Dyadic(x) * Dyadic(powerOfTwo(k)) == Dyadic(x * powerOfTwo(k)))
    }
  }

  @Test(arguments: [6 as UInt64])
  func adjacentDoublesAreDistinct(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<2000 {
      let x = randomDouble(&rng, exponents: -300...300)
      #expect(Dyadic(x) != Dyadic(x.nextUp))
      #expect((Dyadic(x.nextUp) - Dyadic(x)).signum == 1)
    }
  }

  @Test func distinguishesRoundedSums() {
    #expect(Dyadic(0.1) + Dyadic(0.2) != Dyadic(0.1 + 0.2))
    #expect(Dyadic(0.5) + Dyadic(0.25) == Dyadic(0.75))
    #expect(Dyadic(Double.leastNonzeroMagnitude) - Dyadic(Double.leastNonzeroMagnitude) == .zero)
  }
}
