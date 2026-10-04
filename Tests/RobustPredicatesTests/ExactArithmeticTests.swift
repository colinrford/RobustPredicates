//
//  ExactArithmeticTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

/// Nonoverlapping and increasing in magnitude, with zeros eliminated except a
/// lone zero (Shewchuk 1997, §2.1, §2.4).
private func isWellFormed(_ e: [Double]) -> Bool {
  if e.count <= 1 { return true }
  guard !e.contains(0) else { return false }
  return zip(e, e.dropFirst()).allSatisfy { small, large in abs(small) < lowestSetBit(large) }
}

/// The value of the least significant nonzero bit of x.
private func lowestSetBit(_ x: Double) -> Double {
  let m = UInt64(abs(x).significand * 0x1p52)
  return powerOfTwo(x.exponent - 52 + m.trailingZeroBitCount)
}

/// An expansion of up to `terms` twoDiff pairs with widely spread exponents.
private func randomExpansion(_ rng: inout SplitMix64, terms: Int) -> [Double] {
  var e: [Double] = []
  for _ in 0..<Int.random(in: 1...terms, using: &rng) {
    let pair = twoDiffE(randomDouble(&rng, exponents: -60...60), randomDouble(&rng, exponents: -60...60))
    e = expansionSum(e, pair)
  }
  return e
}

@Suite("Exact arithmetic")
struct ExactArithmeticTests {
  
  // Knuth, TAOCP Vol. 2, 3rd ed., §4.2.2, Theorem B; Dekker 1971.
  @Test(arguments: [81 as UInt64, 82])
  func twoSumAndTwoDiffAreExact(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var inexact = 0
    for _ in 0..<5000 {
      let a = randomDouble(&rng, exponents: -80...80)
      let b = randomDouble(&rng, exponents: -80...80)
      
      let s = twoSum(a, b)
      if s.lo != 0 { inexact += 1 }
      #expect(s.hi == a + b)
      #expect(Dyadic(s.hi) + Dyadic(s.lo) == Dyadic(a) + Dyadic(b))
      #expect(abs(s.lo) <= s.hi.ulp / 2)
      
      let d = twoDiff(a, b)
      #expect(d.hi == a - b)
      #expect(Dyadic(d.hi) + Dyadic(d.lo) == Dyadic(a) - Dyadic(b))
      #expect(abs(d.lo) <= d.hi.ulp / 2)
    }
    #expect(inexact > 1000, "most random sums must round for this to test anything")
  }
  
  // Ogita, Rump & Oishi 2005: with a correctly rounded FMA, the product's error is exact.
  @Test(arguments: [83 as UInt64, 84])
  func twoProdIsExact(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var inexact = 0
    for _ in 0..<5000 {
      let a = randomDouble(&rng, exponents: -200...200)
      let b = randomDouble(&rng, exponents: -200...200)
      let p = twoProd(a, b)
      if p.lo != 0 { inexact += 1 }
      #expect(p.hi == a * b)
      #expect(Dyadic(p.hi) + Dyadic(p.lo) == Dyadic(a) * Dyadic(b))
    }
    #expect(inexact > 1000, "most random products must round for this to test anything")
  }
  
  // Fails if the optimizer rewrites (a + b) − a as b, which zeroes the low part.
  @Test func twoSumKeepsTheRoundedOffPart() {
    let s = twoSum(0x1p53, 0.5)
    #expect(s.hi == 0x1p53)
    #expect(s.lo == 0.5)
  }
  
  @Test(arguments: [85 as UInt64, 86])
  func twoDiffEIsAWellFormedExpansion(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<5000 {
      let a = randomDouble(&rng, exponents: -60...60)
      let b = randomDouble(&rng, exponents: -60...60)
      let e = twoDiffE(a, b)
      #expect(isWellFormed(e))
      #expect(Dyadic(sum: e) == Dyadic(a) - Dyadic(b))
    }
  }
  
  @Test(arguments: [87 as UInt64, 88, 89])
  func expansionOperationsAreExactAndWellFormed(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    var longest = 0
    for _ in 0..<300 {
      let e = randomExpansion(&rng, terms: 4)
      let f = randomExpansion(&rng, terms: 4)
      let b = randomDouble(&rng, exponents: -60...60)
      #expect(isWellFormed(e) && isWellFormed(f))
      longest = max(longest, e.count)
      
      let sum = expansionSum(e, f)
      #expect(isWellFormed(sum))
      #expect(Dyadic(sum: sum) == Dyadic(sum: e) + Dyadic(sum: f))
      
      let scaled = expansionScale(e, b)
      #expect(isWellFormed(scaled))
      #expect(Dyadic(sum: scaled) == Dyadic(sum: e) * Dyadic(b))
      
      let product = expansionProduct(e, f)
      #expect(isWellFormed(product))
      #expect(Dyadic(sum: product) == Dyadic(sum: e) * Dyadic(sum: f))
      
      #expect(Dyadic(sum: expansionNegate(e)) == -Dyadic(sum: e))
    }
    #expect(longest >= 4, "inputs must include multi-component expansions")
  }
  
  @Test(arguments: [90 as UInt64])
  func mostSignificantComponentHasTheSignOfTheValue(seed: UInt64) {
    var rng = SplitMix64(seed: seed)
    for _ in 0..<1000 {
      let e = randomExpansion(&rng, terms: 4)
      let f = randomExpansion(&rng, terms: 4)
      let value = expansionSum(e, expansionNegate(f))
      let sign = mostSignificantComponent(value)
      #expect((sign > 0 ? 1 : sign < 0 ? -1 : 0) == (Dyadic(sum: e) - Dyadic(sum: f)).signum)
    }
  }
  
  // (e + f) − e − f regroups the same value through different components, so
  // it only reaches exactly zero if every intermediate sum is exact.
  @Test func cancellationAcrossRegroupingGivesZero() {
    var rng = SplitMix64(seed: 91)
    for _ in 0..<200 {
      let e = randomExpansion(&rng, terms: 4)
      let f = randomExpansion(&rng, terms: 4)
      let g = expansionSum(e, f)
      #expect(mostSignificantComponent(expansionSum(expansionSum(g, expansionNegate(e)), expansionNegate(f))) == 0)
    }
  }
}
