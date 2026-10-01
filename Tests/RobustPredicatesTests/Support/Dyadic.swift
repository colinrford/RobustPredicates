//
//  Dyadic.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

/// An exact dyadic rational, magnitude · 2^exponent, with an arbitrary-precision
/// magnitude. Every finite Double is one (IEEE 754-2019 §3.3), so sums, differences
/// and products of Doubles are exact. Test oracle only: it shares no code with the
/// library's expansion arithmetic.
struct Dyadic: Equatable {
  private var negative = false
  private var limbs: [UInt32] = []  // little-endian magnitude; empty is zero
  private var exponent = 0

  static let zero = Dyadic()

  init() {}

  init(_ x: Double) {
    precondition(x.isFinite)
    guard x != 0 else { return }
    negative = x < 0
    let m = UInt64(x.significand * 0x1p52)
    limbs = normalized([UInt32(truncatingIfNeeded: m), UInt32(truncatingIfNeeded: m >> 32)])
    exponent = x.exponent - 52
  }

  init(_ n: Int128) {
    negative = n < 0
    var m = n.magnitude
    while m != 0 {
      limbs.append(UInt32(truncatingIfNeeded: m))
      m >>= 32
    }
  }

  var signum: Int { limbs.isEmpty ? 0 : (negative ? -1 : 1) }

  static prefix func - (x: Dyadic) -> Dyadic {
    var r = x
    if !r.limbs.isEmpty { r.negative.toggle() }
    return r
  }

  static func + (a: Dyadic, b: Dyadic) -> Dyadic {
    if a.limbs.isEmpty { return b }
    if b.limbs.isEmpty { return a }
    let e = min(a.exponent, b.exponent)
    let am = shiftedLeft(a.limbs, by: a.exponent - e)
    let bm = shiftedLeft(b.limbs, by: b.exponent - e)
    var r = Dyadic()
    r.exponent = e
    if a.negative == b.negative {
      r.limbs = added(am, bm)
      r.negative = a.negative
    } else {
      switch compared(am, bm) {
      case 0: return .zero
      case 1:
        r.limbs = subtracted(am, bm)
        r.negative = a.negative
      default:
        r.limbs = subtracted(bm, am)
        r.negative = b.negative
      }
    }
    return r
  }

  static func - (a: Dyadic, b: Dyadic) -> Dyadic { a + -b }

  static func * (a: Dyadic, b: Dyadic) -> Dyadic {
    if a.limbs.isEmpty || b.limbs.isEmpty { return .zero }
    var r = Dyadic()
    r.limbs = multiplied(a.limbs, b.limbs)
    r.negative = a.negative != b.negative
    r.exponent = a.exponent + b.exponent
    return r
  }

  static func == (a: Dyadic, b: Dyadic) -> Bool { (a - b).signum == 0 }
}

extension Dyadic {
  /// The exact sum of an expansion's components.
  init(sum components: [Double]) {
    self = components.reduce(.zero) { $0 + Dyadic($1) }
  }
}

// MARK: - Magnitudes (little-endian UInt32 limbs, no high zero limbs)

private func normalized(_ a: [UInt32]) -> [UInt32] {
  var a = a
  while a.last == 0 { a.removeLast() }
  return a
}

private func compared(_ a: [UInt32], _ b: [UInt32]) -> Int {
  if a.count != b.count { return a.count < b.count ? -1 : 1 }
  for i in a.indices.reversed() where a[i] != b[i] {
    return a[i] < b[i] ? -1 : 1
  }
  return 0
}

private func added(_ a: [UInt32], _ b: [UInt32]) -> [UInt32] {
  var r: [UInt32] = []
  r.reserveCapacity(max(a.count, b.count) + 1)
  var carry: UInt64 = 0
  for i in 0..<max(a.count, b.count) {
    let t = UInt64(i < a.count ? a[i] : 0) + UInt64(i < b.count ? b[i] : 0) + carry
    r.append(UInt32(truncatingIfNeeded: t))
    carry = t >> 32
  }
  if carry != 0 { r.append(UInt32(carry)) }
  return r
}

/// a − b, requiring a ≥ b.
private func subtracted(_ a: [UInt32], _ b: [UInt32]) -> [UInt32] {
  var r: [UInt32] = []
  r.reserveCapacity(a.count)
  var borrow: UInt64 = 0
  for i in a.indices {
    let d = UInt64(a[i]) &- UInt64(i < b.count ? b[i] : 0) &- borrow
    r.append(UInt32(truncatingIfNeeded: d))
    borrow = d >> 63
  }
  return normalized(r)
}

private func multiplied(_ a: [UInt32], _ b: [UInt32]) -> [UInt32] {
  var r = [UInt32](repeating: 0, count: a.count + b.count)
  for i in a.indices {
    var carry: UInt64 = 0
    for j in b.indices {
      let t = UInt64(a[i]) * UInt64(b[j]) + UInt64(r[i + j]) + carry
      r[i + j] = UInt32(truncatingIfNeeded: t)
      carry = t >> 32
    }
    r[i + b.count] = UInt32(carry)
  }
  return normalized(r)
}

private func shiftedLeft(_ a: [UInt32], by bits: Int) -> [UInt32] {
  guard bits > 0 else { return a }
  var r = [UInt32](repeating: 0, count: bits / 32)
  let s = bits % 32
  if s == 0 {
    r.append(contentsOf: a)
  } else {
    var carry: UInt32 = 0
    for x in a {
      r.append((x << s) | carry)
      carry = x >> (32 - s)
    }
    if carry != 0 { r.append(carry) }
  }
  return r
}
