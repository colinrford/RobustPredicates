//
//  main.swift
//  PredicateBenchmarks
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

// Times RobustPredicates against Shewchuk's predicates.c on the same inputs.
// Run from this folder with `swift run -c release`.

import Foundation
import RobustPredicates
import Shewchuk

exactinitC()

func sign(_ x: Double) -> Int { x > 0 ? 1 : x < 0 ? -1 : 0 }

extension Orientation {
  var sign: Int {
    switch self {
    case .ccw: 1
    case .collinear: 0
    case .cw: -1
    }
  }
}

extension CirclePosition {
  var sign: Int {
    switch self {
    case .inside: 1
    case .on: 0
    case .outside: -1
    }
  }
}

// predicates.c's orient3d is positive when d is below the plane.
extension PlaneSide {
  var sign: Int {
    switch self {
    case .below: 1
    case .on: 0
    case .above: -1
    }
  }
}

/// The sign of call `i` on the coordinates at `p`.
typealias Call = (_ p: UnsafeMutablePointer<Double>, _ i: Int) -> Int

struct Predicate {
  let name: String
  let inputs: [InputSet]
  let shewchuk: Call
  let slow: Call
  let ours: Call
}

let predicates = [
  Predicate(
    name: "orient2d", inputs: orient2dInputs(),
    shewchuk: { p, i in
      let q = p + 6 * i
      return sign(orient2dC(q, q + 2, q + 4))
    },
    slow: { p, i in
      let q = p + 6 * i
      return sign(orient2dSlowC(q, q + 2, q + 4))
    },
    ours: { p, i in
      let q = p + 6 * i
      return orient2d(SIMD2(q[0], q[1]), SIMD2(q[2], q[3]), SIMD2(q[4], q[5])).sign
    }),
  Predicate(
    name: "inCircle", inputs: inCircleInputs(),
    shewchuk: { p, i in
      let q = p + 8 * i
      return sign(inCircleC(q, q + 2, q + 4, q + 6))
    },
    slow: { p, i in
      let q = p + 8 * i
      return sign(inCircleSlowC(q, q + 2, q + 4, q + 6))
    },
    ours: { p, i in
      let q = p + 8 * i
      return inCircle(
        SIMD2(q[0], q[1]), SIMD2(q[2], q[3]), SIMD2(q[4], q[5]), SIMD2(q[6], q[7])
      ).sign
    }),
  Predicate(
    name: "orient3d", inputs: orient3dInputs(),
    shewchuk: { p, i in
      let q = p + 12 * i
      return sign(orient3dC(q, q + 3, q + 6, q + 9))
    },
    slow: { p, i in
      let q = p + 12 * i
      return sign(orient3dSlowC(q, q + 3, q + 6, q + 9))
    },
    ours: { p, i in
      let q = p + 12 * i
      return orient3d(
        SIMD3(q[0], q[1], q[2]), SIMD3(q[3], q[4], q[5]),
        SIMD3(q[6], q[7], q[8]), SIMD3(q[9], q[10], q[11])
      ).sign
    }),
]

// `--profile <predicate> <inputs> <seconds>` runs only RobustPredicates on one
// input set, e.g. `--profile orient2d "Kettner grid" 10`, for a profiler.
let arguments = CommandLine.arguments
if arguments.count == 5, arguments[1] == "--profile",
  let predicate = predicates.first(where: { $0.name == arguments[2] }),
  var set = predicate.inputs.first(where: { $0.name == arguments[3] }),
  let duration = Double(arguments[4])
{
  let count = set.count
  var checksum = 0
  let clock = ContinuousClock()
  let end = clock.now + .seconds(duration)
  set.coordinates.withUnsafeMutableBufferPointer { buffer in
    while clock.now < end {
      for i in 0..<count { checksum &+= predicate.ours(buffer.baseAddress!, i) }
    }
  }
  print("Checksum: \(checksum).")
  exit(0)
}

let samples = 11
let minimumSample = 0.05  // seconds

func seconds(_ d: Duration) -> Double {
  Double(d.components.seconds) + Double(d.components.attoseconds) * 1e-18
}

/// Nanoseconds per call for each sample, and a checksum so the calls can't be
/// optimized away.
func time(_ call: Call, _ p: UnsafeMutablePointer<Double>, count: Int, repeats: Int) -> (Double, Int) {
  var checksum = 0
  let elapsed = ContinuousClock().measure {
    for _ in 0..<repeats {
      for i in 0..<count { checksum &+= call(p, i) }
    }
  }
  return (seconds(elapsed) * 1e9 / Double(count * repeats), checksum)
}

func median(_ xs: [Double]) -> Double { xs.sorted()[xs.count / 2] }

func format(_ xs: [Double]) -> String {
  String(format: "%.1f (%.1f–%.1f)", median(xs), xs.min()!, xs.max()!)
}

#if arch(arm64)
let arch = "arm64"
#elseif arch(x86_64)
let arch = "x86_64"
#else
let arch = "unknown"
#endif

print("\(ProcessInfo.processInfo.operatingSystemVersionString), \(arch)")
print("""
  Median ns per call over \(samples) samples (min–max). predicates.c is the adaptive \
  predicate; slow is its exact routine without filtering or adaptivity. Ratio is \
  RobustPredicates ÷ predicates.c.

  | Predicate | Inputs | Calls | predicates.c | predicates.c slow | RobustPredicates | Ratio |
  | :--- | :--- | ---: | ---: | ---: | ---: | ---: |
  """)

var disagreements = 0
var checksum = 0
for predicate in predicates {
  for var set in predicate.inputs {
    let count = set.count
    set.coordinates.withUnsafeMutableBufferPointer { buffer in
      let p = buffer.baseAddress!
      for i in 0..<count {
        let want = predicate.shewchuk(p, i)
        if predicate.slow(p, i) != want || predicate.ours(p, i) != want { disagreements += 1 }
      }

      let (once, _) = time(predicate.shewchuk, p, count: count, repeats: 1)
      let repeats = max(1, Int(minimumSample / (once * 1e-9 * Double(count))))
      var shewchuk: [Double] = [], slow: [Double] = [], ours: [Double] = []
      for _ in 0..<samples {
        let (s, a) = time(predicate.shewchuk, p, count: count, repeats: repeats)
        let (w, b) = time(predicate.slow, p, count: count, repeats: repeats)
        let (o, c) = time(predicate.ours, p, count: count, repeats: repeats)
        shewchuk.append(s)
        slow.append(w)
        ours.append(o)
        checksum &+= a &+ b &+ c
      }
      let ratio = String(format: "%.2f", median(ours) / median(shewchuk))
      print(
        "| \(predicate.name) | \(set.name) | \(count) | \(format(shewchuk)) | \(format(slow))",
        "| \(format(ours)) | \(ratio) |")
    }
  }
}
print("\nDisagreements: \(disagreements). Checksum: \(checksum).")
if disagreements > 0 { exit(1) }
