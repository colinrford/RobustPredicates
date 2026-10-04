//
//  StagedPredicatesTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Foundation
import Testing
@testable import RobustPredicates

/// The coordinates and expected determinant sign on each line of a file in
/// StagedPredicates/.
private func cases(_ name: String) throws -> [(x: [Double], sign: Double)] {
  let url = try #require(
    Bundle.module.url(forResource: name, withExtension: "txt", subdirectory: "StagedPredicates"))
  return try String(contentsOf: url, encoding: .utf8).split(separator: "\n").map { line in
    let fields = try line.split(separator: " ").dropFirst().map { try #require(Double($0)) }
    return (Array(fields.dropLast()), fields[fields.count - 1])
  }
}

// Nanevski, Blelloch & Harper's random inputs. Each set includes a few the
// naive determinant gets wrong.
@Suite("Staged predicates data")
struct StagedPredicatesTests {

  @Test func orient2dMatches() throws {
    let cases = try cases("orient2d")
    #expect(cases.count == 1000)
    var wrong: [Int] = []
    var naiveWrong = 0
    for (i, (x, sign)) in cases.enumerated() {
      let a = SIMD2(x[0], x[1]), b = SIMD2(x[2], x[3]), c = SIMD2(x[4], x[5])
      let want = Orientation(sign: sign)
      if orient2d(a, b, c) != want { wrong.append(i + 1) }
      if naiveOrient2d(a, b, c) != want { naiveWrong += 1 }
    }
    #expect(wrong.isEmpty, "first failing lines: \(wrong.prefix(3))")
    #expect(naiveWrong > 0)
  }

  @Test func inCircleMatches() throws {
    let cases = try cases("incircle")
    #expect(cases.count == 1000)
    var wrong: [Int] = []
    var naiveWrong = 0
    for (i, (x, sign)) in cases.enumerated() {
      let a = SIMD2(x[0], x[1]), b = SIMD2(x[2], x[3])
      let c = SIMD2(x[4], x[5]), d = SIMD2(x[6], x[7])
      let want = CirclePosition(sign: sign)
      if inCircle(a, b, c, d) != want { wrong.append(i + 1) }
      if naiveInCircle(a, b, c, d) != want { naiveWrong += 1 }
    }
    #expect(wrong.isEmpty, "first failing lines: \(wrong.prefix(3))")
    #expect(naiveWrong > 0)
  }

  @Test func orient3dMatches() throws {
    let cases = try cases("orient3d")
    #expect(cases.count == 1000)
    var wrong: [Int] = []
    var naiveWrong = 0
    for (i, (x, sign)) in cases.enumerated() {
      let a = SIMD3(x[0], x[1], x[2]), b = SIMD3(x[3], x[4], x[5])
      let c = SIMD3(x[6], x[7], x[8]), d = SIMD3(x[9], x[10], x[11])
      let want = PlaneSide(sign: sign)
      if orient3d(a, b, c, d) != want { wrong.append(i + 1) }
      if naiveOrient3d(a, b, c, d) != want { naiveWrong += 1 }
    }
    #expect(wrong.isEmpty, "first failing lines: \(wrong.prefix(3))")
    #expect(naiveWrong > 0)
  }
}
