//
//  ErrorBoundsTests.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Testing
@testable import RobustPredicates

@Suite("Error bounds")
struct ErrorBoundsTests {
  // predicates.c, exactinit(): ε is the largest power of two with 1 + ε == 1.
  @Test func shewchukEpsilonMatchesItsDefinition() {
    #expect(1.0 + shewchukEpsilon == 1.0)
    #expect(1.0 + 2 * shewchukEpsilon != 1.0)
    #expect(shewchukEpsilon.significandBitPattern == 0)
  }
}
