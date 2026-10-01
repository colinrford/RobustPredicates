//
//  OverlapProbe.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

import Synchronization

/// Counts the tasks inside `measure` and records the most seen at once, to show
/// that work split across a task group actually ran in parallel.
///
/// Until a second task arrives, a task waits up to `patience` before running
/// its body. Chunks short enough to finish before the next one is scheduled
/// would otherwise run one at a time, and the overlap would depend on luck.
final class OverlapProbe: Sendable {
  private let active = Atomic<Int>(0)
  private let highWater = Atomic<Int>(0)
  private let patience: Duration

  init(patience: Duration = .seconds(2)) {
    self.patience = patience
  }

  func measure<T>(_ body: () throws -> T) rethrows -> T {
    let now = active.add(1, ordering: .relaxed).newValue
    highWater.max(now, ordering: .relaxed)
    defer { active.subtract(1, ordering: .relaxed) }

    let deadline = ContinuousClock.now + patience
    while highWater.load(ordering: .relaxed) < 2 && ContinuousClock.now < deadline {}
    return try body()
  }

  /// The most tasks that were inside `measure` at the same time.
  var peak: Int { highWater.load(ordering: .relaxed) }
}
