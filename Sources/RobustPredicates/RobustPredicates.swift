//
//  RobustPredicates.swift
//  RobustPredicates
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

// Filtered-exact geometric predicates after Shewchuk, "Adaptive Precision
// Floating-Point Arithmetic and Fast Robust Geometric Predicates" (1997).

/// The winding order or orientation of 3 points.
public enum Orientation: Sendable, Equatable {
  /// The points are in standard or positive orientation, i.e. they wind counterclockwise.
  case ccw
  /// Degenerate case where all points lie on the same line.
  case collinear
  /// Opposite of standard, negative orientation i.e. points wind clockwise.
  case cw
}

/// Position of a point `d` relative to the circle through three positively oriented points `a`, `b`, `c`.
///
/// For negatively-oriented `(a, b, c)` the cases ``inside`` and ``outside`` swap (determinant sign).
public enum CirclePosition: Sendable, Equatable {
  /// `d` lies strictly inside the circle.
  case inside
  /// Degenerate case where `d` lies on the the circle, i.e. the four points are cocircular.
  case on
  /// `d` lies strictly outside the circle.
  case outside
}

/// The side of the plane through three points `a`, `b`, `c`, a fourth point `d` lies.
///
/// Plane oriented by normal `(b−a)×(c−a)`.
public enum PlaneSide: Sendable, Equatable {
  /// The point `d` is on the side the normal points toward.
  case above
  /// The four points are coplanar.
  case on
  /// The point `d` lies on the side opposite the normal.
  case below
}

// Shewchuk's ε: the largest power of two with 1 + ε == 1.
@inlinable
var shewchukEpsilon: Double {
  Double(sign: .plus, exponent: -(Double.significandBitCount + 1), significand: 1)
}
