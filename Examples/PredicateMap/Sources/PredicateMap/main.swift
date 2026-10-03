//
//  main.swift
//  PredicateMap
//
//  Copyright (c) 2026 Colin Ford.
//  SPDX-License-Identifier: MIT
//

// Colours each point of a 256 × 256 grid of consecutive Doubles from (0.5, 0.5)
// by the sign a predicate reports, once evaluated naively and once with
// RobustPredicates. The orient2d setup is from L. Kettner, K. Mehlhorn, S. Pion,
// S. Schirra and C. Yap, "Classroom examples of robustness problems in geometric
// computations", Computational Geometry 40(1):61–78, 2008.
//
// Run from this folder with `swift run`; images are written to docs/images.
// The setup diagrams are TikZ, in Diagrams/.

import CoreGraphics
import Foundation
import ImageIO
import RobustPredicates
import UniformTypeIdentifiers

enum Sign {
  case positive, zero, negative

  init(_ x: Double) { self = x > 0 ? .positive : x < 0 ? .negative : .zero }

  var rgb: [UInt8] {
    switch self {
    case .positive: [0xD9, 0x3A, 0x3A]
    case .zero: [0x2E, 0xA0, 0x4F]
    case .negative: [0x3A, 0x6E, 0xD9]
    }
  }
}

let size = 256
let xs = sequence(first: 0.5, next: \.nextUp).prefix(size)

func map(_ sign: (SIMD2<Double>) -> Sign) -> [Sign] {
  xs.reversed().flatMap { y in xs.map { x in sign(SIMD2(x, y)) } }
}

// MARK: - Predicates

let q = SIMD2(12.0, 12.0), r = SIMD2(24.0, 24.0)
let a = SIMD2(0.5, -22.5), b = SIMD2(23.5, -22.5), c = SIMD2(23.5, 0.5)

// Expanded around p, as in Kettner et al.
func naiveOrient2d(_ p: SIMD2<Double>) -> Sign {
  Sign((q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x))
}

func robustOrient2d(_ p: SIMD2<Double>) -> Sign {
  switch orient2d(p, q, r) {
  case .ccw: .positive
  case .collinear: .zero
  case .cw: .negative
  }
}

func naiveInCircle(_ d: SIMD2<Double>) -> Sign {
  let adx = a.x - d.x, ady = a.y - d.y
  let bdx = b.x - d.x, bdy = b.y - d.y
  let cdx = c.x - d.x, cdy = c.y - d.y
  return Sign((adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
    + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy)
    + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady))
}

func robustInCircle(_ d: SIMD2<Double>) -> Sign {
  switch inCircle(a, b, c, d) {
  case .inside: .positive
  case .on: .zero
  case .outside: .negative
  }
}

// MARK: - Output

let images = URL(filePath: #filePath)
  .deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent()
  .deletingLastPathComponent().deletingLastPathComponent()
  .appending(path: "docs/images")
try FileManager.default.createDirectory(at: images, withIntermediateDirectories: true)

/// Writes `signs` as a PNG, each grid point drawn as a `scale × scale` square.
func writePNG(_ signs: [Sign], to name: String, scale: Int = 2) throws {
  let width = size * scale
  var bytes: [UInt8] = []
  bytes.reserveCapacity(width * width * 3)
  for row in 0..<size {
    let line = signs[row * size..<(row + 1) * size].flatMap { s in
      repeatElement(s.rgb, count: scale).joined()
    }
    for _ in 0..<scale { bytes += line }
  }
  guard
    let provider = CGDataProvider(data: Data(bytes) as CFData),
    let image = CGImage(
      width: width, height: width, bitsPerComponent: 8, bitsPerPixel: 24, bytesPerRow: width * 3,
      space: CGColorSpace(name: CGColorSpace.sRGB)!, bitmapInfo: [], provider: provider,
      decode: nil, shouldInterpolate: false, intent: .defaultIntent),
    let destination = CGImageDestinationCreateWithURL(
      images.appending(path: name) as CFURL, UTType.png.identifier as CFString, 1, nil)
  else { throw CocoaError(.fileWriteUnknown) }
  CGImageDestinationAddImage(destination, image, nil)
  guard CGImageDestinationFinalize(destination) else { throw CocoaError(.fileWriteUnknown) }
}

for (name, sign) in [
  ("orient2d-naive", naiveOrient2d), ("orient2d-robust", robustOrient2d),
  ("incircle-naive", naiveInCircle), ("incircle-robust", robustInCircle),
] as [(String, (SIMD2<Double>) -> Sign)] {
  try writePNG(map(sign), to: "\(name).png")
}
print("Wrote", images.path(percentEncoded: false))
