// swift-tools-version: 6.2

import PackageDescription

let package = Package(
  name: "Benchmarks",
  platforms: [.macOS(.v15)],
  dependencies: [
    .package(path: "..")
  ],
  targets: [
    .target(name: "Shewchuk"),
    .executableTarget(
      name: "PredicateBenchmarks",
      dependencies: [
        .product(name: "RobustPredicates", package: "RobustPredicates"),
        "Shewchuk",
      ]
    ),
  ]
)
