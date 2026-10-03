// swift-tools-version: 6.2

import PackageDescription

let package = Package(
  name: "PredicateMap",
  platforms: [.macOS(.v15)],
  dependencies: [
    .package(path: "../..")
  ],
  targets: [
    .executableTarget(
      name: "PredicateMap",
      dependencies: [.product(name: "RobustPredicates", package: "RobustPredicates")]
    )
  ]
)
