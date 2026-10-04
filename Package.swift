// swift-tools-version: 6.2
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "RobustPredicates",
  platforms: [.iOS(.v18), .macOS(.v15)],
  products: [
    // Products define the executables and libraries a package produces, making them visible to other packages.
    .library(
      name: "RobustPredicates",
      targets: ["RobustPredicates"]
    ),
  ],
  targets: [
    // Targets are the basic building blocks of a package, defining a module or a test suite.
    // Targets can depend on other targets in this package and products from dependencies.
    .target(
      name: "RobustPredicates",
      swiftSettings: [.strictMemorySafety()]
    ),
    .target(
      name: "ShewchukPredicates",
      path: "Tests/ShewchukPredicates",
      exclude: ["predicates.c"]
    ),
    .testTarget(
      name: "RobustPredicatesTests",
      dependencies: ["RobustPredicates", "ShewchukPredicates"],
      resources: [.copy("StagedPredicates")],
      swiftSettings: [.strictMemorySafety()]
    ),
  ]
)
