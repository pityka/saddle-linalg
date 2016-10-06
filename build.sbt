scalaVersion := "2.11.8"

name := "saddle-linalg"

libraryDependencies ++= Seq(
  "org.scala-saddle" %% "saddle-core" % "1.3.4",
  "com.github.fommil.netlib" % "all" % "1.1.2" pomOnly (),
  "org.scalatest" %% "scalatest" % "3.0.0" % "test"
)

reformatOnCompileSettings
