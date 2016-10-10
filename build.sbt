scalaVersion := "2.11.8"

name := "saddle-linalg"

organization := "io.github.pityka"

version := "0.0.2"

libraryDependencies ++= Seq(
  "org.scala-saddle" %% "saddle-core" % "1.3.4",
  "com.github.fommil.netlib" % "all" % "1.1.2" pomOnly (),
  "org.scalatest" %% "scalatest" % "3.0.0" % "test"
)

reformatOnCompileSettings

pomExtra in Global := {
  <url>https://pityka.github.io/saddle-linalg</url>
  <licenses>
    <license>
      <name>MIT</name>
      <url>https://opensource.org/licenses/MIT</url>
    </license>
  </licenses>
  <scm>
    <connection>scm:git:github.com/pityka/saddle-linalg</connection>
    <developerConnection>scm:git:git@github.com:pityka/saddle-linalg</developerConnection>
    <url>github.com/pityka/saddle-linalg</url>
  </scm>
  <developers>
    <developer>
      <id>pityka</id>
      <name>Istvan Bartha</name>
      <url>https://pityka.github.io/saddle-linalg/</url>
    </developer>
  </developers>
}
