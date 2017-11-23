name := "saddle-linalg"

organization := "io.github.pityka"

version := "0.0.23"

libraryDependencies ++= Seq(
  "io.github.pityka" %% "saddle-core-fork" % "1.3.4-fork1",
  "com.github.fommil.netlib" % "all" % "1.1.2" pomOnly (),
  "net.sourceforge.f2j" % "arpack_combined_all" % "0.1",
  "org.scalatest" %% "scalatest" % "3.0.0" % "test"
)

scalaVersion := "2.11.11"

crossScalaVersions := Seq("2.12.4", "2.11.11")

scalafmtOnCompile in ThisBuild := true

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
