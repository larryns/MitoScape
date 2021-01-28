/*
 * Copyright 2021 Larry N. Singh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

name := "MitoScape"

version := "0.1"

scalaVersion := "2.12.12"

val sparkVersion = "3.0.0"
val hadoopVersion = "3.2.1"

/* Spark requirements.
 */
libraryDependencies ++= Seq(
	"org.apache.spark" %% "spark-sql" % sparkVersion excludeAll(
		ExclusionRule(organization = "jakarta.xml.bind"),
		ExclusionRule(organization = "org.glassfish.hk2.external")
	),
	"org.apache.spark" %% "spark-mllib" % sparkVersion excludeAll(
		ExclusionRule(organization = "jakarta.xml.bind"),
		ExclusionRule(organization = "org.glassfish.hk2.external")
	),
	"org.bdgenomics.adam" %% "adam-core-spark3" % "0.32.0" excludeAll(
	),
	"org.apache.hadoop" % "hadoop-common" % hadoopVersion excludeAll(
		ExclusionRule(organization = "com.sun.jersey"),
		ExclusionRule(organization = "org.eclipse.jetty")
	)
)

/* Compiler settings. Use scalac -X for other options and their description.
 * See Here for more info at:
 * http://www.scala-lang.org/files/archive/nightly/docs/manual/html/scalac.html
 */
scalacOptions ++= List("-feature", "-deprecation", "-unchecked", "-Xlint")

/*********************  Dockerizing  ***********************/

mainClass in (Compile, packageBin) := Some("MitoScape.MTClassify")

enablePlugins(DockerPlugin, JavaAppPackaging)

dockerBaseImage := "openjdk:8-jre-slim"
maintainer := "singhln@chop.edu"
exportJars := true

/* NOTE: The COPY command only works if you have picard.jar in the same 
 * directory as the Dockerfile, i.e. MTClassifier/target/docker/stage
 * but this directory is only created if you try to run docker once. This is
 * the only workaround I've found. 
 *
 * Option 1: Add a line:
 * mappings in Universal ++= directory/file(<file>).
 * 
 * Option 2: Just add your files to src/universal and they'll be copied auto-
 * matically. This option is what I chose.
 *
 */
dockerCommands ++= Seq(
)

// For assembly
test in assembly := {} // Don't include tests

// include all scala runtime JARs
assemblyOption in assembly := (assemblyOption in assembly).value.copy(includeScala = true, includeDependency = true)
assemblyJarName in assembly := "MitoScapeClassify.jar"
mainClass in assembly := Some("MitoScape.MTClassify")

assemblyMergeStrategy in assembly := {
	case PathList("module-info.class") => MergeStrategy.discard
	case PathList("org", "apache", "parquet", "avro", xs @ _*) => MergeStrategy.last
	case PathList("org", "apache", "spark", "unused", "UnusedStubClass.class") => MergeStrategy.first
	case PathList("org", "apache", "commons", "logging", xs @_ *) => MergeStrategy.first
	case PathList("com", "esotericsoftware", xs @ _*) => MergeStrategy.first
	case x => 
		val oldStrategy = (assemblyMergeStrategy in assembly).value
		oldStrategy(x)
}

