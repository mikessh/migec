/**
 Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.migec

import java.util.jar.JarFile

def version = (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
        /Implementation-Version: (.+)/)[0][1]

println " MiGEC Pipeline V$version "
println ""
println "Run as \$java -cp migec-${version}.jar SCRIPT_NAME arguments"
println ""
println "where SCRIPT_NAME is one of the following:"
println ""
println "[Main pipeline]"
println "com.antigenomics.migec.Checkout"
println "com.antigenomics.migec.Histogram"
println "com.antigenomics.migec.Assemble"
println "com.antigenomics.migec.CdrBlast"
println "com.antigenomics.migec.FilterCdrBlastResults"
println "com.antigenomics.migec.CreateCdrHypermGraph"
println ""
println "[Miscellaneous]"
println "com.antigenomics.migec.GroupByCdr"
println "com.antigenomics.migec.SplitSpikeClonotypes"
println "com.antigenomics.migec.BacktrackSequence"
println "com.antigenomics.migec.HistQ"
println "com.antigenomics.migec.UmiFrequencyTable"