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

package com.milaboratory.migec

import java.util.jar.JarFile

def version = (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
        /Implementation-Version: (.+)/)[0][1]

println " MiGEC Pipeline V$version "
println ""
println "Run as \$java -cp com.milaboratory.migec-${version}.jar SCRIPT_NAME arguments"
println ""
println "where SCRIPT_NAME is one of the following:"
println ""
println "[Main pipeline]"
println "com.milaboratory.com.milaboratory.migec.Checkout"
println "com.milaboratory.com.milaboratory.migec.Histogram"
println "com.milaboratory.com.milaboratory.migec.Assemble"
println "com.milaboratory.com.milaboratory.migec.CdrBlast"
println "com.milaboratory.com.milaboratory.migec.FilterCdrBlastResults"
println "com.milaboratory.com.milaboratory.migec.CreateCdrHypermGraph"
println ""
println "[Miscellaneous]"
println "com.milaboratory.com.milaboratory.migec.GroupByCdr"
println "com.milaboratory.com.milaboratory.migec.SplitSpikeClonotypes"
println "com.milaboratory.com.milaboratory.migec.BacktrackSequence"
println "com.milaboratory.com.milaboratory.migec.HistQ"
println "com.milaboratory.com.milaboratory.migec.UmiFrequencyTable"