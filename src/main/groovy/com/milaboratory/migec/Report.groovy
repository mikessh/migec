/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.milaboratory.migec

def cli = new CliBuilder(usage: "Report [options] output_path'")
cli.h("usage")
cli.c(longOpt: "checkout-path", args: 1, argName: "path", "Path to Checkout directory.")
cli.i(longOpt: "histogram-path", args: 1, argName: "path", "Path to Histogram directory.")
cli.a(longOpt: "assemble-path", args: 1, argName: "path", "Path to Assemble directory.")
cli.b(longOpt: "cdrblast-path", args: 1, argName: "path", "Path to CdrBlast directory.")
cli.f(longOpt: "cdrfinal-path", args: 1, argName: "path", "Path to FilterCdrBlastResults directory.")

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 1) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

def checkoutPath = opt.c ?: "checkout/", histogramPath = opt.i ?: "histogram/",
    assemblePath = opt.a ?: "assemble/", cdrBlastPath = opt.b ?: "cdrblast/",
    cdrFinalPath = opt.f ?: "cdrfinal/", outputPath = opt.arguments()[0]

def rmdLines = new InputStreamReader(this.class.classLoader.getResourceAsStream("migec_summary.Rmd")).readLines().join("\n")

def replacePath = { String module, String path ->
    new File(path).exists() ?
            rmdLines.replace("__${module}__", path) :
            rmdLines.replace("\"__${module}__\"", "NULL")
}

rmdLines = replacePath("checkout", checkoutPath)
rmdLines = replacePath("histogram", histogramPath)
rmdLines = replacePath("assemble", assemblePath)
rmdLines = replacePath("cdrblast", cdrBlastPath)
rmdLines = replacePath("cdrfinal", cdrFinalPath)

if (new File(outputPath).isDirectory())
    outputPath += "/migec_summary.Rmd"
else if (outputPath.toLowerCase().endsWith(".html"))
    outputPath = outputPath[0..-6]

new File(outputPath).withPrintWriter { pw ->
    pw.println(rmdLines)
}

def scriptOutputPath = outputPath + ".generate.R"

new File(scriptOutputPath).withPrintWriter { pw ->
    pw.println(
            "require(\"rmarkdown\")\n" +
                    "render(\"$outputPath\")"
    )
}

def proc = "Rscript $scriptOutputPath".execute()

proc.in.eachLine {
    println(it)
}

proc.out.close()
proc.waitFor()

if (proc.exitValue()) {
    println "[ERROR] See log below:\n${proc.getErrorStream()}"
    println "[NOTE] Intermediate files were not deleted"
} else {
    new File(scriptOutputPath).delete()
    new File(outputPath).delete()
}