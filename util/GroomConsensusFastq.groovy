import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

def cli = new CliBuilder(usage:
        "GroomConsensusFastq [options] input.fastq[.gz] output.fastq[.gz]")
cli.r("Reverse complement (should be used for R2).")

def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

// Inline here

def getWriter = { String outfile ->
    boolean compressed = outfile.endsWith(".gz")
    new BufferedWriter(new OutputStreamWriter(compressed ?
            new GZIPOutputStream(new FileOutputStream(outfile)) : new FileOutputStream(outfile)))
}

def getReader = { String fname ->
    new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fname)) :
            new FileInputStream(fname)))
}

def revCompl = { String seq ->
    def chars = seq.reverse().toCharArray()
    for (int i = 0; i < chars.length; i++) {
        switch (chars[i]) {
            case ((char) 'A'):
                chars[i] = (char) 'T'
                break
            case ((char) 'T'):
                chars[i] = (char) 'A'
                break
            case ((char) 'G'):
                chars[i] = (char) 'C'
                break
            case ((char) 'C'):
                chars[i] = (char) 'G'
                break
            default:
                chars[i] = (char) 'N'
                break
        }
    }
    return new String(chars)
}

def inputFileName = opt.arguments()[0],
    outputFileName = opt.arguments()[1]
boolean rc = opt.r

def reader = getReader(inputFileName), writer = getWriter(outputFileName)

def header

while ((header = reader.readLine()) != null) {
    writer.writeLine(header.replaceAll(" UMI:", "_UMI:"))
    writer.writeLine(rc ? revCompl(reader.readLine()) : reader.readLine())
    writer.writeLine(reader.readLine())
    writer.writeLine(rc ? reader.readLine().reverse() : reader.readLine())
}

writer.close()