MIGEC="java -Xmx2G -jar ../target/migec-*.jar"

$MIGEC Checkout -ou barcodes.txt trb_R1.fastq.gz trb_R2.fastq.gz checkout/
$MIGEC Histogram checkout/ histogram/
$MIGEC AssembleBatch --force-overseq 5 --force-collision-filter --default-mask 0:1 checkout/ histogram/ assemble/
$MIGEC CdrBlastBatch -R TRB checkout/ assemble/ cdrblast/
$MIGEC FilterCdrBlastResultsBatch cdrblast/ cdrfinal/
$MIGEC Report .

if [[ -s _migec_error.log ]]
   then cat _migec_error.log; exit 1
fi
