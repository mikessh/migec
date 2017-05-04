MIGEC="java -Xmx2G -jar ../target/migec-*.jar"

$MIGEC Checkout -ou barcodes.txt trb_R1.fastq.gz trb_R2.fastq.gz checkout/
head checkout/S2-2-beta_R1.fastq
$MIGEC Histogram checkout/ histogram/
$MIGEC AssembleBatch --force-overseq 5 --force-collision-filter --default-mask 0:1 checkout/ histogram/ assemble/
head assemble/S2-1-beta_R2.t5.cf.fastq
$MIGEC CdrBlastBatch -R TRB checkout/ assemble/ cdrblast/
$MIGEC FilterCdrBlastResultsBatch cdrblast/ cdrfinal/
#$MIGEC Report .

if [[ -s _migec_error.log ]]
   then cat _migec_error.log; exit 1
fi
