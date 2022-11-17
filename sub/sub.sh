## Successful
input="gs://cpg-tob-wgs-test/mt/v7.mt"
outputDir="mt-call-rate-test"
chr=22
analysis-runner --dataset sgs-somatic-mtn --access-level test  --output-dir "$outputDir" --description "pVCF" sub_dataproc.py \
    --script "scripts/mt_vcf.py --dataset '$input' --chrom chr${chr} --output Check_chr${chr}.vcf.bgz" \
    --jobname convert$chr


#### run export pVCF
input="gs://cpg-tob-wgs-test/mt/v7.mt"
outputDir="mt-call-rate-test"
chr=22
analysis-runner --dataset sgs-somatic-mtn --access-level test  --output-dir "$outputDir" --description "test" sub_dataproc.py \
    --script "scripts/mt_vcf2.py --dataset '$input' --chrom chr${chr} --output Check_chr${chr}.mt" \
    --jobname convert$chr


##### run GC
chrom="22"
vcf="/data/pVCF3/chr$chrom.vcf.bgz"
ovcf="/data/pVCF3/chr${chrom}_gc.vcf.gz"
ref="/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
docker="us.gcr.io/mccarroll-mocha"
image="bcftools:1.14-20220112"

analysis-runner --dataset sgs-somatic-mtn --access-level test  --output-dir "pVCF3" --description "output" sub.py \
--cmd "bcftools +mochatools $vcf -Oz -o $ovcf -- -t GC -f $ref" \
--image "$docker/$image" \
--mount "cpg-sgs-somatic-mtn-test => /data; cpg-sgs-somatic-mtn-test-upload => /ref" \
--readonly false \
--jobname GC$chrom
