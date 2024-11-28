use warnings;
use strict;

my $name = "gallus_9d_lung";
my $postfix = "fq.gz";
my @fastq_1s = <"/data02/zhangfenghua/project/cuttag_gallus_11d19d_lung/gallus_9d/rawdata/*/*1.$postfix"> ; # Path to clean fasta
my $refIdx = "/data02/zhangfenghua/project/cuttag_gallus_11d19d_lung/reference/chicken"; # Path to reference index. U need to index reference first (bowtie2-build xxx.fasta.gz name)
my $tssref = "/data02/zhangfenghua/project/cuttag_gallus_11d19d_lung/reference/Gallus.TSS_new.ref.bed";

my $fastp = "/public/home/zhuchenglong/software/fastp/fastp";
my $fastpCMD = "--detect_adapter_for_pe --trim_poly_g --length_required 30 --thread 16 --dedup --dup_calc_accuracy 5";
my $bowtie = "/public/home/zhuchenglong/software/bowtie2/bowtie2-2.3.4.3-linux-x86_64/bowtie2";
my $picard = "/public/home/zhuchenglong/software/java/jdk-11.0.12/bin/java -jar /public/home//xuwenjie/software/picard/picard.jar MarkDuplicates";
my $bamCoverage = "/public/home/zhangfenghua/miniconda3/envs/bio/bin/bamCoverage";
my $bedtools = "/public/home/zhangfenghua/miniconda3/envs/bio/bin/bedtools";
my $macs3 = "/public/home/zhangfenghua/miniconda3/envs/macs3/bin/macs3";
my $gatk = "/public/home/wangkun/software/gatk/gatk-4.1.9.0/gatk CollectInsertSizeMetrics";

my $outdir = "$name.output";
my $qualiControl = "$name.qualiControl";
my $callpeak = "$name.callpeak";
`mkdir -p $outdir` if (!-d $outdir);
`mkdir -p $qualiControl` if (!-d $qualiControl);
`mkdir -p $callpeak` if (!-d $callpeak);
#export R
# `export PATH=$PATH:/public/home/zhangfenghua/miniconda3/envs/R/bin`;

open O , "> 01.$name.fastp.sh";
open F , "> 02.$name.bowtie2.sh";
open K , "> 03.$name.picard.sh";
open L , "> 04.$name.bamCoverage.sh";
open P , "> 05.$name.bedtools.sh";
open D , "> 06.$name.callPeak.sh";
open A , "> 07.$name.fragmentsLength.sh";
open I , "> 08.$name.computeMatrix.sh";
open U , "> 09.$name.plotHeatmap.sh";

foreach my $fq_1 (@fastq_1s){
	my $fq_2 = $fq_1 ;
	$fq_2 =~ s/_1\.$postfix$/_2\.$postfix/;
	$fq_1 =~ /([^\/]+)\/([^\/]+)_1\.$postfix$/;
	my ($tmpdir, $basename) = ($1, $2) ;
	my $tmpoutdir = "$outdir/$tmpdir";
	my $quaCondir = "$qualiControl/$tmpdir";
	my $callpeakdir = "$callpeak/$tmpdir";
	`mkdir -p $tmpoutdir` if (!-d $tmpoutdir);
	`mkdir -p $quaCondir` if (!-d $quaCondir);
	`mkdir -p $callpeakdir` if (!-d $callpeakdir);
	print O "$fastp -i $fq_1 -I $fq_2 -o $tmpoutdir/${basename}_1.fq.gz -O $tmpoutdir/${basename}_2.fq.gz -h $tmpoutdir/${basename}.html $fastpCMD\n";
	#/public/home/xuwenjie/software/chromap/chromap-master/chromap
	print F "$bowtie --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 40 -x $refIdx -1 $tmpoutdir/${basename}_1.fq.gz -2 $tmpoutdir/${basename}_2.fq.gz 2> $tmpoutdir/${basename}.bowtie.log | samtools view -\@10 -F 1804 -h -f 2 -q 25 | samtools sort -\@10 -m 2G -O BAM -o $tmpoutdir/${basename}.bam\n";
	print K "$picard INPUT=$tmpoutdir/${basename}.bam OUTPUT=$tmpoutdir/${basename}.rmdup.bam METRICS_FILE=$tmpoutdir/${basename}.rmdup.bam.dupqc.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true; samtools index -\@10 $tmpoutdir/${basename}.rmdup.bam\n";
	print L "$bamCoverage -b $tmpoutdir/${basename}.rmdup.bam -o $tmpoutdir/${basename}.rmdup.bw -of bigwig --binSize 10 --normalizeUsing RPKM --extendReads \n";
	print P "$bedtools bamtobed -i $tmpoutdir/${basename}.rmdup.bam > $tmpoutdir/${basename}.rmdup.bed \n";
	print D "$macs3 callpeak -t $tmpoutdir/${basename}.rmdup.bam -n $basename --outdir $callpeakdir -B -f BAMPE --broad --broad-cutoff 0.1 -g 1070912093\n"; # -g 1070912093 为鸡基因组大小
	print A "$gatk  -H $quaCondir/${basename}.fragmentslength.pdf -M 0.5 -I $tmpoutdir/${basename}.rmdup.bam -O $quaCondir/${basename}.InsertSize.txt \n";
	print I "computeMatrix reference-point --referencePoint TSS -S $tmpoutdir/${basename}.rmdup.bw -R $tssref -b 3000 -a 3000 -out $quaCondir/${basename}.scale_regions.tab.gz --skipZeros --missingDataAsZero\n";
	print U "plotHeatmap -m $quaCondir/${basename}.scale_regions.tab.gz -out $quaCondir/${basename}.plotHeatmap.png --legendLocation none\n";
}
close O ;
close F ;
close K ;
close L ;
close P ;
close D ;
close A ;
close I ;
close U ;
