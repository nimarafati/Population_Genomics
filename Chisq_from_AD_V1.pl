#!/usr/bin/perl
#Nima Rafati 20170531
#use warnings;
use Getopt::Long;
#use lib '~/PERL_LIB/lib/perl5';
#use Statistics::Descriptive;
#my $stat = Statistics::Descriptive::Full->new();

my $minD=5;
my $maxD_Norm=45;
my $maxD_Filt=100;
my $GQ=0;
my $prop_PASS=0.8;
my @tmp_valArr=();
my $contact="Nima Rafati nimarafati\@gmail.com
V1	20170531	";
my $usage = "This script perform chi-squre test between two groups by providing following info:
1- File including AD generated from GATK:VariantsToTable (Please make sure to use PASS SNPs where you have exluded missing calls. For this use \"--maxNoCall 0\" while running GATK:VariatsToTable)
2- Two files listing samples to compare (one sample per line)
3- Depth info

Chisq_from_AD_V1.script.pl -i AD.file -a group1.txt -b group2.txt -minD 5 -maxD 45
-i AD file (GATK VariantsToTable)
-a A file consisting of samples in group 1
-b A file consisting of samples in group 2
-minD minimum depth (default 5)
-maxD-to-Normalize maximum depth to normalize to(default 45) 
-maxD-to-Filter samples having depth larger than this will be excluded (default 100)
-GQ genotype quality cut-off. Please provide a cut-off if you have extracted genotype quality.
$contact\n";

&GetOptions('h=s' =>\$helpFlag,
'i=s' =>\$inputFile,
'a=s' => \$group1,
'b=s' =>\$group2,
'minD=i' =>\$minD,
'maxD-to-Normalize=i' =>\$maxD_Norm,
'maxD-to-Filter=i' =>\$maxD_Filt,
'GQ=i' =>\$GQ,
'prop-samples=f' =>\$prop_PASS);
if($helpFlag eq "h" || $inputFile eq "" || $group1 eq "" || $group2 eq "")
{
	print $usage;
	exit;
}
##Identify position of samples in vcf AD file
#Group1
open(inFGroup1,$group1);
while(<inFGroup1>)
{
	chomp($_);
	$group1{$_}=0;
	$group1_freq{$_}=-1;
}
close inFGroup1;
#Group2
open(inFGroup2,$group2);
while(<inFGroup2>)
{
        chomp($_);
        $group2{$_}=0;
	 $group2_freq{$_}=-1;
}
close inFGroup2;
my $group1_Size=keys %group1;
my $group2_Size=keys %group2;
##
##
my @lineArr=();
my $cntr=0;
my $order=0;
open(outF1, ">${inputFile}_${group1}_${group2}_tmp.bed");
print outF1 "Chr	Start	End	Ref_Freq_1	Ref_Freq_2	Prop_Samples_1	Prop_Samples_2	Chisq_Sum";
open(inF1,$inputFile) || die print $usage;
while(<inF1>)
{
	$cntr++;
#	print $cntr,"\n";
	chomp($_);
	@lineArr=split("\t",$_);
	###########Identify the position of samples in depth file
	if($lineArr[0] eq "CHROM")
	{
		foreach (my $i=0;$i < scalar(@lineArr);$i++)
		{
			#print "$lineArr[$i]";<STDIN>;
			if($lineArr[$i]=~ m/(.*).AD/)
			{
				$id=$1;
#				print "AD=$id\t$lineArr[$i]";<STDIN>;
#				print "$id\t$group1{$id}\t$group2{$id}";<STDIN>;
				if(exists $group1{$id})
				{
					$group1{$id}=$order;
#					print "1==> $group1{$id}\t$lineArr[$i]";
				}
				if(exists $group2{$id})
				{
					$group2{$id}=$order;
#					print "2==> $group2{$id}\t$id";
				}
			}
			if($lineArr[$i]=~ m/(.*).GQ/)
			{
				$id=$1;
#				print "GQ=$id\t$lineArr[$i]";<STDIN>;
				$id_GQ="$id.GQ";
				$ind_GQ{$id}=$order;
#				print "$id_GQ\t$ind_GQ{$id}\n";
			}
			$order++;
		}
		foreach my $k (keys %group1){
			print outF1 "\t$k";
		}
		foreach my $k (keys %group2){
			print outF1  "\t$k";
		}
		print outF1 "\tSD_group1\tSD_group2\n";
	}
	if ($cntr>=1 && $lineArr[0] ne "CHROM")
	{
		#Set SNP info
		$chr=$lineArr[0];
		$start=$lineArr[1]-1;
		$end=$lineArr[1];
		#Reseting reference and variant caount 
		$sum_Chi=0;
		$group1_R_Obs=0;$group1_V_Obs=0;$cntr_1=0;
		$group2_R_Obs=0;$group2_V_Obs=0;$cntr_2=0;
		#extracting reference and variant counts per group
		#-minD minimum depth (default 5)
		#-maxD-to-Normalize maximum depth to normalize to(default 45) 
		#-maxD-to-Filter samples having depth larger than this will be excluded (default 100)
		#-GQ genotype quality cut-off. Please provide a cut-off if you have extracted genotype quality.
		#'minD=i' =>\$minD,
		#'maxD-to-Normalize=i' =>\$maxD_Norm,
		#'maxD-to-Filter=i' =>\$maxD_Filt,
		#

		chomp($_);
		foreach my $k (keys %group1)
		{
			$tmp_R=0;$tmp_V=0;
			$lineArr[$group1{$k}]=~ m/(\d+),(\d+)/;
			$tmp_R=$1;$tmp_V=$2;
			if($tmp_R+$tmp_V >= $minD && $tmp_R+$tmp_V <= $maxD_Filt && $lineArr[$ind_GQ{$k}]>=$GQ)
			{
				$cntr_1++;
				$tmp_ref = int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V));
				$tmp_var = int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V));
				$group1_R_Obs+= $tmp_ref; #int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V));
				$group1_V_Obs+= $tmp_var; #int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V));
				$group1_freq{$k} = $tmp_ref/($tmp_ref+$tmp_var);
#print "tmp_ref: $tmp_ref
#tmp_var: $tmp_var
#$k: $group1_freq{$k}";<STDIN>;
				$tmp_ref = "";
				$tmp_var = "";
#				print "$k\t$lineArr[1]\t",int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V)),"\t",int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V)),"\n";
#				print "$k\t$lineArr[1]\t(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V))\t(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V))\n";
			}
		}
		foreach my $k (keys %group2)
		{
			$tmp_R=0;$tmp_V=0;
			$lineArr[$group2{$k}]=~ m/(\d+),(\d+)/;
			$tmp_R=$1;$tmp_V=$2;
			if($tmp_R+$tmp_V >= $minD && $tmp_R+$tmp_V <= $maxD_Filt  && $lineArr[$ind_GQ{$k}]>=$GQ)
			{
				$cntr_2++;
				$tmp_ref = int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V));
				$tmp_var = int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V));
				$group2_R_Obs+= $tmp_ref; #int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V));
				$group2_V_Obs+= $tmp_var; #int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V));
				$group2_freq{$k} = $tmp_ref/($tmp_ref+$tmp_var);
				$tmp_ref = "";
				$tmp_var = "";
#				print "$k\t$lineArr[$ind_GQ{$k}]\t$lineArr[1]\tR:$tmp_R\tV:$tmp_V\tmaxD_Norm:$maxD_Norm\ttmp_R+tmp_V:$tmp_R+$tmp_V==>",int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V)),"\t",int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V)),"\n";
#				print "$k\t$group2{$k}\t$lineArr[$group2{$k}]==>$tmp_R,$tmp_V\t$group2_R_Obs,$group2_V_Obs";<STDIN>;
			}
		}

		#Statistical analysis
#		print "$cntr_1/$group1_Size\t$cntr_2/$group2_Size\ngroup1:$group1_R_Obs / $group1_All\ngroup2:$group2_R_Obs / $group2_All";<STDIN>;
		if($cntr_1/$group1_Size>=0.8 && $cntr_2/$group2_Size>=0.8)
		{
			#group1
			$prop1_Sample=sprintf("%.2f",$cntr_1/$group1_Size);
			$group1_All=$group1_R_Obs + $group1_V_Obs;
			$group1_R_Freq=sprintf("%.3f",$group1_R_Obs/$group1_All);
			$group1_V_Freq=$group1_V_Obs/$group1_All;
			if($group1_R_Obs == 0)
			{
				$group1_R_Obs=1;
			}
			if($group1_V_Obs == 0)
			{
				$group1_V_Obs=1;
			}
			#group2
			$prop2_Sample=sprintf("%.2f",$cntr_2/$group2_Size);
			$group2_All=$group2_R_Obs+$group2_V_Obs;
			$group2_R_Freq=sprintf("%.3f",$group2_R_Obs/$group2_All);
			$group2_V_Freq=$group2_V_Obs/$group2_All;
			if($group2_R_Obs == 0)
			{
				$group2_R_Obs=1;
			}
			if($group2_V_Obs == 0)
			{
				$group2_V_Obs=1;
			}
#			print "$cntr_1/$group1_Size\t$cntr_2/$group2_Size==> ($group1_R_Obs+$group1_V_Obs == $grouop1_All  group1:$group1_R_Obs / $group1_All-------- $group2_R_Obs+$group2_V_Obs == $group2_All  $group2:$group2_R_Obs / $group2_All";<STDIN>;			

			#Calculating expeted freq
			$reference_Sum=$group1_R_Obs+$group2_R_Obs;
			$variant_Sum=$group1_V_Obs+$group2_V_Obs;
			$total_Sum=$reference_Sum+$variant_Sum;
			
			$group1_R_Exp=$group1_All*$reference_Sum/$total_Sum;
			$group1_V_Exp=$group1_All*$variant_Sum/$total_Sum;
			
			$group2_R_Exp=$group2_All*$reference_Sum/$total_Sum;
			$group2_V_Exp=$group2_All*$variant_Sum/$total_Sum;
		
			$group1_R_OE=($group1_R_Obs-$group1_R_Exp)**2/$group1_R_Exp;
			$group1_V_OE=($group1_V_Obs-$group1_V_Exp)**2/$group1_V_Exp;
			
			$group2_R_OE=($group2_R_Obs-$group2_R_Exp)**2/$group2_R_Exp;
			$group2_V_OE=($group2_V_Obs-$group2_V_Exp)**2/$group2_V_Exp;

			$sum_Chi=$group1_R_OE+$group1_V_OE+$group2_R_OE+$group2_V_OE;
#print "reference_Sum:$reference_Sum
#variant_Sum: $variant_Sum
#total_Sum: $total_Sum
#group1_R_Exp: $group1_R_Exp
#group1_V_Exp: $group1_V_Exp
#group2_R_Exp: $group2_R_Exp
#group2_V_Exp: $group2_V_Exp
#group1_R_OE: $group1_R_OE
#group1_V_OE: $group1_V_OE
#group2_R_OE: $group2_R_OE
#group2_V_OE: $group2_V_OE
#sum_Chi: $sum_Chi"; <STDIN>;
			print outF1 "$chr\t$start\t$end\t$group1_R_Freq\t$group2_R_Freq\t$prop1_Sample\t$prop2_Sample\t$sum_Chi";
			my @tmp_valArr = ();
			foreach my $k (keys %group1_freq){
				print outF1 "\t$group1_freq{$k}";
				push (@tmp_valArr, $group1_freq{$k});
			}
			$group1_stdev=std_dev(@tmp_valArr);
			my @tmp_valArr = ();
			foreach my $k (keys %group2_freq){
				print outF1 "\t$group2_freq{$k}";
				push (@tmp_valArr, $group2_freq{$k});
			}
			if(scalar @tmp_valArr >1){
				$group2_stdev=std_dev(@tmp_valArr);
			}else{
				$group2_stdev="NA";
			}
		}else{
			$group1_R_Freq= "NA"; #sprintf("%.3f",$group1_R_Obs/$group1_All);
			$group2_R_Freq= "NA"; #sprintf("%.3f",$group2_R_Obs/$group2_All);
			$prop1_Sample= sprintf("%.2f",$cntr_1/$group1_Size);
			$prop2_Sample= sprintf("%.2f",$cntr_2/$group2_Size);
			print outF1 "$chr\t$start\t$end\t$group1_R_Freq\t$group2_R_Freq\t$prop1_Sample\t$prop2_Sample\t0";
			foreach my $k (keys %group1_freq){
				print outF1 "\tNA";
			}
			foreach my $k (keys %group2_freq){
				print outF1 "\tNA";
			}
			$group1_stdev = "NA";
			$group2_stdev = "NA";
		}
		print outF1 "\t$group1_stdev\t$group2_stdev\n";
		$group1_stdev = ".";
		$group2_stdev = ".";
	}
}
close inF1;
close outF1;


### Subroutines
sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}
sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}
sub std_dev {
    my (@data) = @_;
    my ( $sq_dev_sum, $avg ) = ( 0, 0 );
    
    $avg = mean(@data);
    foreach my $elem (@data) {
        $sq_dev_sum += ( $avg - $elem )**2;
    }
    return ( sqrt( $sq_dev_sum / ( @data - 1 ) ) );
}


open(outF2,">${inputFile}_${group1}_${group2}_tmp.Rscript");
print outF2 "
library(data.table)
tmp<-fread(\"${inputFile}_${group1}_${group2}_tmp.bed\",header=T)
chisq_result<-pchisq(tmp\$Chisq_Sum,1,lower.tail=F)
tmp<-cbind.data.frame(tmp,log10_P_Value=-log10(chisq_result))
pdf(\"${inputFile}_chisq_${group1}_${group2}.pdf\",width=10,height=6)
plot(tmp\$log10_P_Value,col=as.factor(tmp\$Chr),main=\"chisq_${group1}_${group2}\",ylab=\"-log10(P_Value\")
dev.off()
write.table(tmp,\"${inputFile}_chisq_${group1}_${group2}.bed\",col.names=T,row.names=F,quote=F,sep=\"\\t\")";
system(`R --vanilla -q < ${inputFile}_${group1}_${group2}_tmp.Rscript`);
#system(`rm -f ${inputFile}_tmp.Rscript ${inputFile}_${group1}_${group2}_tmp.bed`)
