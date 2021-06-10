#!/usr/bin/perl
########################################################
# nimarafati@gmail.com	                               #
# Please cite the script by referencing to github      #
# repository 					                                 #
########################################################
use Getopt::Long;
my $minD=5;
my $maxD_Norm=45;
my $maxD_Filt=100;
my $GQ=0;
my $prop_PASS=0.8;

my $contact="Nima Rafati nimarafati\@gmail.com
V1	20170531	";
my $usage = "This script perform chi-squre test between two groups by providing following info:
1- File including AD generated from GATK:VariantsToTable (Please make sure to use PASS SNPs where you have exluded missing calls. For this use \"--maxNoCall 0\" while running GATK:VariatsToTable)
2- Two files listing samples to compare (one sample per line)
3- Depth info

Chisq_from_AD_V1.script.pl -i AD.file -a group1.txt -b group2.txt -minD 5 -maxD-to-Normalize 45 -maxD-to-Filter 100
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
'GQ=i' =>\$GQ);
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
}
close inFGroup1;
#Group2
open(inFGroup2,$group2);
while(<inFGroup2>)
{
        chomp($_);
        $group2{$_}=0;
}
close inFGroup2;
my $group1_Size=keys %group1;
my $group2_Size=keys %group2;
##
##
my @lineArr=();
my $cntr=0;
my $order=0;
open(outF1, ">${inputFile}_tmp.bed");
print outF1 "Chr	Start	End	Ref_Freq_1	Ref_Freq_2	Prop_Samples_1	Prop_Samples_2	Chisq_Sum\n";
open(inF1,$inputFile) || die print $usage;
while(<inF1>)
{
	$cntr++;
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
				if(exists $group1{$id})
				{
					$group1{$id}=$order;
				}
				if(exists $group2{$id})
				{
					$group2{$id}=$order;
				}
			}
			if($lineArr[$i]=~ m/(.*).GQ/)
			{
				$id=$1;
				$id_GQ="$id.GQ";
				$ind_GQ{$id}=$order;
			}
			$order++;
		}
	}
	if ($cntr>=1 && $lineArr[0] ne "CHROM")
	{
		chomp($_);
		foreach my $k (keys %group1)
		{
			$tmp_R=0;$tmp_V=0;
			$lineArr[$group1{$k}]=~ m/(\d+),(\d+)/;
			$tmp_R=$1;$tmp_V=$2;
			if($tmp_R+$tmp_V >= $minD && $tmp_R+$tmp_V <= $maxD_Filt && $lineArr[$ind_GQ{$k}]>=$GQ)
			{
				if($tmp_R+$tmp_V>=$maxD_Norm)
				{
					$cntr_1++;
					$group1_R_Obs+=int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V));
					$group1_V_Obs+=int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V));
				}
				else
				{
					$cntr_1++;
					$group1_R_Obs+=$tmp_R;
					$group1_V_Obs+=$tmp_V;
				}
			}
		}
		foreach my $k (keys %group2)
		{
			$tmp_R=0;$tmp_V=0;
			$lineArr[$group2{$k}]=~ m/(\d+),(\d+)/;
			$tmp_R=$1;$tmp_V=$2;
			if($tmp_R+$tmp_V >= $minD && $tmp_R+$tmp_V <= $maxD_Filt  && $lineArr[$ind_GQ{$k}]>=$GQ)
			{
				if($tmp_R+$tmp_V>=$maxD_Norm)
				{
					$cntr_2++;
					$group2_R_Obs+=int(($tmp_R*$maxD_Norm)/($tmp_R+$tmp_V));
					$group2_V_Obs+=int(($tmp_V*$maxD_Norm)/($tmp_R+$tmp_V));
				}
				else
				{
					$cntr_2++;
					$group2_R_Obs+=$tmp_R;
					$group2_V_Obs+=$tmp_V;
				}
			#	print "$_\n$k\t$group2{$k}\t$lineArr[$group2{$k}]==>$tmp_R,$tmp_V\t$group2_R_Obs,$group2_V_Obs";<STDIN>;
			}
		}

		#Statistical analysis
		if($cntr_1/$group1_Size>=0.8 && $cntr_2/$group2_Size>=0.8)
		{
			#group1
			$prop1_Sample=sprintf("%.2f",$cntr_1/$group1_Size);
			$group1_All=$group1_R_Obs+$group1_V_Obs;
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

			print outF1 "$chr\t$start\t$end\t$group1_R_Freq\t$group2_R_Freq\t$prop1_Sample\t$prop2_Sample\t$sum_Chi\n";
		}
	}
}
close inF1;
close outF1;
open(outF2,">${inputFile}_tmp.Rscript");
print outF2 "
tmp<-read.table(\"${inputFile}_tmp.bed\",header=T)
chisq_result<-pchisq(tmp\$Chisq_Sum,1,lower.tail=F)
tmp<-cbind.data.frame(tmp,log10_P_Value=-log10(chisq_result))
pdf(\"${inputFile}_chisq_${group1}_${group2}.pdf\",width=10,height=6)
plot(tmp\$log10_P_Value,col=as.factor(tmp\$Chr),main=\"chisq_${group1}_${group2}\",ylab=\"-log10(P_Value\")
dev.off()
write.table(tmp,\"${inputFile}_chisq.bed\",col.names=T,row.names=F,quote=F,sep=\"\t\")";
system(`R --vanilla -q < ${inputFile}_tmp.Rscript`);
