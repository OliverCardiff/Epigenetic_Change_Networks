use strict;
use warnings;

my %miIntra;
my %meIntra;
my %geIntra;

my $stub = $ARGV[0];
my $limit = $ARGV[1];

my @tList = ();

my $clcnt = 0;

my $oldFLD = 0;

open MI_IN, "mi_" . $stub . "_clust.txt" or die $!;
open my $netOut, '>', $stub . "_net.csv" or die $!;
open my $nodeOut, '>', $stub . "_node.csv" or die $!;
print {$netOut} "Source,Target,Attri\n";
print {$nodeOut} "Node,Label,Count,Fold,Typeof,Desc,Log10\n";

while(<MI_IN>)
{
	chomp;
	my $subs = substr $_,0,1;
	if($_ eq "" || $subs eq "\t")
	{
		$clcnt++;
		
		if($clcnt < $limit)
		{
			my $ln = scalar @tList;
			
			for(my $i = 0; $i < $ln-1; $i++)
			{
				for(my $j = $i + 1; $j < $ln; $j++)
				{
					my $t1 = $tList[$i]; my $t2 = $tList[$j];
					
					print {$netOut} "mi_" . $t1 . ",mi_" . $t2 . ",2\n";
				}
			}
		}
		@tList = ();
	}
	else
	{
		my @sps = split(/\t/);
		
		my @spX = split(/\s+/, $sps[1]);
		
		if($spX[0] eq "Enrichment")
		{
			$oldFLD = $spX[2];
		}
		
		my $lx = scalar @sps;
		if($lx > 3 && $sps[0] ne "Category" && $sps[2] ne "" && $clcnt < $limit - 1)
		{
			my @sp2 = split(/~/, $sps[1]);
			
			push(@tList, $sp2[0]);
			$sp2[1] =~ tr/,//d;
			my $nlog = log10($sps[4]) * -1;
			print {$nodeOut} "mi_" . $sp2[0] . ",miRNA," . $sps[2] . "," . $oldFLD . ",1," . $sp2[1] . "," . $nlog ."\n";
			
			$miIntra{$sp2[0]} = 1;
		}
	}
	
}

close MI_IN;

@tList = ();

$clcnt = 0;

open ME_IN, "me_" . $stub . "_clust.txt" or die $!;

while(<ME_IN>)
{
	chomp;
	my $subs = substr $_,0,1;
	
	if($_ eq "" || $subs eq "\t")
	{
		$clcnt++;
		
		if($clcnt < $limit)
		{
			my $ln = scalar @tList;
			
			for(my $i = 0; $i < $ln-1; $i++)
			{
				for(my $j = $i + 1; $j < $ln; $j++)
				{
					my $t1 = $tList[$i]; my $t2 = $tList[$j];
					
					print {$netOut} "me_" . $t1 . ",me_" . $t2 . ",2\n";
				}
			}
		}
		@tList = ();
	}
	else
	{
		my @sps = split(/\t/);
		
		my @spX = split(/\s+/, $sps[1]);
		
		if($spX[0] eq "Enrichment")
		{
			$oldFLD = $spX[2];
		}
		
		my $lx = scalar @sps;
		if($lx > 3 && $sps[0] ne "Category" && $sps[2] ne "" && $clcnt < $limit - 1)
		{
			my @sp2 = split(/~/, $sps[1]);
			
			push(@tList, $sp2[0]);
			$sp2[1] =~ tr/,//d;
			my $nlog = log10($sps[4]) * -1;
			print {$nodeOut} "me_" . $sp2[0] . ",MeDIPs," . $sps[2] . "," . $oldFLD  . ",2," . $sp2[1] . "," . $nlog ."\n";
			
			$meIntra{$sp2[0]} = 1;
		}
	}
	
}

close ME_IN;

@tList = ();

$clcnt = 0;

open GE_IN, "g_" . $stub . "_clust.txt" or die $!;

while(<GE_IN>)
{
	chomp;
	my $subs = substr $_,0,1;
	if($_ eq "" || $subs eq "\t")
	{
		$clcnt++;
		
		if($clcnt < $limit)
		{
			my $ln = scalar @tList;
			
			for(my $i = 0; $i < $ln-1; $i++)
			{
				for(my $j = $i + 1; $j < $ln; $j++)
				{
					my $t1 = $tList[$i]; my $t2 = $tList[$j];
					
					print {$netOut} "ge_" . $t1 . ",ge_" . $t2 . ",2\n";
				}
			}
		}
		@tList = ();
	}
	else
	{
		my @sps = split(/\t/);
		
		my @spX = split(/\s+/, $sps[1]);
		
		if($spX[0] eq "Enrichment")
		{
			$oldFLD = $spX[2];
		}
		
		my $lx = scalar @sps;
		if($lx > 3 && $sps[0] ne "Category" && $sps[2] ne "" && $clcnt < $limit - 1)
		{
			my @sp2 = split(/~/, $sps[1]);
			
			push(@tList, $sp2[0]);
			$sp2[1] =~ tr/,//d;
			my $nlog = log10($sps[4]) * -1;
			print {$nodeOut} "ge_" . $sp2[0] . ",RNASeq," . $sps[2] . "," . $oldFLD  . ",3," . $sp2[1] . "," . $nlog . "\n";;
			
			$geIntra{$sp2[0]} = 1;
		}
	}
	
}

close GE_IN;

my @miKey = keys %miIntra;

foreach(@miKey)
{
	my $mi = $_;
	
	if(exists $geIntra{$mi})
	{
		print {$netOut} "mi_" . $mi . ",ge_" . $mi . ",1\n";
	}
	if(exists $meIntra{$mi})
	{
		print {$netOut} "mi_" . $mi . ",me_" . $mi . ",1\n";
	}
}

my @geKey = keys %geIntra;

foreach(@geKey)
{
	my $ge = $_;
	
	if(exists $meIntra{$ge})
	{
		print {$netOut} "ge_" . $ge . ",me_" . $ge . ",1\n";
	}
}

close $netOut;
close $nodeOut;

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}