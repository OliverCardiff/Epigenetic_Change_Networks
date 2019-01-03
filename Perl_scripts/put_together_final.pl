#usage: perl put_together_final.pl <annolist.txt> <genfx.blast> <deseq_out_gene.txt> <deseq_out_meth.txt> <miRNA_net.txt> <deseq_out_miRNA> <cc_david_chart> <pval>

use strict;
use warnings;

my %annos;

my %keepG;
my %keepMi;

my %geneMeth;
my %geneFC;
my %miFC;

my $pval = $ARGV[7];

open ALIST, $ARGV[0] or die $!;

while(<ALIST>)
{
	chomp;
	
	#chop($_);
	
	$annos{$_} = 1;
}

close ALIST;

my %blHash;
my %blHash_a;


open BLAST, $ARGV[1] or die $!;

while(<BLAST>)
{
	chomp;
	
	my @sps = split(/\t/);
	
	my @sp = split(/\|/, $sps[1]);
	
	$blHash{$sps[0]} = $sp[1];
	push(@{$blHash_a{$sp[1]}}, $sps[0]);
}

close BLAST;

open DE_G, $ARGV[2] or die $!;

while(<DE_G>)
{
	chomp;
	my $line = $_;
	chop($line);
	my @sps = split(/\t/);
	if($. > 1)
	{
		if(exists $blHash{$sps[0]})
		{
			my $bl = $blHash{$sps[0]};
			if(defined $annos{$bl})
			{
				if($sps[5] ne "NA" && $sps[5] < $pval)
				{
					$keepG{$sps[0]} = 1;
					$geneFC{$sps[0]} = $sps[2];
				}
			}
		}
	}
}

close DE_G;

open DE_M, $ARGV[3] or die $!;

while(<DE_M>)
{
	chomp;
	my $line = $_;
	chop($line);
	my @sps = split(/\t/);
	if($. > 1)
	{
		if(exists $keepG{$sps[0]})
		{
			$geneMeth{$sps[0]} = $sps[2];
		}
	}
}

close DE_M;

open NET, $ARGV[4] or die $!;
#open my $minOut, '>', "miRNA_net.csv" or die $!;

#print {$minOut} "Source,Attri,Target\n";

while(<NET>)
{
	chomp;
	
	my @sps = split(/,/);
	
	chop($sps[1]);
	
	if(defined $keepG{$sps[1]})
	{
		#$keepMi{$sps[0]} = 1;
		my $st = $sps[0] . ",2," . $sps[1] . "," . "none" . "\n";
		push(@{$keepMi{$sps[0]}}, $st);
	}
}

close NET;

open DE_R, $ARGV[5] or die $!;

while(<DE_R>)
{
	chomp;
	my $line = $_;
	chop($line);
	my @sps = split(/\t/);
	if($. > 1)
	{
		if(exists $keepMi{$sps[0]})
		{
			$miFC{$sps[0]} = $sps[2];
		}
	}
}

close DE_R;

open ANNO, $ARGV[6] or die $!;
open my $fout, '>', "compartment_net.csv" or die $!;
print {$fout} "Source,Attri,Target,detail\n";

while(<ANNO>)
{
	chomp;
	
	if($. > 1 && $. < 5000)
	{
		my @sps = split(/\t/);
		my $ll = scalar @sps;
		if($ll > 5)
		{
			my @sp2 = split(/\s/, $sps[5]);
			my @genList;
			
			foreach(@sp2)
			{
				my $ano = $_;
				$ano =~ tr/,//d;
				
				if(exists $blHash_a{$ano})
				{
					#print "here $ano\n";
					my @genes = @{$blHash_a{$ano}};
					
					my $ln = scalar @genes;
					
					if($ln > 0)
					{
						push(@genList, @genes);
					}
				}
			}
			
			my $ln = scalar @genList;
			for(my $i = 0; $i < $ln-1; $i++)
			{
				for(my $j = $i+1; $j < $ln; $j++)
				{
					my $g1 = $genList[$i];
					my $g2 = $genList[$j];
					
					if(exists $keepG{$g1} && exists $keepG{$g2})
					{
						print {$fout} $g1 . ",1," . $g2 . "," . $sps[1] . "\n";
					}
				}
			}
		}
	}
}

close ANNO;

open my $nodeOut, '>', "Node_data.csv" or die $!;

print {$nodeOut} "Node,lfcSE,MethFC,Sort\n";

my @keyset = keys %keepG;

foreach(@keyset)
{
	my $gn = $_;
	my $gfc = $geneFC{$gn};
	my $mfc = $geneMeth{$gn};
	if($mfc eq "NA")
	{
		$mfc = 0;
	}
	my $ty = 1;
	
	print {$nodeOut} $gn . "," . $gfc . "," . $mfc . "," . $ty . "\n";
}

my @keys2 = keys %keepMi;

foreach(@keys2)
{
	my $mi = $_;
	if(exists $miFC{$mi})
	{
		my $gfc = $miFC{$mi};
		my $mfc = 0;
		my $ty = 2;
		
		print {$nodeOut} $mi . "," . $gfc . "," . $mfc . "," . $ty . "\n";
		my @entries = @{$keepMi{$mi}};
		foreach(@entries)
		{
			print {$fout} $_;
		}
	}
}

close $nodeOut;
close $fout;
#close $minOut;