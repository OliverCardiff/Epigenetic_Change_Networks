use strict;
use warnings;

my %clusters;

open ANNOT, $ARGV[0] or die $!;

my $firstHeader = 0;
my $inc = 0;

my %fileHands;

my %clNams;
my %clLens;

while(<ANNOT>)
{
	chomp;
	
	my @sps = split(/\t/);
	
	my $ln = scalar @sps;
	
	if($ln == 2)
	{
		$firstHeader = 1;
		$inc++;
	}
	elsif($firstHeader != 1 && $ln > 1 && $inc < 21)
	{
		my @sp = split(/,/, $sps[5]);
		
		if(!exists $fileHands{$inc} && $inc < 21)
		{
			$clNams{$inc} = $sps[1];
			$clLens{$inc} = $sps[2];
			open my $fout, '>', "clust_" . $inc . ".txt" or die $!;
			
			$fileHands{$inc} = $fout;
		}
		
		foreach(@sp)
		{
			my $gn = $_;
			if(length($gn) == 7)
			{
				$gn = substr $gn, 1;
			}
			$clusters{$inc}{$gn} = 2;
		}
	}
	else
	{
		$firstHeader = 0;
	}
}

close ANNOT;

my %blHash;

open BLAST, $ARGV[1] or die $!;

while(<BLAST>)
{
	chomp;
	
	my @sps = split(/\t/);
	
	my @sp = split(/\|/, $sps[1]);
	
	$blHash{$sps[0]} = $sp[1];
}

close BLAST;


my %changeHash;
open GLIST, $ARGV[2] or die $!;

while(<GLIST>)
{
	chomp;
	
	my $ln = substr $_, 0, (length($_)-1);
	
	if(exists $blHash{$ln})
	{
		my $gn = $blHash{$ln};
		
		for(my $i = 1; $i <= 20; $i++)
		{
			if(defined $clusters{$i}{$gn})
			{
				if($clusters{$i}{$gn} > 0)
				{
					$clusters{$i}{$gn}--;
					print {$fileHands{$i}} $gn . "\n";
				}
			}
		}
	}
	
}
close GLIST;
=cut