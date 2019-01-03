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
	elsif($firstHeader != 1 && $ln > 1)
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
			$clusters{$inc}{$gn} = 1;
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
	#print $sps[0] . "\t" . $sp[1] . "\n";
}

close BLAST;

my %changeHash;

open GLIST, $ARGV[2] or die $!;

while(<GLIST>)
{
	chomp;
	
	my $ln = substr $_, 0, (length($_)-1);
	
	$changeHash{$ln} = 1;
	
}

close GLIST;

open MATRIX, $ARGV[3] or die $!;

while(<MATRIX>)
{
	chomp;
	my $line = $_;
	my @sps = split(/\t/);
	
	if(exists $blHash{$sps[0]})
	{
		my $gen = $blHash{$sps[0]};
		
		if(exists $changeHash{$sps[0]})
		{
			for(my $i = 1; $i < 21; $i++)
			{
				print $i . "\t" . $clNams{$i} . "\t" . $clLens{$i} . "\n";
				if(exists $clusters{$i}{$gen})
				{
					my $hand = $fileHands{$i};
					
					print {$hand} $line . "\n";
				}
			}
		}
	}
}

close MATRIX;

for(my $i = 1; $i < 21; $i++)
{
	close $fileHands{$i};
}