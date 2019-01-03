use strict;
use warnings;

my %miIntra;
my %meIntra;
my %geIntra;

my $stub = $ARGV[0];

my $clcnt = 0;

my $oldFLD = 0;

print "rn,1,2,3,4,5,6,7,8,9,10\n";

open MI_IN, "mi_" . $stub . "_clust.txt" or die $!;

print "miRNA";

while(<MI_IN>)
{
	chomp;
	my @sps = split(/\t/);
	my $ln = scalar @sps;
	
	if($ln > 1)
	{
		my @spX = split(/\s+/, $sps[1]);
		
		if(defined $spX[0] && $spX[0] eq "Enrichment" && $clcnt < 10)
		{
			$clcnt++;
			$oldFLD = $spX[2];
			print "," . $oldFLD;
		}
	}
}

close MI_IN;

print "\nmeth";

$clcnt = 0;

open ME_IN, "me_" . $stub . "_clust.txt" or die $!;

while(<ME_IN>)
{
	chomp;
	my @sps = split(/\t/);
	my $ln = scalar @sps;
	
	if($ln > 1)
	{
		my @spX = split(/\s+/, $sps[1]);
		
		if(defined $spX[0] && $spX[0] eq "Enrichment" && $clcnt < 10)
		{
			$clcnt++;
			$oldFLD = $spX[2];
			print "," . $oldFLD;
		}
	}
}

close ME_IN;

print"\ngene";

$clcnt = 0;

open GE_IN, "g_" . $stub . "_clust.txt" or die $!;

while(<GE_IN>)
{
	chomp;
	my @sps = split(/\t/);
	my $ln = scalar @sps;
	
	if($ln > 1)
	{
		my @spX = split(/\s+/, $sps[1]);
		
		if(defined $spX[0] && $spX[0] eq "Enrichment" && $clcnt < 10)
		{
			$clcnt++;
			$oldFLD = $spX[2];
			print "," . $oldFLD;
		}
	}
}

close GE_IN;

print "\n";

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}