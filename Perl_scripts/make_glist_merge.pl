use strict;
use warnings;

my %keepG;

my @args = (0,1,2);

foreach(@args)
{
	my $i = $_;
	
	open DVD, $ARGV[$i] or die $!;

	while(<DVD>)
	{
		chomp;
		my @sps = split(/,/);
		
		$keepG{$sps[0]} = 1;
	}

	close DVD;
}

open BLAST, $ARGV[3] or die $!;

while(<BLAST>)
{
	chomp;
	
	my @sps = split(/\t/);
	
	my @sp = split(/\|/, $sps[1]);
	
	if(exists $keepG{$sps[0]})
	{
		print $sp[1] . "\n";
	}
}

close BLAST;