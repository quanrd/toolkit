use strict;
use warnings;

my($list1,$list2) = @ARGV;
my $usage = "USAGE:\nperl $0 <list1> <list2>\n";
$usage .= "Lists are tab-limit and the first column is used to match two lists(locus ID or someting else).\n";
$usage .= "A test.\n";

die $usage unless(@ARGV == 2);

my %hash;

open(IN,"<$list1") or die $!;
while(<IN>){
	chomp;
	my($locus,$other) = split/\t/, $_, 2;
	unless(defined $other){
		$other = "-";
	}
	$hash{$locus} = $other;
}
close IN;

open(IN,"<$list2") or die $!;
while(<IN>){
	chomp;
	my($locus,$other) = split/\t/, $_, 2;
	next unless(exists $hash{$locus});
	unless(defined $other){
		$other = "-";
	}
	print "$locus\t$hash{$locus}\t$other\n"
}
close IN;
