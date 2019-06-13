use strict;
use warnings;


my($input,$outfile,$threshold) = @ARGV;
my $usage = "USAGE:\nperl $0 <input fasta> <output fasta> <threshold>\n";
$usage .= "Slices shorter than <threshold> will be thrown. Default=100\n";

die $usage unless(@ARGV >= 2);
unless(defined $threshold){
	$threshold = 100; #change here
}

open(FA,"<$input") or die "$!";

local $/ = "\n>";
my %hash_seq;
while(<FA>){
	$_ =~ s/^>//g;
	chomp;
	my($seqhead,$seq) = split/\n/,$_,2;
	my ($seqname) = $seqhead =~ /^(\S+)/;
	$seq =~ s/\s+//g;
	$hash_seq{$seqname}{seq} = $seq;
	#print $seqname."\t".length($seq)."\n";
}
close FA;
local $/ = "\n";

open(OUT,">$outfile");
my $splitcount = 0;
my $throwcount = 0;
foreach my $seqname(keys %hash_seq){
	my $seq = $hash_seq{$seqname}{seq};
	$seq = uc($seq);
	$seq =~ s/^(N+)//ig;
	$seq =~ s/N+/\|/ig;
	my @seqs = split/\|/,$seq;
	my $count = 1;
	foreach my $cut (@seqs){
		if(length($cut) < $threshold){
			$throwcount ++;
			#print length($cut)."\n";
			next;
		}
		print OUT ">$seqname:$count\n";
		#print "$seqname:$count\t".length($cut)."\n";
		$cut =~ s/(.{50})/$1\n/g;
		print OUT "$cut\n";
		$count++;
		$splitcount ++;
	}
}
close OUT;
my $inputseqnum = keys %hash_seq;
print "Input sequence num: $inputseqnum\n";
print "Output sequence num: $splitcount\n";
print "Small cuts less than $threshold bases are thrown: $throwcount\n";
