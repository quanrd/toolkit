use strict;
use warnings;

my($bam,$bed,$out) = @ARGV;
my $usage = "USAGE:\nperl $0 <bam> <bed> <out>\n";
$usage .= "<bam> is bam file sorted by position.\n";
$usage .= "<bed> is tab-limited bed file. format: chr start end tag\n";

die $usage unless(@ARGV == 3);

my %hash_bed;

open(IN,"<$bed") or die $!;
while(<IN>){
	chomp;
	my($chr,$start,$end,$tag) = split/\t/;
	$hash_bed{$chr}{reg}{$start}{$end}{$tag} = "";
}
close IN;

foreach my $chr(keys %hash_bed){
	foreach my $start(sort {$a <=> $b} keys %{$hash_bed{$chr}{reg}}){
		foreach my $end(sort {$a <=> $b} keys %{$hash_bed{$chr}{reg}{$start}}){
			push @{$hash_bed{$chr}{arr}}, "$start,$end";
		}
	}
}

open(IN,"samtools view $bam -L $bed|") or die $!;
my %hash_tag;

while(<IN>){
	chomp;
	my($name,undef,$chr,$ctgstart,undef,$cigar,@others) = split/\t/;
	my $len = range_alignment($cigar);
	my $ctgend = $ctgstart + $len - 1;
	next unless(exists $hash_bed{$chr});
	for(my $i = 0; $i < @{$hash_bed{$chr}{arr}}; $i++){
		my($start,$end) = split/,/,${$hash_bed{$chr}{arr}}[$i];
		if($start >= $ctgend){
			last;
		}
		if($end < $ctgstart){
			splice(@{$hash_bed{$chr}{arr}}, $i ,1);
			$i--;
			next;
		}
		my @tags = sort keys %{$hash_bed{$chr}{reg}{$start}{$end}};
		foreach my $tag(@tags){
			$hash_tag{$tag}{ctg}{$name}++;
		}
	}
}
close IN;

open(OUT,">$out");
foreach my $tag(sort keys %hash_tag){
	my @ctgs = sort keys %{$hash_tag{$tag}{ctg}};
	print OUT "$tag\t".join("\t",@ctgs)."\n";
}
close OUT;

sub range_alignment{
	my($CIGAR) = @_;
	my @arr = $CIGAR =~ /(\d+[A-Z])/g;
	my $len = 0;
	foreach(@arr){
		my($num,$symbol) = $_ =~ /(\d+)(\S+)/;
		if($symbol eq "M" or $symbol eq "D"){
			$len += $num;
		}
	}
	return($len);
}