use strict;
use warnings;
use Getopt::Long;

my($genofile,$phenofile,$locus,$msn,$outfile);

my @phenos = ();

my $usage = "USAGE:\nperl $0 --gfile <genotype file> --pfile <phenotype file> --locus <locus> --pheno <phenotype name> --msn <minor sample number of haplotype> --out <outfile>\n";
$usage .= "<genotype file> is the genotype file. Tab-limited. Format: sample marker1 marker2 ...\n";
$usage .= "<phenotype file> is the phenotype file. Tab-limited. Format: sample pheno1 pheno2 ...";
$usage .= "<locus> locus name.\n";
$usage .= "<msn> minor sample number of a haplotype to be picked. Default=10\n";
$usage .= "<outfile> output file.\n";

GetOptions(
	"gfile=s" => \$genofile,
	"pfile=s" => \$phenofile,
	"locus=s" => \$locus,
	"msn=s" => \$msn,
	"pheno=s" => \@phenos,
	"out=s" => \$outfile,
) or die $usage;

die $usage unless(defined $genofile and defined $phenofile and defined $locus and @phenos > 0 and defined $outfile);

unless(defined $msn){
	$msn = 10;
}

my %hash_sample;
my %hash_geno;
my %hash_pheno = map {$phenos[$_], ""} 0..$#phenos;

my @tmps = `head -1 $genofile|sed 's/\\t/\\n/g'`;
my $locusrank;
for(my $i = 0; $i < @tmps; $i++){
	chomp($tmps[$i]);
	if($tmps[$i] eq $locus){
		$locusrank = $i+1;
	}
}
unless(defined $locusrank){
	die "# WARNING: locus:$locus does not exist in the genotype file.\n";
}

open(IN,"cat $genofile|awk '{print \$1\"\\t\"\$$locusrank}'|") or die $!;
while(<IN>){
	chomp;
	my($sample,$haprank) = split/\t/;
	$sample =~ s/Sample_//;
	$sample =~ s/Line//;
	if($. == 1){
		die unless($haprank eq $locus);
		next;
	}
	push @{$hash_geno{$haprank}{samples}},$sample;
	$hash_sample{$sample}{haprank} = $haprank;
}
close IN;

open(IN,"<$phenofile") or die $!;
open(OUT,">$outfile");

my @colranks = ();

while(<IN>){
	chomp;
	my($sample,@phenonames) = split/\t/;
	if($. == 1){
		my $outhead = "sample\thaprank";
		for(my $i = 0; $i < @phenonames; $i++){
			if(exists $hash_pheno{$phenonames[$i]}){
				$hash_pheno{$phenonames[$i]} = $i;
				push @colranks, $i;
				$outhead .= "\t$phenonames[$i]";
			}
		}
		for(my $i = 0; $i < @phenos; $i++){
			if($hash_pheno{$phenos[$i]} eq ""){
				die "# WARNING: phenotype:$phenos[$i] does not exist in the phenotype file. Please check.\n";
			}
		}
		print OUT "$outhead\n";
		next;
	}
	$sample =~ s/Sample_//;
	$sample =~ s/Line//;
	my $outline = "$sample";
	next unless(exists $hash_sample{$sample});
	my $haprank = $hash_sample{$sample}{haprank};
	if(@{$hash_geno{$haprank}{samples}} < $msn or $haprank eq "-"){
		$haprank = "other";
	}
	$outline .= "\t$haprank";
	foreach(@colranks){
		$outline .= "\t$phenonames[$_]";
	}
	print OUT "$outline\n";
}
close IN;
close OUT;

# ploting
system("Rscript /public/share/zhuochen/rrBLUP_GWAS/box_beeswarm.r $outfile $outfile");