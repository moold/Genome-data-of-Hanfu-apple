#!/usr/bin/perl -w
use strict;

use Getopt::Long;


my ($Genome,@Gff,$Threshold,$Help,%seq,%rep);

GetOptions(
	"genome=s"=>\$Genome,
	"gff:s"=>\@Gff,
	"help"=>\$Help,
	"threshold:s"=>\$Threshold
);

die "perl $0 -genome *.fa -gff *.gff -gff *.gff -threshold 500 \n" if ((not defined $Genome) || $Help);

$Threshold ||= 0; 


Read_fasta($Genome,\%seq);

print STDERR "read sequence done\n";

Read_gff(\@Gff,\%rep);

print STDERR "read gff done\n";

foreach my $seq_name (sort keys %seq) {
	my $seq_head = $seq{$seq_name}{head};
	my $seq_str = $seq{$seq_name}{seq};	
	my @pos;
	if (exists $rep{$seq_name}){
		@pos=@{$rep{$seq_name}};
		Remask_seq(\$seq_str,\@pos);
	}
	Display_seq(\$seq_str);
	print ">$seq_head\n$seq_str";
}

print STDERR "Mask task done\n";


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name}{head} =  $head;
		$hash_p->{$name}{len} = length($seq);
		$hash_p->{$name}{seq} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}


##read repeat gff file
#usage: Read_gff($file,\%hash);
############################################
sub Read_gff{
	my $file=shift;
	my $hash_p=shift; 
	for my $infile (@{$file}){
		open (IN,$infile) || die ("fail open $infile\n");
		while (<IN>) {
			chomp;
			s/^\s+//;
			next if (/^#/);
			my (@t,$tname,$start,$end);
			@t=split;
			$tname = $t[0];
			($start,$end)=($t[3]<$t[4])?( $t[3],$t[4]):($t[4],$t[3]);
			next if ($end-$start+1< $Threshold);
			push @{$hash_p->{$tname}},[$start,$end];
		}
		close(IN);
	}
}




#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = $disp;
}
#############################################


##remask a sequence with repeat positon 
##usage: Remask_seq(\$seq,\@pos);
###########################################
sub Remask_seq{
	my $seq_p = shift; ##sequence pointer
	my $rep_ap = shift; ##array pointer for repeat position
	my $maskC = 'n'; ##set the base for masking result
	
	$$seq_p =~ s/\s//g;
	foreach my $p (@$rep_ap) {
		my ($start,$end) = ($p->[0] <= $p->[1]) ? ($p->[0] , $p->[1]) : ($p->[1], $p->[0]);
		substr($$seq_p,$start-1,$end-$start+1) = $maskC x ($end-$start+1) if($start-1<=length $$seq_p);
	}
}
