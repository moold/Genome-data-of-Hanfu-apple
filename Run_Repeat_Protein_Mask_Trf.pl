#!/usr/bin/perl -w
use strict;

=head1 ########################################USAGE#####################################################

=head1 Name

 Run_Repeat_Protein_Mask_Trf.pl  --The pipeline of finding tandem repeats and transposons in sequences.

=head1 DESCRIPTION

 This script include trf, RepeatMasker, and RepeatProteinMask.

 For RepeatMasker, you can either use inside TE libary by not seting -lib parameter,or give self made TE 
 libary by -lib option,but the given libary format should be  fasta and the sequence identifier must be less than 50 bp. 

 The total time depends on RepeatProteinMask and RepeatMasker.

 The final output will write to $outdir,the *.gff is the final result.


=head1 Version

 Author: Jiang Hu.
 Version: 1.0.

=head1 Usage

	perl Run_Repeat_Protein_Mask_Trf.pl -genome genome_path -trf -repeatMasker -repeatProteinMask -pvalue 1e-5 -lib *.fa	-cpu 5

	-genome <str>           set the genome path file,absolute paths,relative paths are allowed,needed.
	-trf                    set run Trf.
	-repeatMasker           set run RepeatMasker.
	-species <str>          set the RepeatMasker species,mammal, carnivore, rodentia, rat, cow, pig, et.al. default=all
	-repeatProteinMask      set run RepeatProteinMask.
	-engine                 set repeatMasker serch engine,[crossmatch|wublast|abblast|ncbi|decypher],default=wublast.
	-pvalue <str>           set the pvalue ,default=1e-4.
	-lib <str>              set the lib file of RepeatMasker,default RepeatMaskerLib.
	-cpu <int>              set the cpu number to use in parallel, default=3.
	-queue <str>            set the queue to use,default = <all.q>
	-outdir <str>           set the result directory , default="$Genome_basename.repeat.result".
	-help                   output help information to screen.


=head1 #########################################THE END###########################################################

=cut



use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);


my ($Genome,$Cpu,$Outdir,$Help,$Species,$Split,$Trf,$RepeatMasker,$RepeatProteinMask,$Lib,$Pvalue,$Engine,$queue);

GetOptions(
	"genome=s"=>\$Genome,
	"cpu:s"=>\$Cpu,
	"trf"=>\$Trf,
	"repeatMasker"=>\$RepeatMasker,
	"repeatProteinMask"=>\$RepeatProteinMask,
	"engine:s"=>\$Engine,
	"pvalue:s"=>\$Pvalue,
    "species:s"=>\$Species,
	"queue:s" =>\$queue,
#	"split"=>\$Split,
	"lib:s"=>\$Lib,
	"outdir:s"=>\$Outdir,
	"help"=>\$Help
);

$Cpu ||=3;
$Pvalue ||=1e-4;
$Engine ||='wublast';
$queue ||='all.q';
$Split = 1;

die "You can choose only one species option (including -species and -lib) at a time." if ((defined $Species) && (defined $Lib));

$Species ||= 'all';
$Lib ||='default';
$Lib = rel2abs $Lib;
#$Outdir ||= './Mask_result';

die `pod2text $0` if ($Help || ((not defined $RepeatMasker ) && (not defined $RepeatProteinMask) && (not defined $Trf)) || (not defined $Genome));
#die "The default output file:./Mask_result has existed,please define -outdir.\n" if (-e './Mask_result');

my $Genomefile=&absolute_path($Genome);
my $Genome_basename=basename $Genome;
$Genome_basename =~s/\.[^.]*$//;

if (length($Genome_basename) >= 30){
        print STDERR "Warning:The input file name is so long,use the top 25bp name.\n";
        $Genome_basename = substr($Genome_basename,0,25);
}

my $work_file=$Genome_basename . '.repeat.work';
my $processRepeats = '/pipeline/genomeGPON/RepeatMasker/ProcessRepeats';

$Outdir ||= "./$Genome_basename.repeat.result";
die "The default output file: $Outdir has existed,please define -outdir.\n" if (-e "$Outdir");

mkdir $work_file or die "cannot make $work_file directory:$!";
chdir $work_file;

###############SPLIT############
my (@split_path,$split_file);
$split_file = 'Split_fasta';
mkdir $split_file or die "cannot make $split_file directory:$!";
if ($Split){
	@split_path = &split_fasta($Genomefile,$Genome_basename,$split_file);
}else{
	push $Genomefile,@split_path;
}
#@split_path = &split_fasta($Genomefile,$Genome_basename,$split_file) if ($Split);

############RUN MASK#############
my $run_file="$Genome_basename.repeat";
my $Trf_file="$run_file/Trf";
my $RepeatMasker_file="$run_file/RepeatMasker";
my $RepeatProteinMask_file="$run_file/RepeatProteinMask";
mkdir $run_file || die "cannot make $run_file directory:$!";
mkdir $Trf_file || die "cannot make $Trf_file directory:$!" if ($Trf);
mkdir $RepeatMasker_file || die "cannot make $RepeatMasker_file directory:$!" if($RepeatMasker);
mkdir $RepeatProteinMask_file || die "cannot make $RepeatProteinMask_file directory:$!" if($RepeatProteinMask);



my $cmd;
my $work_shell_file="Run_mask.shell";
#open MASK,">Run_mask.shell" || die "Cannot creat file Run_mask.shell:$!";

for (my $i=0;$i<=$#split_path;$i++){
	my $split_path_basename=basename $split_path[$i];
	my $input_file=&absolute_path($split_path[$i]);
#	my $input_file_basename=basename $input_file;
	if($Trf){
		my $output_file="$Trf_file/$split_path_basename";
		mkdir "$output_file" || die "cannot make $output_file directory$!";
		print "$output_file\n$input_file\n$split_path_basename\n";
		$cmd .="cd $output_file;ln -s $input_file $split_path_basename;".
			"trf $split_path_basename 2 7 7 80 10 50 2000 -h -d\n";
	}

	if($RepeatMasker){
		my $output_file="$RepeatMasker_file/$split_path_basename";
		mkdir "$output_file" || die "cannot make $output_file directory$!";
		$cmd .="cd $output_file;ln -s $input_file $split_path_basename;";
        $cmd .="RepeatMasker -nolow -s -no_is -gff -norna -parallel 1 -engine $Engine";
        $cmd .=" -species $Species" if ($Lib =~ m/default/);
        $cmd .=" -lib $Lib " unless ($Lib =~ m/default/);
#		$cmd .=" -dir . $split_path_basename\n";
		$cmd .=" $split_path_basename\n";
	}
	if($RepeatProteinMask){
		my $output_file="$RepeatProteinMask_file/$split_path_basename";
		mkdir "$output_file" || die "cannot make $output_file directory:$!";
		$cmd .="cd $output_file;ln -s $input_file $split_path_basename;RepeatProteinMask ".
			"-noLowSimple -pvalue $Pvalue $split_path_basename\n";
	}
}
open OUT,">Run_mask.shell" || die "Cannot creat file Run_mask.shell:$!";
print OUT $cmd;
close OUT;

#system "multi_cpu.pl $Cpu $work_shell_file";
system "qsubSge -q $queue -m $Cpu $work_shell_file";


################Convert GFF format##############
chdir "..";
mkdir "$Outdir" || "cannot make $Outdir directory:$!";
my $Trf_result="$Outdir/$Genome_basename.trf.result";
my $Trf_result_gff="$Outdir/$Genome_basename.trf.gff";
my $RepeatMasker_result="$Outdir/$Genome_basename.RepeatMasker.result";
my $RepeatMasker_result_gff="$Outdir/$Genome_basename.RepeatMasker.gff";
my $RepeatMasker_result_sum="$Outdir/$Genome_basename.RepeatMasker.tbl";
my $RepeatMasker_result_mask = "$Outdir/$Genome_basename.RepeatMasker.masked";
my $RepeatProteinMask_result="$Outdir/$Genome_basename.RepeatProteinMask.result";
my $RepeatProteinMask_result_gff="$Outdir/$Genome_basename.RepeatProteinMask.gff";

`for i in $work_file/$Trf_file/*;do for j in \$i/*.dat;do cat \$j >> $Trf_result;done;done;` if ($Trf);

#`for i in $work_file/$RepeatMasker_file/*;do for j in \$i/*.fa.out;do cat \$j >> $RepeatMasker_result;done;done;` if($RepeatMasker);
&ProcessRepeats($RepeatMasker_result,$RepeatMasker_result_sum,$RepeatMasker_result_mask,$RepeatMasker_file,\@split_path) if ($RepeatMasker);

`for i in $work_file/$RepeatProteinMask_file/*;do for j in \$i/*.annot;do cat \$j >> $RepeatProteinMask_result;done;done;` if ($RepeatProteinMask);


&Trf_to_gff($Trf_result,$Trf_result_gff) if ($Trf);
&RepeatMasker_to_gff($RepeatMasker_result,$RepeatMasker_result_gff) if ($RepeatMasker);
&RepeatProteinMask_to_gff($RepeatProteinMask_result,$RepeatProteinMask_result_gff) if ($RepeatProteinMask);


print  STDERR "Congratulations,all tasks have finished!\n";


#################SUB###############################

sub ProcessRepeats {
    my ($outFile1,$outFile2,$outFile3,$inFile,$inName) = @_;
    my ($allFile,$query,$linkGenomeFile);
    $inFile = $work_file . '/' . $inFile;
    $linkGenomeFile = $inFile . '/' . $Genome_basename. '.fa';
    $allFile = $linkGenomeFile . '.RepeatMasker.cat';



    for my $file (@{$inName}){
        my $baseName = basename $file;
        my $tblFile = "$inFile/$baseName/$baseName.tbl";
        my $catFile = "$inFile/$baseName/$baseName.cat";
       
        if ((not defined $query) && (-e $tblFile)){
        	open IN,"$tblFile" || die "canot open $tblFile!:$!";
        	while (<IN>){
        		next unless (/species.*?(\S+)\s*$/);
        		$query = $1;
        	}
        	close IN;
        }
            
        if (-e $catFile){
        	`cat $catFile >> $allFile`;        
        }elsif(-e "$catFile.gz"){
        	`gzip -cd $catFile.gz >> $allFile`;
        }else{
        	 print STDERR "\nWarning!!!:Maybe the subSplitFastaFile $baseName does not run RepeatMasker successfully,
             you need to run it again by finding the $baseName line from Run_mask.shell.\n";
        }
    }
    `ln -s $Genomefile $linkGenomeFile`;
    my $cmd = "$processRepeats -gff -species \"$query\" -orifile $linkGenomeFile -maskSource $linkGenomeFile $allFile";
    open OUT,">>$work_file/$work_shell_file" || die "cannot open $work_file/$work_shell_file:$!";
    print OUT "$cmd";
    close OUT;

    system  ($cmd);

    `cp $linkGenomeFile.RepeatMasker.out $outFile1`;
    `cp $linkGenomeFile.RepeatMasker.tbl $outFile2`;
    `cp $linkGenomeFile.masked $outFile3`;

}

            
            
            

sub Trf_to_gff {
	my ($input,$output)=@_;
	my $trf_count = 0;
	my $name_root;
	open IN,$input || die "Cannot open file $input:$!";
	open OUT,">$output" || die "Cannot creat file $output:$!";
	while (<IN>) {
		if(/^Sequence: (\S+)/){$name_root=$1};
		my @trf_parts = split;
		my $num_trf_parts = @trf_parts;

		if ($num_trf_parts == 15) {
			$trf_count++;
			my $trf_count_pad = $trf_count;

			my $start = $trf_parts[0];
			my $end = $trf_parts[1];
			($start,$end)=($start<$end)?($start,$end):($end,$start);
			my $repeat_word = $trf_parts[13];
			my $wordlen = length($repeat_word);
			my $len_alias;
			if ($wordlen == 1) {
			 	$len_alias = "monomer";
			}elsif($wordlen == 2) {
				$len_alias = "dimer";
			}elsif($wordlen == 3) {
				$len_alias = "trimer";
			}elsif($wordlen == 4) {
				$len_alias = "tetramer";
			}elsif($wordlen == 5) {
				$len_alias = "pentamer";
			}elsif($wordlen == 6) {
				$len_alias = "hexamer";
			}elsif($wordlen == 7) {
				$len_alias = "heptamer";
			}elsif($wordlen == 8) {
				$len_alias = "octamer";
			}elsif ($wordlen == 9) {
				$len_alias = "nonamer";
			}elsif ($wordlen == 10) {
				$len_alias = "decamer";
			}else{
				$len_alias = $wordlen."mer";
			}

			my $attribute = "TRF$trf_count_pad";
			$attribute = "ID=TRF$trf_count_pad".
				"; Alias=".$len_alias.
				"; Name=$repeat_word";
			my $gff_str = "$name_root\t".
				"TRFv4.04\t".
				"tandem_repeat\t".
				"$start\t".
				"$end\t".
				".\t".
				".\t".
				".\t".
				"$attribute\n";
			print OUT "$gff_str";
		}
	}
	close IN;
	close OUT;
}	






sub RepeatMasker_to_gff {
	my ($input,$output)=@_;
	open IN,$input || die "Cannot open file $input:$!";
	open OUT,">$output" || die "Cannot creat file $output:$!";
	while (<IN>){
		s/^\s+//g;
		next if (/^#/ || /^SW/ || /^score/ || /^\s*$/ || /^There/);
		my @split=split /\s+/;
		my ($start,$end,$strand);
		($start,$end)=($split[5]<$split[6])?($split[5],$split[6]):($split[6],$split[5]);
		$strand=($split[8] eq '+') ? '+' : '-';
        ($split[11],$split[12]) = ($strand eq '+') ? ($split[11],$split[12]) : ($split[12],$split[13]);
		print OUT "$split[4]\tRepeatMasker\t$split[10]\t$start\t$end\t$split[1]\t$strand\t.\tTarget \"Motif:$split[9]\" $split[11] $split[12]\n";
#		print OUT $_;
	}
	close IN;
	close OUT;
}

sub RepeatProteinMask_to_gff {
	my ($input,$output)=@_;
	open IN,$input || die "Cannot open file $input:$!";
	open OUT,">$output" || die "Cannot creat file $output:$!";
	while (<IN>){
		next if (/^(#|pValue)/);
		chomp;
		my @split=split;
        if (@split == 11){
                print OUT "$split[3]\tRepeatProteinMask\t$split[8]\t$split[4]\t$split[5]\t$split[1]";
        		print OUT "\t$split[6]\t.\tTarget $split[7] $split[9] $split[10]\n";
        }
        else{
                print OUT "$split[3]\tRepeatProteinMask\tNon\t$split[4]\t$split[5]\t$split[1]";
                print OUT "\t$split[6]\t.\tTarget $split[7] $split[8] $split[9]\n";
        }
	}
	close IN;
	close OUT;
}


	
sub split_fasta {
	my ($input_file,$output_name,$output_file)=@_;
	my @id;
	my $subdir = "000";
	my $loop = 0;
	my $total_length=0;
#open IN,$input_file || die "Cannot open $input_file :$!\n";
#	if($loop % 200 == 0){
#		$subdir++;
#		mkdir("$output_file/$subdir");
#	}
	open IN,$input_file || die "Cannot open $input_file :$!\n";
	$/='>';<IN>;$/="\n";
	while (<IN>){
		my ($id,$seq,$length);
		if (/^(\S+)/){
			$id=$1;
		}else{
			die "No access number found in header line of fasta file:$input_file!\n";
		}
		$/=">";
		$seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$length=length $seq;
		$total_length+=$length;
		if ($total_length ==$length || $total_length>=10000000){
			if ($total_length>$length){close OUT;$total_length=1;}
			if($loop % 200 == 0){$subdir++;mkdir("$output_file/$subdir");}
			my $output_file_path="$output_file/$subdir/$output_name.split.$loop.fa";
			open OUT,">$output_file_path" || die "Cannot open $output_file_path :$!\n";
			push @id,$output_file_path;
			$loop++;
		}
		$/="\n";
		print OUT "\>$id\n$seq\n";
	}
	close IN;
	close OUT;
	return @id;
}


sub absolute_path{
	my $file=shift;
	my $current_path=`pwd`;
	chomp $current_path;
	my $shell_absolute_path;
	if ($file=~/^\//) {
		$shell_absolute_path = $file;
	}else{
		$shell_absolute_path = $current_path."/".$file;
	}
	return $shell_absolute_path;
}

