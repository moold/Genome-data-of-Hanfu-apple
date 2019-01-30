#!/usr/bin/perl -w
use strict;

=head1 ########################################USAGE#####################################################

=head1 Name

    Run_Repeatscout_Piler_Ltrfinder.pl --The pipeline for RepeatSquence annotation.

=head1 DESCRIPTION

    This script include Repeatscout , Piler , Ltrfinder , RepeatModeler and than output each software result with fasta forma    t,also output the merge result to be used the denovo lib for RepeatMask.

    Please Note:
    1.For large Genome ,you can not run Piler,because piler has RAM BUG for large genome.
    2.The RepeatModeler includes Repeatscout and ReCon ,so when you run RepeatModeler ,you can not run Repeatscout.
    3.For >500M genome ,only you need to run RepeatModeler && Ltrfinder.
    4.You can specially define -funtion to filter the denovo lib ,which can filter the repeat squences with funtion basing
      on -funlib database , but it need enormous time, the total time of this step Trembl >> Kegg >> Swissprot,recommend Kegg

    You can defined Ltrfinder lib by refering to [/home/huj/program/gene_predict/LTR_FINDER.x86_64-1.0.5/tRNA.list.txt].

=head1 Version

    Author: Jiang Hu.
    Version: 1.0.

=head1 Usage

    perl Run_Repeatscout_Piler_Ltrfinder.pl -genome genome_path -repeatscout -piler -ltrfinder -function -cpu 10 
    -ltrlib Athal-tRNAs.fa
    
    -genome <str>           set the genome path file,absolute paths,relative paths are allowed,needed.
    -repeatscout            set run RepeatScout.
    -piler                  set run Piler.
    -repeatmodeler          set run RepeatModeler.
    -ltrfinder              set run Ltrfinder.
    -function               set filter basing on detabase,recommend.
    -funlib <str>           set the library to filter repeat regions with gene,accept [Kegg<default>|Swissprot|Trembl].
    -ltrlib <str>           set the Ltrfinder lib ,default Athal-tRNAs.fa.
    -cpu <int>              set the cpu number to use in parallel, default=5.
    -workdir <str>          set the work directory,default $genome_basename.Run.Repeatscout_Piler_Ltrfinder.
    -outdir <str>           set the output directory,default $genome_basename.Repeat.lip.Result.
    -help                   output help information to screen.



=head1 #########################################THE END###########################################################

=cut




use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use Cwd;

my ($Genome,$Cpu,$Outdir,$Help,$Repeatscout,$Function,$Workdir,$Piler,$Ltrfinder,$Ltrlib,$RepeatModeler,$Functionlib);
my ($genome_basename,$genome_abso_path,$libpath,$Ps_scan,$Cwd,$RepeatModelerPath);

GetOptions(
	"genome=s"=>\$Genome,
	"repeatscout"=>\$Repeatscout,
	"piler"=>\$Piler,
	"ltrfinder"=>\$Ltrfinder,
    "repeatmodeler"=>\$RepeatModeler,
	"function"=>\$Function,
    "funlib:s"=>\$Functionlib,
	"cpu:s"=>\$Cpu,
	"ltrlib:s"=>\$Ltrlib,
	"help"=>\$Help,
	"outdir:s"=>\$Outdir,
	"workdir:s"=>\$Workdir
);

die `pod2text $0` if ($Help || (not defined $Genome) || ((not defined $Repeatscout) && (not defined $Piler) && (not defined $Ltrfinder) && (not defined $RepeatModeler)));

$genome_basename = basename ($Genome);
$genome_abso_path = rel2abs ($Genome);
$Cwd = cwd();
$Cpu ||=5;
#$path = '/home/huj/program/gene_predict/LTR_FINDER.x86_64-1.0.5/tRNAdb/';
$Ltrlib ||= 'Athal-tRNAs.fa';
$Ps_scan = "/home/huj/program/gene_predict/LTR_FINDER.x86_64-1.0.5";
$libpath = "$Ps_scan/tRNAdb/$Ltrlib";
$Function = 1 if (defined $Functionlib);
$Functionlib ||= 'Kegg';
$RepeatModelerPath = '/home/huj/program/gene_predict/RepeatModeler/';

die "Not exists the lib $Ltrlib,please check!!!\n" unless (-e "$libpath");

$Workdir ||= "$genome_basename.Run.RepeatModeler_scout_Piler_Ltrfinder";
die "The default workdir $Workdir has existed,please define -workdir.\n" if (-e "$Workdir");
$Outdir ||= "$genome_basename.Repeat.lip.Result";
die "The default output file: $Outdir has existed,please define -outdir.\n" if (-e "$Outdir");
mkdir $Workdir || die "cannot make $Workdir directory:$!";
mkdir $Outdir || die "cannot make $Outdir directory:$!";
chdir $Workdir;


#my ($Cmd,$Shell_file);
#my $Shell_file="$genome_basename.repeatlib.shell";
#open OUT,">$Shell_file" || die "Cannot creat file $Shell_file:$!.\n";
############read fasta####################
my %fasta;
&read_fasta($genome_abso_path,\%fasta);

############split genome##################
my $split_file="$genome_basename.split";
mkdir $split_file || "cannot make $split_file directory:$!";
my @split_file=&split_fasta($genome_abso_path,$genome_basename,$split_file);

###############Repeatscout##################
if($Repeatscout){

	my ($Cmd,$Shell_file,$Repeatscout_file,$Repeatscout_fasta);
	$Shell_file="$genome_basename.Repeatscout.shell";
	$Repeatscout_file="$genome_basename.Repeatscout";
	$Repeatscout_fasta = "../$Outdir/$genome_basename.Repeatscout.fa";
	
	mkdir $Repeatscout_file || die "cannot make $Repeatscout_file directory:$!";
	chdir $Repeatscout_file;
	$Cmd .= "build_lmer_table -sequence $genome_abso_path -freq $genome_basename.Repeatscout.fre;\n";
	$Cmd .= "RepeatScout -sequence $genome_abso_path -freq $genome_basename.Repeatscout.fre -output $genome_basename.Repeatscout;\n";
	$Cmd .= "filter-stage-1.prl $genome_basename.Repeatscout > $genome_basename.Repeatscout.stage-1.output;\n";
	$Cmd .= "RepeatMasker $genome_abso_path -lib $genome_basename.Repeatscout.stage-1.output -dir . ;\n";
	$Cmd .= "cat $genome_basename.Repeatscout.stage-1.output |filter-stage-2.prl --cat $genome_basename.out --thresh=20 > $genome_basename.Repeatscout.stage-2.output;\n";

	open OUT,">$Shell_file" || die "Cannot creat file $Shell_file:$!.\n";
	print OUT "$Cmd";
	close OUT;
	system "sh $Shell_file";
	`cp $genome_basename.Repeatscout.stage-2.output ../$Repeatscout_fasta;`;
	chdir "..";
}

if ($RepeatModeler){

        my ($Cmd,$Shell_file,$RepeatModeler_file,$RepeatModeler_fasta);
        $Shell_file="$genome_basename.RepeatModeler.shell";
        $RepeatModeler_file="$genome_basename.RepeatModeler";
        $RepeatModeler_fasta="../$Outdir/$genome_basename.RepeatModeler.fa";

        mkdir $RepeatModeler_file || die "cannot make $RepeatModeler_file directory:$!";
        chdir $RepeatModeler_file;

        $Cmd .= "$RepeatModelerPath/BuildDatabase  -name $genome_basename -engine wublast $genome_abso_path\n";
        $Cmd .= "$RepeatModelerPath/RepeatModeler -database $genome_basename -engine wublast\n";

        open OUT,">$Shell_file" || die "Cannot creat file $Shell_file:$!.\n";
        print OUT "$Cmd";
        close OUT;
        system "sh $Shell_file";
        `cp */consensi.fa.masked ../$RepeatModeler_fasta`;
        chdir "..";
}


##################Piler#####################
if($Piler) {
	
	my($Cmd,$Shell_file,$Piler_file,$Piler_split_result,$Piler_fasta);
	$Shell_file="$genome_basename.Piler.shell";
	$Piler_file="$genome_basename.Piler";
	$Piler_split_result="$Piler_file.split.result";
	$Piler_fasta = "../$Outdir/$genome_basename.piler.fa";
	

	mkdir $Piler_file || die "cannot make $Piler_file directory:$!";
	chdir $Piler_file;
	mkdir $Piler_split_result || die "cannot make $Piler_split_result directory:$!";

	my $subdir = "000";
	my $loop = 0;
	for (my $i=0;$i<=$#split_file;$i++){
		for (my $j=$i;$j<=$#split_file;$j++){
			if($loop % 500 == 0){$subdir++;mkdir("$Piler_split_result/$subdir");}
			my $target=basename ($split_file[$i]);
			my $query =basename ($split_file[$j]);
			my $outfile="$Piler_split_result/$subdir/$target.$query.pals.gff";
			$Cmd .= "pals -target $split_file[$i] -query $split_file[$j] > $outfile;\n";
			$loop++;
		}
	}
	
	open OUT, ">$Shell_file" || die "fail creat $Shell_file:$!;\n";
	print OUT $Cmd;
	close OUT;

	`multi_cpu.pl $Cpu $Shell_file;`;
	`for i in $Piler_split_result/*;do for j in \$i/*.pals.gff;do cat \$j >> $genome_basename.pals.gff;done;done;`;
	`piler2 -trs $genome_basename.pals.gff -out $genome_basename.trs.gff;`;
	mkdir "fams" || die "cannot make fams directory:$!";
	`piler2 -trs2fasta $genome_basename.trs.gff -seq $genome_abso_path -path fams;`;
        mkdir "aligned_fams" || die "cannot make aligned_fams directory:$!";
	chdir "fams";
	`for fam in *;do muscle -in \$fam -out ../aligned_fams/\$fam -maxiters 1 -diags1; done;`;
	chdir "..";
	mkdir "cons" || die "cannot make cons directory:$!";
	chdir "aligned_fams";
	`for fam in *;do piler2 -cons \$fam -out ../cons/\$fam -label \$fam;done;`;
	chdir "../cons";
	`cat * > ../$genome_basename.piler.fa;`;
	chdir "..";
	`cp $genome_basename.piler.fa ../$Piler_fasta;`;
    chdir "..";	
}


if($Ltrfinder) {

	my($Cmd,$Shell_file,$Ltrfinder_file,$Ltrfinder_fasta);
	$Shell_file="$genome_basename.Ltrfinder.shell";
	$Ltrfinder_file="$genome_basename.Ltrfinder";
	$Ltrfinder_fasta = "../$Outdir/$genome_basename.Ltrfinder.fa";

	mkdir $Ltrfinder_file || die "cannot make $Ltrfinder_file directory:$!";
	chdir $Ltrfinder_file;
	$Cmd .= "ltr_finder -w 2 -s $libpath -a $Ps_scan $genome_abso_path > $genome_basename.Ltrfinder.result;\n"; 

	open OUT,">$Shell_file" || die "Cannot creat file $Shell_file:$!.\n";
	print OUT "$Cmd";
	close OUT;
	system "sh $Shell_file";
	&tran_fasta("$genome_basename.Ltrfinder.result","../$Ltrfinder_fasta");
	chdir "..";
}



my $denovolib_cmd="Get_repeatmask_devovolib.pl -cpu $Cpu ";
$denovolib_cmd .="-repeatscout ../$Outdir/$genome_basename.Repeatscout.fa " if (defined $Repeatscout);
$denovolib_cmd .="-repeatmodeler ../$Outdir/$genome_basename.RepeatModeler.fa " if (defined $RepeatModeler);
$denovolib_cmd .="-piler ../$Outdir/$genome_basename.piler.fa " if (defined $Piler);
$denovolib_cmd .="-ltrfinder ../$Outdir/$genome_basename.Ltrfinder.fa " if (defined $Ltrfinder);
$denovolib_cmd .="-function -lib $Functionlib " if (defined $Function);
$denovolib_cmd .="-outdir ../$Outdir";
system $denovolib_cmd;
#print "done";	

print  STDERR "Congratulations,all tasks have finished!\n";	


####################sub######################
sub tran_fasta {
	my ($infile,$outfile)=@_;
    my %hash_id;
	open IN,"$infile" || die "Cannot open $infile:$!\n";
	open OUT,">$outfile" || die "Cannot make $outfile:$!\n";
	while (<IN>){
		chomp;
		next unless (/^\[/);
		my (@split,$id,$start,$end,$head,$len,$sqe);
		@split=split /\t/;
		$id=$split[1];
		$split[2]=~/(\d+)-(\d+)/;
		($start,$end)=$1<$2 ? ($1,$2):($2,$1);
		$len=$end-$start+1;
		$head=$id.'_'.$start.'_'.$end;

        
        if (exists $hash_id{$head}){
            next;
        }else{
            $hash_id{$head} = 1;
        }
            


		$sqe=substr($fasta{$id}[0],$start-1,$len);
#		if ($split[12] eq '-'){
#			$sqe=~tr/ATGCNatgcn/TACGNtacgn/;
#			my @temp=split //,$sqe;
#			$sqe=join '',(reverse @temp);
#		}
		print OUT ">$head\n$sqe\n";
	}
	close IN;
	close OUT;
}


sub read_fasta {
	my ($infile,$hash)=@_;
	open IN,$infile || die "Cannot open $infile:$!\n";
	$/='>';<IN>;$/="\n";
	while(<IN>){
		my ($id,$seq,$length);
		if (/^(\S+)/){
			$id=$1;
		}else{
			die "No access number found in header line of fasta file:$infile!\n";
		}
		if ( $id=~/\|/ ) {
			die "No '|' allowed in the access number of fasta file:$infile!\n";
		}
		$/=">";
		$seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$length=length $seq;
		$hash->{$id}=[$seq,$length];
		$/="\n";
	}
	close IN;
}



sub split_fasta {
	my ($input_file,$output_name,$output_file)=@_;
	my @id;
	my $subdir = "000";
	my $loop = 0;
	my $total_length=0;
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
			push @id,rel2abs($output_file_path);
			$loop++;
		}
		$/="\n";
		print OUT "\>$id\n$seq\n";
	}
	close OUT;
	close IN;
	return @id;
}
	
