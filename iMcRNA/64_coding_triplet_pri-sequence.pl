#!/usr/bin/perl -w

use strict;

my $input_file_1 = $ARGV[0];
open(IN_1, "$input_file_1")
	or die "can't open the input file : $!";
my $output_file = $ARGV[1];
open OUT_1, ">$output_file"
	or die "Can not open $output_file : $!";

my $num_total=0;
my $title; my $seq_letter;

while(<IN_1>){
	if($_ =~ />()/){
		$num_total++;
		$title=$_;
		#print $title;
		#print OUT_1 $title;
	}
	elsif($_ =~ /([ACUG]+)/){
		#print "$_";
		$seq_letter = $1;
		
		my @coding_table=translate_to_coding($seq_letter);			
		unite_probability(\@coding_table);
		
		#print feature
		my $len_coding_table = @coding_table;

		for(my $i=0; $i<$len_coding_table; $i++){
			print OUT_1 "$coding_table[$i] ";
		}
		print OUT_1 "\n";
	}
}

sub translate_to_coding{
	
	my ($seq_letter)=@_;
	
	my @table;
	for(my $i=0; $i<64; $i++){
		push(@table, 0);
	}
	
	#statistic %XYZ
	my $len=length($seq_letter);
	for(my $i=0;$i<$len-2;$i++){
		my $triplet_seq=substr($seq_letter,$i,3);
		#print $triplet_seq."\n";
		
		if(length($triplet_seq)!=3){
			die "$triplet_seq is not local triplet elements";
		}
		
		my $first_letter=substr($triplet_seq,0,1);
		my $second_letter=substr($triplet_seq,1,1);
		my $third_letter=substr($triplet_seq,2,1);
		#calculate the combin value
		my $first_letter_value=get_letter_value($first_letter);
		my $second_letter_value=get_letter_value($second_letter);
		my $third_letter_value=get_letter_value($third_letter);				
		my $comb_value=$first_letter_value*16+$second_letter_value*4+$third_letter_value*1;		
		$table[$comb_value]++;		
	}	
	return @table;
}

sub get_letter_value{
	my ($char) = @_;
	my ($l_char, $ret);
	
	$l_char = lc($char);
	$ret = -1;
	if($l_char eq 'a'){
		$ret = 0;	
	}
	elsif($l_char eq 'g'){
		$ret = 1;	
	}
	elsif($l_char eq 'c'){
		$ret = 2;	
	}
	elsif($l_char eq 'u'){
		$ret = 3;	
	}
	elsif($l_char eq 't'){
		$ret = 3;	
	}
	else{
		print "ERROR: has not A G C U \n";
	}
	
	return $ret;
}

sub unite_probability{
	my ($table) = @_;
	my ($len_table, $i, $sum);
#	local($show_table);
	
	$len_table = @$table;
	$sum = 0;
	for($i=0; $i<$len_table; $i++){
		$sum = $sum + $$table[$i];			
	}
#	print OUT_1 "sum is $sum\n";
	for($i=0; $i<$len_table; $i++){
		$$table[$i] = $$table[$i] / $sum;	
	}
}