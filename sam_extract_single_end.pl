    Bayesian epiallele detection
    Copyright (C) 2019 James E. Barrett (regmjeb@ucl.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

# sam file extraction
# based on /Users/James/ UCL/TRACERx/ENCODE/R_June15/data_extraction.R

use List::Util qw[min max];

open(my $in, "<", "test3.sam");
open(my $chr1, ">", "chr1.txt");
open(my $chr2, ">", "chr2.txt");

while(<$in>){

	my @sam = split /\t/;

	my $chr = $sam[2];
 	my $pos = $sam[3];
	my $CIGAR = $sam[5];
 	my $read_length = length($sam[13]) - 5;
 	my $m1 = substr($sam[13], 5, $read_length);
 	# call the absolute.reference function to obtain the coordinates of each nucleotide
	my @ar1 = abs_ref($CIGAR, $pos, $read_length);
	
	# Pull out the indices where there is a 'z' or 'Z'
	my @z_index1 = ();
	my @Z_index1 = ();
	while ($m1 =~ /z/g) {push @z_index1, pos $m1}
	while ($m1 =~ /Z/g) {push @Z_index1, pos $m1}

	# Pull out the corresponding genome coordinates of the 'z's and 'Z's
	my @z_loc1 = ();
	my @Z_loc1 = ();
	
	push @z_loc1, @ar1[$_-1] foreach @z_index1;
	push @Z_loc1, @ar1[$_-1] foreach @Z_index1;
	
	# Take the set union of these coordinates (this combines the paired reads into one)
	my @z_loc = @z_loc1;
	my @Z_loc = @Z_loc1;
	
	# Save to file
	if ($chr == 1) {$out = $chr1}
	elsif ($chr == 2) {$out = $chr2};
	print $out join("\t",
		$sam[0],						# name
		$sam[3],						# start
		max(max(@z_loc),max(@Z_loc)),	# end
		join(":",@z_loc),				# unmeth coords
		join(":",@Z_loc)), "\n";		# meth coords
}

close $in;
close $out;


#=============================================================================#
# function: abs_ref
#=============================================================================#	

sub abs_ref{

	my $CIGAR = $_[0];
	my $pos = $_[1];
	my $read_length = $_[2];	
 	my $n_cigar = length($CIGAR);
 	
 	# An empty vector of length READ.LENGTH
	my @absolute_reference = (0) x $read_length;

 	# Pull out any M, I, or D from CIGAR string
  	my @matches = $CIGAR =~ /[MID]/g;
	
	# Pull out corresponding locations in CIGAR string
	my @location = ();
	while ($CIGAR =~ /[MID]/g){
		push @location, pos $CIGAR;
	};

	# d = total number of letters in CIGAR
	my $d = $#location;

	# An index to keep track of our current position on the read
	my $read_position = 0;
   
	# An index to keep track of our current position on the reference genome
	my $reference_position = 0;
	
	# Begin loop over the d letters
	for ($mu=0; $mu <= $d; $mu++){

      # Here we want to extract the number corresponding to the mu-th letter
      # The start position of the number in the CIGAR string
		if ($mu > 0) {$start = @location[$mu-1]} else {$start = 0}

		# The end position of the number in the CIGAR string
		my $stop = @location[$mu]-1;
		
		# len = the number itself corresponding to the d-th letter
		my $len = substr($CIGAR, $start, $stop-$start);
		
		# If mu-th letter is 'M' continue to add indices
		if (@matches[$mu] eq 'M'){
			@absolute_reference[$read_position..($read_position+$len-1)] = ($reference_position..($reference_position+$len-1));
			$read_position = $read_position + $len;
			$reference_position = $reference_position + $len;
		}
		# If mu-th letter is 'I' add NA and advance read.position
		if (@matches[$mu] eq 'I'){
			@absolute_reference[$read_position..($read_position+$len-1)] = ('NA') x $len;
			$read_position = $read_position + $len;
		}

		# If mu-th letter is 'D' advance reference.position only
		if (@matches[$mu] eq 'D'){
			$reference_position = $reference_position + $len;
		}

	}
	
	@absolute_reference[$_] = @absolute_reference[$_] + $pos foreach 0..($read_length-1);	
	 
	return @absolute_reference
}


