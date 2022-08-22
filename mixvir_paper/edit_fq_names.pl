#! /usr/bin/perl
use warnings;
use strict;


my @levels = ("10x", "50x", "100x");
#my @levels = ("100x");

my $numlevels = @levels;


for my $levnum (0..$numlevels-1) {
	
	opendir FQ, "subset_fastqs/$levels[$levnum]";
	my @R1Names = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' && $_ =~ /_1.fastq/} readdir(FQ);
	#my @R1Names = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' && $_ =~ /ERR221664_1.fastq/} readdir(FQ);
	close FQ;
	
	my $numR1files = @R1Names;
	
	foreach my $num (0..$numR1files-1) {
		
		my $currentfile = $R1Names[$num];
		my $linecounter = 1;
		my $namebase = $R1Names[$num];
		$namebase =~ s/_..fastq//;
		my $numseqs = 0;
		
		print "$currentfile\n";
		
		open IN, "subset_fastqs/$levels[$levnum]/$currentfile" or die$!;
		open OUT, ">subset_fastqs/$levels[$levnum]/$namebase-$levels[$levnum]_S30_L001_R1_001.fastq" or die$!;
		
		while(<IN>) {
			if ($linecounter == 1) { #on name line
				$numseqs++;
				print OUT "\@VH00443:11:AAC7MF7M5:1:1101:30009:$numseqs 1:N:0:1\n";
				$linecounter++;
				next;
			}	
			
			elsif ($linecounter == 2) {
				print OUT "$_";
				$linecounter++;
				next;
			}
			
			elsif ($linecounter == 3) {
				print OUT "+\n";
				$linecounter++;
				next;
			}
			
			else {
				print OUT "$_";
				$linecounter = 1;
				next;
			}
			
		}
		
		close IN;
		close OUT;
		
		my $R2file = $currentfile;
		$R2file =~ s/_1.fastq/_2.fastq/g;
		
		open IN, "subset_fastqs/$levels[$levnum]/$R2file" or die$!;
		open OUT, ">subset_fastqs/$levels[$levnum]/$namebase-$levels[$levnum]_S30_L001_R2_001.fastq" or die$!;
		
		$linecounter = 1;
		$numseqs = 0;
		
		while(<IN>) {
			if ($linecounter == 1) { #on name line
				$numseqs++;
				print OUT "\@VH00443:11:AAC7MF7M5:1:1101:30009:$numseqs 2:N:0:1\n";
				$linecounter++;
				next;
			}	
			
			elsif ($linecounter == 2) {
				print OUT "$_";
				$linecounter++;
				next;
			}
			
			elsif ($linecounter == 3) {
				print OUT "+\n";
				$linecounter++;
				next;
			}
			
			else {
				print OUT "$_";
				$linecounter = 1;
				next;
			}
			
		}
		
		close IN;
		close OUT;
		
		
		
	}	

}




#replace read descriptors

#name format:
#ERR221664-10x_S30_L001_R1_001.fastq.gz
#ERR221664-10x_S30_L001_R2_001.fastq.gz

#@VH00443:11:AAC7MF7M5:1:1101:30009:1000 1:N:0:TGTGTTAGTA+TAACAGAGTA
#@VH00443:11:AAC7MF7M5:1:1101:30009:1000 2:N:0:TGTGTTAGTA+TAACAGAGTA


#read in R1 files in dir

#Substitute R1 with R2 to get array with R2 file names.

#Create array to store samp names created for R1 files

