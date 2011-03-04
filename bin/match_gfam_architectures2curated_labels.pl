#!/usr/bin/perl -w
#match_gfam_architectures2curated_labels.pl
#The Arabidopsis Information Resource (TAIR)
#author: Chris Wilks
#date: 1/24/2011
#description: 	
#		associates the current run of genes and their domain 
#		architectures to curated labels which are already assigned to 
#		existing domain architectures, using the gene name and domain 
#		architecture as the association key

#NOTE: this script requires that the common Linux utility 'sort' be installed and accessible through the current path.

use strict;

#used to get the Locus part of the gene ID
my $GENE_PATTERN='^([^\.]+)\.\d+$';
#used to as part of a substitution to "NULL" out the current NOVEL##### domain ID
my $NOVEL_DOMAIN_PATTERN='NOVEL(\d+)';
#used in the same substitution as the "NULL" NOVEL##### domain ID
my $NULL_NOVEL_DOMAIN='NOVEL00000';

#used to delineate when there is an architecture which exists in the previous release but does not have a curated label
#NOTE: "desc" refers to functional description/label and "arc" refers to the domain architecture
#throughout the script
my $UNDEFINED_DESC=1;

#mapping file between genes and the previous set of domain architectures, also includes the curated labels for those domain
#architectures which have them.
#format: Gene Name<TAB>Domain Architecture<TAB>New Label
my $previous_descs_file = shift; 
#output from GFAM for newest set of proteins (domain_architectures.tab)
#format: Gene Name<TAB>Sequence Length<TAB>Domain Architecture (rest of the columns unused)
my $domain_arcs_file = shift;

if(!$previous_descs_file || !$domain_arcs_file)
{
	print_usage();
	exit(-1);
}
	
#we output three files (exclusive: each gene is output to ONE of the three files):
#1) those genes with domain architectures but with no matching label in the previous label set
open(NONE,">$domain_arcs_file.not_matching");
#2) those genes with domain architectures and with a matching label in the previous label set
open(LABEL,">$domain_arcs_file.label");
#3) those genes without a domain architecture
open(NO_LABEL,">$domain_arcs_file.no_label");
	
main();

close(NONE);
close(LABEL);
close(NO_LABEL);

sub print_usage
{
	print "Usage: match_gfam_architectures2curated_labels.pl previous_descriptions_file new_domain_architectures_file\n";
}

sub main
{
	my $stat = system("sort -s -k2,2 $previous_descs_file > $previous_descs_file.tmp_sorted");
	if($stat != 0)
	{
		print "Could not sort $previous_descs_file, exiting\n";
		exit(-1);
	}
	my ($locus2func_href,$arc2desc_href) = load_previous_descriptions("$previous_descs_file.tmp_sorted");
	match_genes2descriptions($locus2func_href,$arc2desc_href,$domain_arcs_file);
	`rm -rf $previous_descs_file.tmp_sorted`;
}

#main category determination routine, 
#splits the input into the three exclusive output files
sub match_genes2descriptions
{

	my ($locus2func_href,$arc2desc_href,$domain_arcs_file) = @_;
	
	open(IN,"<$domain_arcs_file");
	
	my %seen;
	while(my $line = <IN>)
	{
		chomp($line);
		my ($gene,$size,$arc) = split(/\t/,$line);
	
		my $desc = $arc2desc_href->{$arc};
		
		if($arc eq "NO_ASSIGNMENT")	
		{
			print NONE "$gene\t$arc\n";
		}
		elsif($arc =~ /NOVEL/i)
		{
			process_novel($line,$locus2func_href);
		}
		elsif(!$desc)
		{
			print NO_LABEL "KNOWN_UPDATED_ARCHITECTURE\t$gene\t$arc\n";
		}
		else
		{
			if($desc ne $UNDEFINED_DESC)
			{
				print LABEL "KNOWN_WITH_LABEL\t$gene\t$arc\t$desc\n";
			}
			else
			{
				print NO_LABEL "KNOWN_NON_CURATED_LABEL\t$gene\t$arc\n";
			}
		}	
	}
	close(IN);
}

#define additional subcategories for those domain architectures with NOVEL domains in them
sub process_novel
{
	my ($line,$locus2func_href) = @_;
	
	my ($gene,$size,$arc) = split(/\t/,$line);
	$gene =~ s/$GENE_PATTERN/$1/;
	my $locus = $1;
	my $newarc = $arc;
	
	my $info_aref = $locus2func_href->{$locus};	
	if(defined($info_aref))
	{
		my $oldarc = $info_aref->[0];
		$newarc=~s/$NOVEL_DOMAIN_PATTERN/$NULL_NOVEL_DOMAIN/g;

		if($newarc eq $oldarc)
		{
			print LABEL "NOVEL_WITH_LABEL\t$gene\t$arc\t".$info_aref->[1]."\n" if(defined($info_aref->[1]));
			print NO_LABEL "NOVEL_NON_CURATED_LABEL\t$gene\t$arc\n" if(!defined($info_aref->[1]));
		}
		else
		{
			print NO_LABEL "NOVEL_UPDATED_ARCHITECTURE\t$gene\t$arc\n";
		}
	}
	else
	{
		print NO_LABEL "NOVEL_NEW_LOCUS\t$gene\t$arc\n";
	}
}

		
sub load_previous_descriptions
{
	my $file = shift;
	my %locus2func;
	my %arc2desc;
	open(IN,"<$file");
	while(my $line = <IN>)
	{
		chomp($line);
		my ($gene,$arc,$desc) = split(/\t/,$line);
	
		#replace gene name with locus name
		$gene =~ s/$GENE_PATTERN/$1/;
		#replace whatever novel domain ID ##### we have with the NULL version
		#this is done since the old novel domain ID #'s are obsolete at this point and DO NOT match the 
		#new set of domain ID #'s
		$arc=~s/$NOVEL_DOMAIN_PATTERN/$NULL_NOVEL_DOMAIN/g;
		
		#even if there is no label defined for this architecture, we still need to 
		#record that there was an architecture last time, otherwise just record the label/description itself
		my $desc_status = (length($desc) < 1?$UNDEFINED_DESC:$desc);
		$arc2desc{$arc}=$desc_status;
		
		$desc_status = (length($desc) < 1?undef:$desc);
		$locus2func{$gene}=[$arc,$desc_status];
	}
	close(IN);
	return (\%locus2func,\%arc2desc);
}

#SDG
