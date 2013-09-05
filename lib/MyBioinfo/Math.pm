package MyBioinfo::Math;

use 5.006;
use strict;
use warnings;
use Math::CDF qw(pchisq);
use Data::Dumper;
use MyBioinfo::Common qw(sum);
require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use MyBioinfo::Common ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	gtest_gof chisqtest_gof
) ] );

our @EXPORT_OK = @{ $EXPORT_TAGS{'all'} };

our @EXPORT = ();

our $VERSION = '0.11';

############### CHANGES: ###################
# Mar. 1 2012  15:25  Added protection in G-test when round-off error causes
#                     G-statistic to be a small negative.
#

######## Preloaded methods go here. ##############

# Use G-test to perform Goodness-of-fit test.
sub gtest_gof{
	my $ra = shift;	# reference to array of numbers.
	my $rp;	# reference to array of probabilities.
	# If vector p is not specified, assume uniform probability.
	if(@_ > 0){
		$rp = shift;
		if(!&p_sanity($rp)){
			die "Probability vector contains elements that are zero or negative. G-test was not performed. Abort.\n";
		}
	}else{
		$rp = [ (1/@{$ra}) x @{$ra} ];
	}
	# Fill out the vector of expectated values.
	my @a_exp = fill_exp($ra, $rp);
	# Calcualte G-statistic.
	my $g_stat = 0;
	for(my $i = 0; $i < @a_exp; $i++){
		if($ra->[$i] != 0){	# skip zero entry in observation vector.
			$g_stat += $ra->[$i] * log( $ra->[$i] / $a_exp[$i] );
		}
	}
	$g_stat *= 2.0;
	# Use Chi-square distribution to calculate p-value.
	my $df = @{$ra} - 1;
	my $p;	# p-value.
	if($g_stat <= 0){	# this may happen due to round-off errors.
		$p = 1;
	}else{
		$p = 1 - pchisq($g_stat, $df);
	}
	return ($g_stat, $df, $p);
}

# Use Chi-square test to perform Goodness-of-fit test.
sub chisqtest_gof{
	my $ra = shift;	# reference to array of numbers.
	my $rp;	# reference to array of probabilities.
	# If vector p is not specified, assume uniform probability.
	if(@_ > 0){
		$rp = shift;
		if(!&p_sanity($rp)){
			die "Probability vector contains elements that are zero or negative. Chi-square test was not performed. Abort.\n";
		}
	}else{
		$rp = [ (1/@{$ra}) x @{$ra} ];
	}
	# Fill out the vector of expectated values.
	my @a_exp = fill_exp($ra, $rp);
	# Calcualte Chisq-statistic.
	my $ch_stat = 0;
	for(my $i = 0; $i < @a_exp; $i++){
		$ch_stat += ( $ra->[$i] - $a_exp[$i] ) ** 2 / $a_exp[$i];
	}
	# Use Chi-square distribution to calculate p-value.
	my $df = @{$ra} - 1;
	my $p = 1 - pchisq($ch_stat, $df);
	return ($ch_stat, $df, $p);
}

# Fill out a vector of expected values. A common subroutine used by both G-test 
# and Chi-square test.
sub fill_exp{
	my($ra, $rp) = @_;
	my $v_sum = sum(@{$ra});
	my @a_exp;
	for(my $i = 0; $i < @{$ra}; $i++){
		push @a_exp, $v_sum * $rp->[$i];
	}
	return @a_exp;
}

# Check the sanity for probability vector: all elements should be larger than zero.
sub p_sanity{
	my $p_ref = shift;
	foreach my $p(@{$p_ref}){
		return 0 if $p <= 0;
	}
	return 1;
}



1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyBioinfo::Math - Math routines used for bioinformatics.

=head1 SYNOPSIS

  use MyBioinfo::Common qw( gtest_gof chisqtest_gof );

  gtest_gof
  chisqtest_gof


=head1 DESCRIPTION

  Just some small functions that I found to be useful for my bioinformatics
  work but could not find in Perl release. I put them here whenever I thought
  about sth. It's just a start...

=head2 EXPORT

  None.

=head1 SEE ALSO

  MyBioinfo::Common

=head1 AUTHOR

Li Shen, E<lt>li.shen@mssm.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Li Shen

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.6.0 or,
at your option, any later version of Perl 5 you may have available.


=cut
