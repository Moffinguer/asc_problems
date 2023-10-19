#/usr/bin/perl

package Window;

use strict;
use warnings;

use Data::Dumper;

sub euclidean_distance    #( $x_00, $x_10, $x_01, $x_11 ) -> (dist)
{
	my ( $x_00, $x_01, $x_10, $x_11 ) = @_;
	my $dist = sqrt ( ( $x_00 - $x_01 )**2 + ( $x_10 - $x_11 )**2 );
	return $dist;
}

sub window                #(population, vecinity) -> (subproblems_vecinity)
{
	my $population    = shift;
	my $vecinity      = shift;
	my $step          = 1 / ( $population - 1 );
	my @lambda_window = ();
	my $lambd         = 0;
	for ( 0 .. $population - 1 )
	{
		## Add a list of N values, by each subproblem + 1 with the index of the individual
		push @lambda_window, [ $lambd, 1 - $lambd, $_ ];
		$lambd += $step;
	}

	my %subproblems_vecinity;

	my $number_near = int ( $population * $vecinity );
	for my $index ( 0 .. $population - 1 )
	{
		my @new_window = sort {
			euclidean_distance ( $lambda_window[$index][0], $a->[0],
								 $lambda_window[$index][1], $a->[1]
			) <=> euclidean_distance ( $lambda_window[$index][0],
										   $b->[0], $lambda_window[$index][1],
										   $b->[1] );
		} @lambda_window;

		@new_window = splice ( @new_window, 0, $number_near );
		$subproblems_vecinity{$index} = \@new_window;
	}

	#foreach my $clave( keys %subproblems_vecinity )
	#{
	#print("Individuo $clave =>\n". Dumper($subproblems_vecinity{$clave})."\n");
	#}
	return %subproblems_vecinity;
}

1;
