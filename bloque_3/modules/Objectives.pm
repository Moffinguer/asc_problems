#/usr/bin/perl

package Objectives;

use strict;
use warnings;

use Math::Trig;
use Data::Dumper;

sub best_evaluations    #(f1, f2, z) -> (new_z)
{
	my ( $f1, $f2, $z ) = @_;
	my @new_z = ( $z->[0], $z->[1] );

	#print("Z1, Z2, f1, f2:\n ".Dumper($z->[0], $z->[1], $f1, $f2)."\n");

	$new_z[0] = $f1 if not defined $z->[0] or $f1 < $z->[0];
	$new_z[1] = $f2 if not defined $z->[1] or $f2 < $z->[1];

	return @new_z;
}

sub select_objective_function    #( algorithm, x ) -> ( f1, f2 )
{
	my ( $problem, $x ) = @_;

	if ( $problem eq "zdt3" )
	{
		return zdt3 (@$x);
	}
	else
	{
		return cf6 (@$x);
	}
}

sub zdt3    #(x, n) -> (f1(x), f2(x))
{
	my (@x) = @_;

	my $f_1 = $x[0];

	my $g = 0;

	my $n = scalar (@x);

	$g += $x[$_] for ( 1 .. $n - 1 );
	$g = 9 * $g / ( $n - 1 ) + 1;

	my $h = $f_1 / $g;
	$h = 1 - sqrt ($h) - $h * sin ( 10 * $f_1 * pi );

	my $f_2 = $g * $h;

	return ( $f_1, $f_2 );
}

sub cf6    #(x, n) -> (f1(x), f2(x))
{

	my ( $n, @x ) = @_;

	@x = @{ $x[0][0] };

	my ( $f_1, $f_2 ) = ( $x[0], ( 1 - $x[0] )**2 );

	for ( 1 .. $n - 1 )
	{
		if ( $_ % 2 == 1 )
		{
			$f_1 += ( $x[$_] - .8 * $x[0] * cos ( 6 * pi * $x[0] + $_ * pi / $n ) )**2;
		}
		else
		{
			$f_2 += ( $x[$_] - .8 * $x[0] * sin ( 6 * pi * $x[0] + $_ * pi / $n ) )**2;
		}
	}
	return ( $f_1, $f_2 );
}

sub init_population #( n, upper_limit, lower_limit, dimensions, algorithm ) -> (population_list)
{
	my ( $n, $upper_limit, $lower_limit, $dimensions, $algorithm ) = @_;

	my @population_list;
	if ( $algorithm eq "zdt3" )
	{
		for my $i ( 0 .. $n - 1 )
		{
			my @element;
			for my $j ( 0 .. $dimensions - 1 )
			{
				my $value = $lower_limit + rand ($upper_limit);
				push @element, $value;
			}
			push @population_list, [@element];
		}
	}
	else
	{
		for my $i ( 0 .. $n - 1 )
		{
			my @element;
			for my $j ( 0 .. $dimensions - 1 )
			{
				if ( $j == 0 )
				{
					my $value = $lower_limit + rand ($upper_limit);
					push @element, $value;
				}
				else
				{
					my $value = -2 + rand (4);
					push @element, $value;
				}
			}
			push @population_list, [@element];
		}
	}

	#print("POPULATION LIST:\n ". Dumper(@population_list)."\n");

	return @population_list;
}

sub constraints    #(x, algorithm, gen) -> ( violated )
{
	my ( $x, $algorithm, $gen ) = @_;

	my @x = @$x;
	return 0 if $algorithm eq "zdt3";

	my $n            = scalar (@x);
	my $restriction1 = $x[1] - .8 * $x[0] * sin ( 6 * pi * $x[0] + 2 * pi / $n );

	my $restriction2 = $x[3] - .8 * $x[0] * sin ( 6 * pi * $x[0] + 4 * pi / $n );

	my $sg1 = .5 * ( 1 - $x[0] ) - ( 1 - $x[0] )**2;
	my $sg2 = .25 * sqrt ( 1 - $x[0] ) - .5 * ( 1 - $x[0] );

	if ( $sg1 != 0 )
	{
		if ( $sg1 > 0 )
		{
			$restriction1 -= sqrt ( abs ( .5 * ( 1 - $x[0] ) - ( 1 - $x[0] )**2 ) );
		}
		else
		{
			$restriction1 += sqrt ( abs ( .5 * ( 1 - $x[0] ) - ( 1 - $x[0] )**2 ) );
		}
	}
	if ( $sg2 != 0 )
	{
		if ( $sg2 > 0 )
		{
			$restriction2
				-= sqrt ( abs ( .25 * sqrt ( 1 - $x[0] ) - .5 * ( 1 - $x[0] ) ) );
		}
		else
		{
			$restriction2
				+= sqrt ( abs ( .25 * sqrt ( 1 - $x[0] ) - .5 * ( 1 - $x[0] ) ) );
		}
	}

	my $violated = 0;
	if ( $restriction1 < 0 )
	{
		$violated += ( $gen + 1 ) * 3;
	}
	if ( $restriction2 < 0 )
	{
		$violated += ( $gen + 1 ) * 3;
	}
	return $violated;
}

sub reproduce #( neighbours, population, CR, F, SIG, individual, dimensions, algorithm ) -> ( son )
{
	my ( $neighborhood, $population, $CR, $F, $SIG, $individual, $dimensions,
		 $algorithm )
		= @_;

	my @population   = @{$population};
	my @neighborhood = @{ $neighborhood->[0] };

	my $PR = 1 / scalar (@population);

	my @prev_individual = @{ $population[$individual] };

	my @son = ();
	my $alele;

	my $sigma = 1 / $SIG;

	my $mother = $neighborhood[ rand @neighborhood ][2];

	#while ( not defined $mother)
	#{
	#    $mother =  $neighborhood[rand @neighborhood ][2];
	#}

	my $parent1 = $neighborhood[ rand @neighborhood ][2];

	#while ( not defined $parent1 or $parent1 == $mother)
	#{
	#    $parent1 = $neighborhood[rand @neighborhood ][2];
	#}
	my $parent2 = $neighborhood[ rand @neighborhood ][2];

	#while ( not defined $parent2 or $parent2 == $mother or $parent2 == $parent1)
	#{
	#    $parent2 = $neighborhood[rand @neighborhood ][2];
	#}

	print
		"Se han elegido a los individos ($mother,$parent1,$parent2) como padres\n";

	my $random_index = int ( rand ( $dimensions - 1 ) );
	for ( 0 .. $dimensions - 1 )
	{
		# Cross

		my $prob = rand (1);

		if ( $prob <= $CR or $_ == $random_index )
		{
			$alele = $population[$mother][$_]
				+ $F * ( $population[$parent1][$_] - $population[$parent2][$_] );
		}
		elsif ( $prob > $CR )
		{
			$alele = $prev_individual[$_];
		}

		if ( $_ > 0 and $algorithm ne "zdt3" )
		{
			$sigma = 4 / $SIG;
		}

		# Gauss
		$alele += random_gauss ( 0, $sigma ) if rand (1) <= $PR;

		# Meter en rango
		if ( $algorithm eq "zdt3" or $_ == 0 )
		{
			while ( $alele < 0 or $alele > 1 )
			{
				$alele = -$alele    if $alele < 0;
				$alele = 2 - $alele if $alele > 1;
			}
		}
		else
		{
			while ( $alele < -2 or $alele > 2 )
			{
				$alele = -2 - $alele if $alele < -2;
				$alele = 4 - $alele  if $alele > 2;
			}
		}
		push @son, $alele;
	}

	return @son;
}

sub random_gauss    #( media, sigma ) -> ( gauss )
{
	my ( $media, $sigma ) = @_;

	my $u1 = rand;
	my $u2 = rand;

	my $z1 = sqrt ( -2 * log ($u1) ) * cos ( 2 * pi * $u2 );
	my $z2 = sqrt ( -2 * log ($u1) ) * sin ( 2 * pi * $u2 );

	my $gauss = $media + $sigma * $z1;
	return $gauss;

}

sub tchebycheff #( son, z, neighbour, population, F, algorithm, individual) -> ( is_upgradeable)
{
	my ( $son, $z, $neighbour, $population, $F, $algorithm, $individual ) = @_;

	my @neighbour = @{$neighbour};

	my @son = @{$son};
	my @z   = @$z;

	my @population = @$population;
	my @F          = @$F;

	my $tchebycheff_son_0 = $neighbour[0] * abs ( $F[0] - $z[0] );
	my $tchebycheff_son_1 = $neighbour[1] * abs ( $F[1] - $z[1] );

	my ( $f1, $f2 )
		= select_objective_function ( $algorithm,
									  [ @{ $population[ $neighbour[2] ] } ] );

	my $tchebycheff_neighbour_0 = $neighbour[0] * abs ( $f1 - $z[0] );
	my $tchebycheff_neighbour_1 = $neighbour[1] * abs ( $f2 - $z[1] );

	my $tchebycheff_son
		= $tchebycheff_son_0 > $tchebycheff_son_1
		? $tchebycheff_son_0
		: $tchebycheff_son_1;
	my $tchebycheff_neighbour
		= $tchebycheff_neighbour_0 > $tchebycheff_neighbour_1
		? $tchebycheff_neighbour_0
		: $tchebycheff_neighbour_1;

	my $is_upgradeable = $tchebycheff_son <= $tchebycheff_neighbour;

	print
		"Se harán estas cuentas:\nmax($neighbour[0] * abs( $F[0] - $z[0]),$neighbour[1] * abs( $F[1] - $z[1])\nmax( $neighbour[0] * abs( $f1 - $z[0]), $neighbour[1] * abs( $f2 - $z[1]))\n";
	if ($is_upgradeable)
	{
		print ("Individuo $individual está actualizando al vecino $neighbour[2]\n");
		print "\n";
	}
	return $is_upgradeable;
}

1;
