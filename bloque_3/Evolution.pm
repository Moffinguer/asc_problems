#/usr/bin/perl

package Evolution;

use strict;
use warnings;
use Data::Dumper;
use File::Path qw(make_path);

use Math::Trig;

sub initialization #(population, upper_limit, lower_limit, dimensions) -> (population_list)
{
	my ( $n, $upper_limit, $lower_limit, $dimensions, $algorithm ) = @_;
	my $population_list = [];
	my $chromosome      = [];
	my $alelo           = ();
	if ( $algorithm eq "zdt3" )
	{
		for ( 0 .. $n - 1 )
		{
			for ( 0 .. $dimensions - 1 )
			{
				$alelo = $lower_limit + rand ($upper_limit);
				push @{$chromosome}, $alelo;
			}
			push @{$population_list}, $chromosome;
			$chromosome = [];
		}
	}
	else
	{
		for ( 0 .. $n - 1 )
		{
			$alelo = $lower_limit + rand ($upper_limit);
			push @{$chromosome}, $alelo;
			for ( 1 .. $dimensions - 1 )
			{
				$alelo = rand (4) - 2;
				push @{$chromosome}, $alelo;
			}
			push @$population_list, $chromosome;
			$chromosome = [];
		}
	}

	#print "GENERATED $n INDIVIDUALS WITH $dimensions CHROMOSOMES\n";

	return $population_list;
}

sub window    #(population) -> (lambda_window)
{
	my $p             = shift;
	my $step          = 1 / ( $p - 1 );
	my $lambda_window = ();
	my $lambd         = 0;
	for ( 0 .. $p - 1 )
	{
		## Add a list of N values, by each subproblem + 1 with the index of the individual
		push @{$lambda_window}, [ $lambd, 1 - $lambd, $_ ];
		$lambd += $step;
	}
	return $lambda_window;
}

sub euclidean_distance    #( $x_00, $x_10, $x_01, $x_11 )
{
	my ( $x_00, $x_01, $x_10, $x_11 ) = @_;
	my $dist = sqrt ( ( $x_00 - $x_01 )**2 + ( $x_10 - $x_11 )**2 );
	return $dist;
}

sub select_subproblems    #( index, vecinity, population, lambda_window )
{
	my ( $index, $vecinity, $population, $lambda_window ) = @_;

	my $number_near = int ( $population * $vecinity );

#	print
#		"SELECT $number_near NEIGHBOURS NEAR $index LAMBDA[$lambda_window->[$index]->[0], $lambda_window->[$index]->[1]]\n";

	my @new_window = sort {
		euclidean_distance ( $lambda_window->[$index][0], $a->[0],
							 $lambda_window->[$index][1], $a->[1]
		) <=> euclidean_distance ( $lambda_window->[$index][0],
									   $b->[0], $lambda_window->[$index][1],
									   $b->[1] );
	} @{$lambda_window};

	#print( Dumper(@new_window) . "\n" );

	@new_window = splice ( @new_window, 0, $number_near );

	return \@new_window;
}

sub zdt3    #(x, n) -> (f1(x), f2(x))
{
	my ( $x, $n ) = @_;

	my $f_1 = $x->[0];
	my $g;

	$g += $x->[$_] for ( 1 .. $n - 1 );

	$g = 9 * $g / ( $n - 1 ) + 1;
	my $h = $f_1 / $g;
	$h = 1 - sqrt ($h) - $h * sin ( 10 * $f_1 * pi );

	my $f_2 = $g * $h;

	return ( $f_1, $f_2 );
}

sub cf6    #(x, n) -> (f1(x), f2(x))
{

	my ( $x,   $n )   = @_;
	my ( $f_1, $f_2 ) = ( $x->[0], ( 1 - $x->[0] )**2 );

	for ( 1 .. $n - 1 )
	{
		if ( $_ % 2 == 1 )
		{
			$f_1
				+= ( $x->[$_] - .8 * $x->[0] * cos ( 6 * pi * $x->[0] + $_ * pi / $n ) )
				**2;
		}
		else
		{
			$f_2
				+= ( $x->[$_] - .8 * $x->[0] * sin ( 6 * pi * $x->[0] + $_ * pi / $n ) )
				**2;
		}
	}
	return ( $f_1, $f_2 );
}

sub problem_function    #(x, n, problem) -> (f1(x), f2(x))
{
	my ( $x, $n, $problem ) = @_;

	if ( $problem eq "zdt3" )
	{
		return zdt3 ( $x, $n );
	}
	else
	{
		return cf6 ( $x, $n );
	}
}

sub constraint_sigma_evolution    #($gen, $growth_factor) -> ( $sigma )
{
	my ( $gen, $growth_factor ) = @_;
	return 1 / ( 1 + exp ( -$growth_factor * $gen ) );

}

sub constraints # ($problem, $sol, population, generation) -> (f_1_restricted, f_2_restricted)
{
	my ( $problem, $sol, $population, $generation ) = @_;

	return 0 if $problem eq "zdt3";

	my $restriction1
		= $sol->[1]
		- .8 * $sol->[0] *
		sin ( 6 * pi * $sol->[0] + 2 * pi / scalar (@$population) );
	my $restriction2
		= $sol->[3]
		- .8 * $sol->[0] *
		sin ( 6 * pi * $sol->[0] + 4 * pi / scalar (@$population) );

	my $sg1 = .5 * ( 1 - $sol->[0] ) - ( 1 - $sol->[0] )**2;
	my $sg2 = .25 * sqrt ( 1 - $sol->[0] ) - .5 * ( 1 - $sol->[0] );

	if ( $sg1 != 0 )
	{
		$restriction1
			-= $sg1 / abs ($sg1) *
			sqrt ( abs ( .5 * ( 1 - $sol->[0] ) - ( 1 - $sol->[0] )**2 ) );
	}
	if ( $sg2 != 0 )
	{
		$restriction2
			-= $sg2 / abs ($sg2) *
			sqrt ( abs ( .25 * sqrt ( 1 - $sol->[0] ) - .5 * ( 1 - $sol->[0] ) ) );
	}

	my $restricted = 0;
	if ( $restriction1 < 0 )
	{
		$restricted += constraint_sigma_evolution ( $generation, 0.1 );
	}
	if ( $restriction2 < 0 )
	{
		$restricted += constraint_sigma_evolution ( $generation, 0.2 );
	}

	return $restricted;
}

sub read_input # (filename) -> ( population, generations, neighborhood, inferior_limit, upper_limit, dimensions )
{
	my $filename = shift;
	my ( $population,  $generations, $neighborhood, $inferior_limit,
		 $upper_limit, $dimensions,  $experiments,  $algorithm
	) = ( 100, 100, 0.25, 0, 1, 30, 1, "zdt3" );

	$filename = "./INPUT_FILES/$filename.in";
	unless ( -e $filename )
	{
		print "$filename does not exists\n";
		return ( $population,  $generations, $neighborhood, $inferior_limit,
				 $upper_limit, $dimensions,  $experiments,  $algorithm );
	}
	my $folder = open ( my $input, '<', $filename ) or do
	{
		print "Error opening $filename: $!";
		exit 1;
	};

	while (<$input>)
	{
		unless ( $_
			=~ /^\s*(algorithm|experiments|population|generations|dimensions|upperLimit|inferiorLimit|neighborhood)\=([+-]?(?:\d*\.\d+|\d+\.\d*|\d+)|zdt3|cf6)\s*$/x
			)
		{
			print "\'$_\' bad formatted\n";
			next;
		}

		if ( $1 eq "population" )
		{
			$population = 0 + $2;
			next;
		}
		if ( $1 eq "experiments" )
		{
			$experiments = 0 + $2;
			next;
		}
		if ( $1 eq "generations" )
		{
			$generations = 0 + $2;
			next;
		}
		if ( $1 eq "dimensions" )
		{
			$dimensions = 0 + $2;
			next;
		}
		if ( $1 eq "upperLimit" )
		{
			$upper_limit = 0 + $2;
			next;
		}
		if ( $1 eq "inferiorLimit" )
		{
			$inferior_limit = 0 + $2;
			next;
		}
		if ( $1 eq "neighborhood" )
		{
			$neighborhood = 0 + $2;
			next;
		}
		if ( $1 eq "algorithm" )
		{
			$algorithm = $2;
		}
	}
	close ($input);
	return ( $population,  $generations, $neighborhood, $inferior_limit,
			 $upper_limit, $dimensions,  $experiments,  $algorithm );
}

sub mutate #($x, $window, $upper_limit, $lower_limit, $dimensions, $population)
{
	my ( $x, $window, $upper_limit, $lower_limit, $dimensions, $population,
		 $algorithm )
		= @_;
	my ( $CR, $F, $PR, $SIG )
		= ( .5, .5, 1 / $population, 20 )
		;    #F and SIG should be added as a parameter to try things
	my $chromosome_mutated;
	my $info;

	my $mother = $window->[ int rand @$window ]->[2];
	my $parent1;
	my $parent2;
	$parent1 = $window->[ int rand @$window ]->[2]
		while not defined $parent1
		or $parent1 == $mother;
	$parent2 = $window->[ int rand @$window ]->[2]
		while not defined $parent2
		or $parent2 == $mother
		or $parent2 == $parent1;

	my $parents = [ $mother, $parent1, $parent2 ];
	my $sigma   = $upper_limit - $lower_limit;
	$sigma /= $SIG;

	my $forced_mutant = int ( rand ( $dimensions - 1 ) );

	for ( 0 .. $dimensions - 1 )
	{
		# Cross
		if ( rand () <= $CR or $_ == $forced_mutant )
		{
			$info
				= $population->[ $parents->[0] ]->[$_]
				+ $F
				* (
				$population->[ $parents->[1] ]->[$_] - $population->[ $parents->[2] ]->[$_] );
		}
		else
		{
			$info = $x->[$_];
		}

		# Gauss
		#	$info += exp ( -.5 * ( $info / $sigma )**2 ) / ( $sigma * sqrt ( 2 * pi ) )
		#	if rand () <= $PR;

		push @$chromosome_mutated, $info;
	}

	if ( $algorithm eq "zdt3" )
	{
		for my $i ( 0 .. scalar @$chromosome_mutated - 1 )
		{
			while (    $chromosome_mutated->[$i] < $lower_limit
					or $chromosome_mutated->[$i] > $upper_limit )
			{
				if ( $chromosome_mutated->[$i] < $lower_limit )
				{
					$chromosome_mutated->[$i]
						= $lower_limit + ( $lower_limit - $chromosome_mutated->[$i] );
				}
				if ( $chromosome_mutated->[$i] > $upper_limit )
				{
					$chromosome_mutated->[$i]
						= $upper_limit - ( $chromosome_mutated->[$i] - $upper_limit );
				}
			}
		}
	}
	else
	{
		while (    $chromosome_mutated->[0] < $lower_limit
				or $chromosome_mutated->[0] > $upper_limit )
		{
			if ( $chromosome_mutated->[0] < $lower_limit )
			{
				$chromosome_mutated->[0]
					= $lower_limit + ( $lower_limit - $chromosome_mutated->[0] );
			}
			if ( $chromosome_mutated->[0] > $upper_limit )
			{
				$chromosome_mutated->[0]
					= $upper_limit - ( $chromosome_mutated->[0] - $upper_limit );
			}
		}
		my $cf6_lower_limit = -2;
		my $cf6_upper_limit = 2;
		for my $i ( 1 .. scalar @$chromosome_mutated - 1 )
		{
			while (    $chromosome_mutated->[$i] < $cf6_lower_limit
					or $chromosome_mutated->[$i] > $cf6_upper_limit )
			{
				if ( $chromosome_mutated->[$i] < $cf6_lower_limit )
				{
					$chromosome_mutated->[$i]
						= $cf6_lower_limit + ( $cf6_lower_limit - $chromosome_mutated->[$i] );
				}
				if ( $chromosome_mutated->[$i] > $cf6_upper_limit )
				{
					$chromosome_mutated->[$i]
						= $cf6_upper_limit - ( $chromosome_mutated->[$i] - $cf6_upper_limit );
				}
			}
		}
	}
	return $chromosome_mutated;
}

my ( $population,  $generations, $neighborhood, $inferior_limit,
	 $upper_limit, $dimensions,  $experiments,  $algorithm
) = read_input $ARGV[0];
print "####################################\n";
print "Variable\tValor\n";
print "ALGORITHM\t$algorithm\n";
print "POPULATION\t$population\n";
print "GENERATIONS\t$generations\n";
print "NEIGHBORHOOD\t$neighborhood\n";
print "INFERIOR_LIMIT\t$inferior_limit\n";
print "UPPER_LIMIT\t$upper_limit\n";
print "DIMENSIONS\t$dimensions\n";
print "EXPERIMENTS\t$experiments\n";
print "###################################\n";

my ( $lambda_window,       $population_list,
	 $f_1,                 $f_2,
	 $evaluated_functions, $z_1,
	 $z_2,                 $times,
	 $lambda,              $reproducted_evaluated_functions,
	 $new_sols,            $time_required,
	 $f_1t,                $f_2t,
	 $tchebycheff1,        $tchebycheff2,
	 $violations
);

my $directory_path = "./results/ZDT3";
$directory_path = "./results/CF6" if $algorithm eq "cf6";
my $output_file;

for my $exec ( 0 .. $experiments - 1 )
{
	open my $gnuplot, '|-', 'gnuplot -persistent' or do
	{
		print "Graphic for exec $exec not created $!";
		next;
	};
	if ( $algorithm eq "zdt3" )
	{
		print $gnuplot <<'GNUPLOT_SCRIPT';
		set title "Real time graph"
		set xlabel "f1"
		set ylabel "f2"
		set style data points
		set grid
		plot "-" using 1:2 title "ZDT3" with points pointtype 7 linecolor rgb 'blue'
GNUPLOT_SCRIPT
	}
	else
	{
		print $gnuplot <<'GNUPLOT_SCRIPT';
		set title "Real time graph"
		set xlabel "f1"
		set ylabel "f2"
		set style data points
		set grid
		plot "-" using 1:2 title "CF6" with points pointtype 7 linecolor rgb 'blue'
GNUPLOT_SCRIPT
	}

	print "...................\nEXECUTION $exec...\n";

	unless ( -d $directory_path )
	{
		make_path ($directory_path) or do
		{
			print "Unable to create $directory_path: $!\n";
			exit 1;
		};
	}
	$output_file = "./results/ZDT3/$ARGV[0]$exec.out";
	unless ( -e $output_file )
	{
		open ( my $create_file, '>', $output_file ) or do
		{
			print "Error opening $output_file: $!\n";
			exit 1;
		};
		close $create_file;
	}
	else
	{
		open ( my $clear_file, '>', $output_file ) or do
		{
			print "Error cleaning $output_file: $!\n";
			exit 1;
		};
		close $clear_file;
	}

	$times  = ();
	$lambda = window $population;
	( $z_1, $z_2 ) = ( undef, undef );

	#print ( "Lambdas:\n" . Dumper ($lambda) . "\n" );

	#X[i]
	$population_list
		= initialization ( $population, $upper_limit, $inferior_limit,
						   $dimensions, $algorithm );

	for my $gen ( 0 .. $generations - 1 )
	{
		print "...................\nGENERATION $gen...\n";
		$time_required = time;

		## Initial Evaluation
		for my $x ( @{$population_list} )
		{
			( $f_1, $f_2 ) = problem_function ( $x, $dimensions, $algorithm );
			$violations = constraints ( $algorithm, $x, $population_list, $gen );
			( $f_1, $f_2 ) = ( $f_1 - $violations, $f_2 - $violations );
			$z_1 = $f_1 if not defined $z_1 or $f_1 < $z_1;
			$z_2 = $f_2 if not defined $z_2 or $f_2 < $z_2;
		}

		#print ( "Population:\n" . Dumper ($population_list) . "\n" );

		print ( "Functions:\n" . Dumper($evaluated_functions). "\n");

		## Guardar soluciones no nominadas (Z?)
		print "Z_best=($z_1,$z_2)\n\n";

		# Time to evolve
		for my $individual ( 0 .. $population - 1 )
		{
			## Reproduction
			#lambda
			$lambda_window->{$individual}
				= select_subproblems ( $individual, $neighborhood, $population, $lambda )
				if not exists $lambda_window->{$individual};

			# print ( "Subproblems:\n" . Dumper($lambda_window). "\n");

			$new_sols = mutate $population_list->[$individual],
				$lambda_window->{$individual},
				$upper_limit, $inferior_limit, $dimensions, $population_list, $algorithm;

			print ( "Mutated individual $individual:\n" . Dumper ($new_sols) . "\n" );

			# Restrictions violated

			##Evaluation
			( $f_1, $f_2 ) = problem_function ( $new_sols, $dimensions, $algorithm );
			$violations = constraints ( $algorithm, $new_sols, $population_list, $gen );

			( $f_1, $f_2 ) = ( $f_1 - $violations, $f_2 - $violations );
			print "F(y_$individual) GEN:$gen =($f_1,$f_2)\n";

			open ( my $file_handle, ">>", $output_file ) or do
			{
				print "Error opening $output_file: $!\n";
				exit 1;
			};
			print $file_handle "$f_1\t$f_2\t$violations\n";
			close $file_handle;

			##Update Best Sol
			$z_1 = $f_1 if $f_1 < $z_1;
			$z_2 = $f_2 if $f_2 < $z_2;
			print "Z_best GEN:$gen =($z_1,$z_2)\n";

			##Update Neigbours
			#	0->[
			#		A_0,A_1 => $j->[0], $j->[1]
			#	(0,1        ,0),
			#	(0.5,0.5    ,1),
			#	(1,0        ,2)
			#]

			for my $j ( @{ $lambda_window->{$individual} } )
			{
				$tchebycheff1 = $j->[0] * abs ( $f_1 - $z_1 );
				my $tch11 = $j->[1] * abs ( $f_2 - $z_2 );
				$tchebycheff1 = $tch11 if $tchebycheff1 < $tch11;

				( $f_1t, $f_2t )
					= problem_function ( $population_list->[ $j->[2] ],
										 $dimensions, $algorithm );

				$tchebycheff2 = $j->[0] * abs ( $f_1t - $z_1 );
				my $tch22 = $j->[1] * abs ( $f_2t - $z_2 );
				$tchebycheff2 = $tch22 if $tchebycheff2 < $tch22;

				next if $tchebycheff1 > $tchebycheff2;

				print(Dumper($individual,$j). "\n");
				$population_list->[ $j->[2] ] = $new_sols;
				print ("Actualizado el vecino $j->[2] del individuo $individual");
			}
			$tchebycheff1 = undef;
			if ( $gen == $generations - 1 )
			{
				print $gnuplot "$f_1 $f_2\n";
				#print ("Ultima generación: $f_1\t$f_2\n");
			}
			## Update EP
		}

		#close $gnuplot;
		#sleep(5);

		$time_required = time - $time_required;

		push @{$times}, $time_required;
		$reproducted_evaluated_functions = ();

		#	print "GEN $gen TIME: $time_required\n";

		$evaluated_functions             = [];
		$reproducted_evaluated_functions = [];

		#$z_1                             = undef;
		#$z_2                             = undef;
	}

	my $average_time = 0;
	$average_time += $_ for (@$times);
	$average_time /= scalar @$times;

	#	print "AVERAGE TIME TO EXECUTE EACH GEN IS $average_time\n";

	(  $lambda_window,       $population_list,
	   $f_1,                 $f_2,
	   $evaluated_functions, $z_1,
	   $z_2,                 $times,
	   $lambda,              $reproducted_evaluated_functions,
	   $new_sols,            $time_required
		)
		= ( undef, undef, undef, undef, undef, undef, undef, undef,
			undef, undef, undef, undef, undef, undef );

	( $f_1t, $f_2t ) = ( undef, undef );

	close $gnuplot;
}

1;
