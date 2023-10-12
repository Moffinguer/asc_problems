#/usr/bin/perl

package Evolution;

use strict;
use warnings;
use Data::Dumper;
use File::Path qw(make_path);

use Math::Trig;

sub initialization #(population, upper_limit, lower_limit, dimensions) -> (population_list)
{
	my ( $n, $upper_limit, $lower_limit, $dimensions ) = @_;
	my $population_list = [];
	my $chromosome      = [];
	my $alelo           = ();
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
	print "GENERATED $n INDIVIDUALS WITH $dimensions CHROMOSOMES\n";

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

sub read_input # (filename) -> ( population, generations, neighborhood, inferior_limit, upper_limit, dimensions )
{
	my $filename = shift;
	my ( $population,  $generations, $neighborhood, $inferior_limit,
		 $upper_limit, $dimensions,  $experiments
	) = ( 100, 100, 0.25, 0, 1, 30, 1 );

	$filename = "./INPUT_FILES/$filename.in";
	unless ( -e $filename )
	{
		print "$filename does not exists\n";
		return ( $population,     $generations, $neighborhood,
				 $inferior_limit, $upper_limit, $dimensions );
	}
	my $folder = open ( my $input, '<', $filename ) or do
	{
		print "Error opening $filename: $!";
		exit 1;
	};

	while (<$input>)
	{
		unless ( $_
			=~ /^\s*(experiments|population|generations|dimensions|upperLimit|inferiorLimit|neighborhood)\=([+-]?(?:\d*\.\d+|\d+\.\d*|\d+))\s*$/x
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
	}
	close ($input);
	return ( $population,  $generations, $neighborhood, $inferior_limit,
			 $upper_limit, $dimensions,  $experiments );
}

sub mutate #($x, $window, $upper_limit, $lower_limit, $dimensions, $population)
{
	my ( $x, $window, $upper_limit, $lower_limit, $dimensions, $population ) = @_;
	my ( $CR, $F, $PR, $SIG )
		= ( .5, .5, 1 / $population, 20 )
		;    #F and SIG should be added as a parameter to try things
	my $chromosome_mutated;
	my $info;
	my $parents = [ $window->[ rand @$window ]->[2],
					$window->[ rand @$window ]->[2],
					$window->[ rand @$window ]->[2]
	];
	my $sigma = $upper_limit - $lower_limit;
	$sigma /= $SIG;

	for ( 0 .. $dimensions - 1 )
	{
		# Cross
		if ( rand () <= $CR )
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
		$info += exp ( -.5 * ( $info / $sigma )**2 ) / ( $sigma * sqrt ( 2 * pi ) )
			if rand () <= $PR;

		push @$chromosome_mutated, $info;
	}

	for my $i ( 0 .. scalar @$chromosome_mutated - 1 )
	{
		while (    $chromosome_mutated->[$i] < $lower_limit
				or $chromosome_mutated->[$i] > $upper_limit )
		{
			if ( $chromosome_mutated->[$i] <= $lower_limit )
			{
				$chromosome_mutated->[$i]
					= $lower_limit + ( $lower_limit - $chromosome_mutated->[$i] );
			}
			if ( $chromosome_mutated->[$i] >= $upper_limit )
			{
				$chromosome_mutated->[$i]
					= $upper_limit - ( $chromosome_mutated->[$i] - $upper_limit );
			}
		}
	}

	return $chromosome_mutated;
}

my ( $population,  $generations, $neighborhood, $inferior_limit,
	 $upper_limit, $dimensions,  $experiments
) = read_input $ARGV[0];
print "####################################\n";
print "Variable\tValor\n";
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
my $output_file;

for my $exec ( 0 .. $experiments - 1 )
{
	open my $gnuplot, '|-', 'gnuplot -persistent' or do
	{
		print "Graphic for exec $exec not created $!";
		next;
	};
	print $gnuplot <<'GNUPLOT_SCRIPT';
set title "Real time graph"
set xlabel "f1"
set ylabel "f2"
set style data points
set grid
plot "-" using 1:2 title "ZDT3" with points pointtype 7 linecolor rgb 'blue'
GNUPLOT_SCRIPT

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

	#print ("Lambdas:\n" . Dumper($lambda) . "\n");

	#X[i]
	$population_list
		= initialization ( $population, $upper_limit, $inferior_limit,
						   $dimensions );

	for my $gen ( 0 .. $generations - 1 )
	{
		print "...................\nGENERATION $gen...\n";
		$time_required = time;

		## Initial Evaluation
		for my $x ( @{$population_list} )
		{
			( $f_1, $f_2 ) = zdt3 ( $x, $dimensions );
			push @{$evaluated_functions}, [ $f_1, $f_2 ];
		}

		#print ( "Population:\n" . Dumper($population_list) . "\n");
		#print ( "Functions:\n" . Dumper($evaluated_functions). "\n");

		for my $f_minimum ( @{$evaluated_functions} )
		{
			$z_1 = $f_minimum->[0] if not defined $z_1 or $f_minimum->[0] < $z_1;
			$z_2 = $f_minimum->[1] if not defined $z_2 or $f_minimum->[1] < $z_2;
		}

		## Guardar soluciones no nominadas (Z?)
		#print "Z=($z_1,$z_2)\n\n";

		# Time to evolve
		for my $individual ( 0 .. $population - 1 )
		{
			## Reproduction
			#lambda
			$lambda_window
				= select_subproblems ( $individual, $neighborhood, $population, $lambda );

			#print ( "Subproblems:\n" . Dumper($lambda_window). "\n");

			$new_sols = mutate $population_list->[$individual], $lambda_window,
				$upper_limit, $inferior_limit, $dimensions, $population_list;

			#print("Mutated individual:\n" . Dumper( $new_sols ) . "\n");

			# Restrictions violated
			$violations = 0;

			##Evaluation
			( $f_1, $f_2 ) = zdt3 ( $new_sols, $dimensions );

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

			#print "Z_best=($z_1,$z_2)\n\n";

			##Update Neigbours
			for my $j (@$lambda_window)
			{
				( $f_1t, $f_2t ) = zdt3 ( $population_list->[ $j->[2] ], $dimensions );
				$tchebycheff1 = $j->[0] * abs ( $f_1 - $z_1 );
				$tchebycheff2 = $j->[0] * abs ( $f_1t - $z_1 );
				$tchebycheff1 = $j->[1] * abs ( $f_2 - $z_2 )
					if $tchebycheff1 < $j->[1] * abs ( $f_2 - $z_2 );
				$tchebycheff2 = $j->[1] * abs ( $f_2t - $z_2 )
					if $tchebycheff2 < $j->[1] * abs ( $f_2t - $z_2 );

				$population_list->[ $j->[2] ] = $new_sols if $tchebycheff1 <= $tchebycheff2;
			}
			print $gnuplot "$f_1 $f_2\n" if $gen == $generations - 1;
			## Update EP
		}

		$time_required = time - $time_required;

		push @{$times}, $time_required;
		$reproducted_evaluated_functions = ();

		print "GEN $gen TIME: $time_required\n";

		$evaluated_functions             = [];
		$reproducted_evaluated_functions = [];

		#$z_1                             = undef;
		#$z_2                             = undef;
	}

	my $average_time = 0;
	$average_time += $_ for (@$times);
	$average_time /= scalar @$times;
	print "AVERAGE TIME TO EXECUTE EACH GEN IS $average_time\n";

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
