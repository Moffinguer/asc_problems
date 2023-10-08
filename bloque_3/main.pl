#/usr/bin/perl

sub read_input # (filename) -> ( population, generations, neighborhood, inferior_limit, upper_limit, dimensions )
{
	my $filename = shift;
	my ( $population, $generations, $neighborhood, $inferior_limit, $upper_limit, $dimensions ) = ( 100, 100, 0.25, 0, 1, 30 );

	$filename = "./INPUT_FILES/$filename.in";
	unless( -e $filename )
	{
		print "$filename does not exists\n";
		return ($population, $generations, $neighborhood, $inferior_limit, $upper_limit, $dimensions);
	}
	open( my $input, '<', $filename ) or die "Error opening $filename: $!";
	while ( <$input> ) {
		unless ( $_ =~ /^\s*(population|generations|dimensions|upperLimit|inferiorLimit|neighborhood)\=([+-]?(?:\d*\.\d+|\d+\.\d*|\d+))\s*$/ )
		{
			print "\'$_\' bad formatted\n";
			next;
		}

		if ( $1 eq "population" )
		{
			$population = 0 + $2;
		}
		elsif ( $1 eq "generations" )
		{
			$generations = 0 + $2;
		}
		elsif ( $1 eq "dimensions" )
		{
			$dimensions = 0 + $2;
		}
		elsif( $1 eq "upperLimit" )
		{
			$upper_limit = 0 + $2;
		}
		elsif( $1 eq "inferiorLimit" )
		{
			$inferior_limit = 0 + $2;
		}
		elsif( $1 eq "neighborhood" )
		{
			$generations = 0 + $2;
		}
        }
        close($input);
	return ($population, $generations, $neighborhood, $inferior_limit, $upper_limit, $dimensions);
}


my ($population, $generations, $neighborhood, $inferior_limit, $upper_limit, $dimensions) = read_input $ARGV[0];
print "####################################\n";
print "Variable\tValor\n";
print "POPULATION\t$population\n";
print "GENERATIONS\t$generations\n";
print "NEIGHBORHOOD\t$neighborhood\n";
print "INFERIOR_LIMIT\t$inferior_limit\n";
print "UPPER_LIMIT\t$upper_limit\n";
print "DIMENSIONS\t$dimensions\n";
print "###################################\n";


sub initialization #(population, upper_limit, lower_limit, dimensions) -> (population_list)
{
	my ($n, $upper_limit, $lower_limit, $dimensions) = @_;
	my $population_list = ();
	my $chromosome = ();
	my $alelo = ();
	for (0..$n-1)
	{
		for (0..$dimensions-1)
		{
			$alelo = $lower_limit + rand($upper_limit);
			push @{$chromosome}, $alelo;
		}
		push @{$population_list}, $chromosome;
	}
	print "GENERATED $population INDIVIDUALS WITH $dimensions CHROMOSOMES\n";

	return $population_list;
}
sub window #(population) -> (lambda_window)
{
	my $p = shift;
	my $step = 1 / ($p - 1);
	my $lambda_window = ();
	my $lambd = 0;
	for (0..$p-1)
	{
		push @{$lambda_window}, ( $lambd, 1 - $lambd );
		$lambd += $step;
	}
	return $lambda_window;
}
sub select_subproblems #( index, vecinity, population, lambda_window )
{
	my ($index, $vecinity, $population, $lambda_window) = @_;

	my $number_near = int( $population * $vecinity );
	print "SELECT $number_near NEIGHBOURS NEAR $index LAMBDA\n";

	my @new_winddow = sort {
		my $distancia_a = sqrt(
		($lambda_window[$indice][0] - $a[0]) ** 2 +
		($lambda_window[$indice][1] - $a[1]) ** 2
		);

		my $distancia_b = sqrt(
		($lambda_window[$indice][0] - $b[0]) ** 2 +
		($lambda_window[$indice][1] - $b[1]) ** 2
		);

		$distancia_a <=> $distancia_b;
	} @{$lambda_window};
	@new_window = splice(@new_window, 0, $number_near);
	return \@new_window;
}

sub zdt3 #(x, n, upperLimit, lowerLimmit) -> (f1(x), f2(x))
{
	my ( $x, $n, $index, $upper_limit, $lower_limit ) = @_;


	my $f_1 = $x[0];
	my $g = 0;
	for (1..$n-1)
	{
		$g += $x[$_]
	}
	$g = 1 + 9/($n-1) * $g;
	my $h = $f_1/$g;

	import Math::Trig;
	$h = 1 - ( sqrt($h) + $h * sin(10 * Math::Trig::pi * $f_1) );
	my $f_2 = $g * $h;

	#	$f_1 += ( $f_1 > $upper_limit ) * ( -$upper_limit + $f_1 ) + ( $f_1 < $lower_limit ) * ( -$f_1 + $lower_limit ) while $f_1 > $upper_limit or $f_1 < $lower_limit;

	#       $f_2 += ( $f_2 > $upper_limit ) * ( -$upper_limit + $f_2 ) + ( $f_2 < $lower_limit ) * ( -$f_2 + $lower_limit ) while $f_2 > $upper_limit or $f_2 < $lower_limit;

	return ($f_1, $f_2);
}


my ( $lambda_window, $population_list, $f_1, $f_2, $evaluated_functions, $z_1, $z_2, $times, $lambda, $reproducted_evaluated_functions, $new_sols );
$z_1 = undef;
$z_2 = undef;
my $time_required;
$times = ();

for my $gen (0.. $generations-1){
	print "...................\nGENERATION $gen\n\n";
	$time_required = time;
	$lambda = window $population;

	#X[i]
	$population_list = initialization ($population, $upper_limit, $lower_limit, $dimensions);

	## Initial Evaluation
	for my $x (@{$population_list})
	{
		($f_1, $f_2) = zdt3($x, $dimensions, $upper_limit, $lower_limit);
		push @{evaluated_functions}, ($f_1, $f_2);
	}
	for my $f_minimum (@{$evaluated_functions})
	{
		$z_1 = $f_minimum if $z_1 == undef or $f_minimum[0] < $z_1;
		$z_2 = $f_minimum if $z_2 == undef or $f_minimum[1] < $z_2;
	}
	## Guardar soluciones no nominadas (Z?)

	print "Z=($z_1,$z_2)\n";

	# Time to evolve
	for my $individual (0..$population-1)
	{
		## Reproduction
		#lambda
		$lambda_window = select_subproblems($individual, $neighborhood, $population, $lambda);
		## $new_sols =;

		##Evaluation
		for my $x (@{$new_sols})
		{
			($f_1, $f_2) = zdt3($x, $dimensions, $upper_limit, $lower_limit);
			push @{$reproducted_evaluated_functions}, ($f_1, $f_2);
		}

		##Update Best Sol
		for my $f_minimum (@{$reproducted_evaluated_functions})
		{
			$z_1 = $f_minimum if $f_minimum[0] < $z_1;
			$z_2 = $f_minimum if $f_minimum[1] < $z_2;
		}

		##Update Neigbours

		## Update EP

	}

	$time_required = time - $time_required;

	push @{$times}, $time_required;
	($z_1,$z_2) = (undef, undef);
	$reproducted_evaluated_functions=();

	print "GEN TIME: $time_required";

}


sub reproduction
{



	## Mutate
}

