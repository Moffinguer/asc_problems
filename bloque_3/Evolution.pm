#/usr/bin/perl

package Evolution;

use strict;
use warnings;

use Data::Dumper;

use File::Path qw(make_path);

use lib "modules";
use Objectives;
use Window;

sub create_folder_struct #( directory_path, output_file, filename, exec, algorithm) -> ()
{
	my( $directory_path, $output_file, $filename, $exec, $algorithm) = @_;

	unless ( -d $directory_path )
	{
		make_path ($directory_path) or do
		{
			print "Unable to create $directory_path: $!\n";
			exit 1;
		};
	}

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

	return;

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

## Leer fichero texto
my ( $population,  $generations, $neighborhood, $inferior_limit,
	 $upper_limit, $dimensions,  $experiments,  $algorithm
) = read_input $ARGV[0];


my( $CR, $F, $SIG) = (.5, .5, 20);

my $directory_path = "./results/ZDT3";
$directory_path = "./results/CF6" if $algorithm eq "cf6";



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


my $gnuplot;

## Iterar sobre experimentos
for my $exec ( 0 .. $experiments - 1 )
{
	my $output_file = "./results/ZDT3/$ARGV[0]$exec.out";
	$output_file = "./results/CF6/$ARGV[0]$exec.out" if $algorithm eq "cf6";

	create_folder_struct $directory_path, $output_file, $ARGV[0], $exec;
	open ( my $file_handle, ">>", $output_file ) or do
			{
				print "Error opening $output_file: $!\n";
				exit 1;
			};


	print("EXPERIMENT ".($exec+1). "...\n");

	# Crear Lambdas + Crear Ventanas por subproblema
	my %lambda_window = Window::window $population, $neighborhood;

	my @z = (undef, undef); my ( $f1, $f2);
	
	# Inicializar poblacion
	my @population_list = Objectives::init_population($population, $upper_limit, $inferior_limit, $dimensions, $algorithm);
	for my $fenotype (@population_list)
	{
		# Evaluar Prestaciones
		($f1, $f2) = Objectives::select_objective_function($algorithm, [@$fenotype]);
		#print("Valores de F=($f1,$f2)\n");
		# Generar Z
		@z = Objectives::best_evaluations($f1, $f2, \@z);
		#print("Z_best=($z[0],$z[1])\n");
	}
		print("\n");

	## Iterar sobre generaciones
	for my $gen ( 0 .. $generations - 1 )
	{
		print("GENERATION ". ($gen + 1)."....\n");

		## Pintar Grafica Ultima Generacion
		if ( $gen + 1 == $generations )
		{
			open $gnuplot, '|-', 'gnuplot -persistent' or do
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
		}
		
		## Iterar sobre subproblemas
		for my $individual ( 0 .. $population - 1 )
		{
			my @individual_window = $lambda_window{$individual};
			my @son = Objectives::reproduce(\@individual_window, \@population_list, $CR, $F, $SIG, $individual, $dimensions, $algorithm);
			#print("Se ha creado un hijo\n".Dumper(@son)."\n a partir del individuo $individual. Usando el algoritmo $algorithm y una ventana de: ". Dumper(@individual_window)."\n");

			my @F = Objectives::select_objective_function($algorithm, \@son);
			my $violations = Objectives::constraints(\@son, $algorithm, $gen);
			@F = ( $F[0] + $violations, $F[1] + $violations);
			#print("Valores de F=($F[0],$F[1])\n");

			print $file_handle "$F[0]\t$F[1]\t$violations\n";

			@z = Objectives::best_evaluations($F[0],$F[1] , \@z);
			print("Z_best=($z[0],$z[1])\n");

			print $gnuplot "$F[0] $F[1]\n" if  $gen + 1 == $generations;

			for my $neighbour ( @{$individual_window[0]} )
			{
				print "Evaluando al individuo $individual y su vecino $neighbour->[2]\n";
				if (Objectives::tchebycheff(\@son, \@z, \@{$neighbour}, \@population_list,\@F, $algorithm, $individual))
				{
					@{$population_list[$individual]} = @son;
				}
				else
				{
					print "Individuo $individual no pudo actualizar al vecino $neighbour->[2]\n\n"
				}
			}
		}
		print(Dumper(@population_list)."\n");
	}
	close $file_handle;
	close($gnuplot);
}

1;
