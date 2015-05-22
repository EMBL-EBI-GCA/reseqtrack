package ReseqTrack::Hive::Process::ConvertBedToBigBed;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists delete_file);
use ReseqTrack::Tools::GeneralUtils   qw(get_open_file_handle execute_system_command);
use Math::Trig;
use Scalar::Util::Numeric qw(isint isnum);
use File::Basename;

sub param_defaults {
  return {
     min_score   => 370,
     max_score   => 1000,
     add_names   => 1,
     rescale_scores => undef,
     name_prefix => undef,
     verbose     => 1,
     in_columns  => { seq    => 1,
                      start  => 2,
                      end    => 3,
                      name   => undef,
                      score  => 4,
                      strand => undef,
                    },
     out_columns => { seq    => 1,
                      start  => 2,
                      end    => 3,
                      name   => 4,
                      score  => 5,
                      strand => 6,
                    }, 
  };
}


sub run {
  my $self = shift @_;

  $self->param_required( 'bed' );
  my $chr_file        = $self->param_required( 'chr_file' );
  my $bedToBigBedPath = $self->param_required( 'bedToBigBedPath' );
  my $min_score       = $self->param('min_score');
  my $max_score       = $self->param('max_score');
  my $rescale_scores  = $self->param('rescale_scores');
  my $in_columns      = $self->param('in_columns');
  my $out_columns     = $self->param('out_columns');
  my $add_names       = $self->param('add_names');
  my $name_prefix     = $self->param('name_prefix');
  my $verbose         = $self->param('verbose');
  my $beds            = $self->param_as_array( 'bed' );

  foreach my $bed ( @$beds ) {
    check_file_exists( $bed );
  }
  
  check_file_exists(  $chr_file );
 
  throw("Rescale scores must be a number, or iqr")
  if ( defined $rescale_scores
    && !( isnum($rescale_scores) || $rescale_scores eq 'iqr' ) );
  
  
  _check_column_mapping_indices( $in_columns );
  _check_column_mapping_indices( $out_columns );
  _check_output_mapping_possible( $in_columns, $out_columns, $add_names );
  
  my @columns_in  = _column_mappings( $in_columns );
  my @columns_out = _column_mappings( $out_columns );

  my @output_files;

  foreach my $in_file ( @$beds ) {

    unless ($name_prefix) {
      $name_prefix = basename($in_file);
      $name_prefix =~ s/.gz$//;
      $name_prefix =~ s/.bed$//;
      $name_prefix .= '_';
    }

    my $out_file = $in_file;
    $out_file =~ s/\.bed(?:\.gz)?$/\.bb/;
    push ( @output_files, $out_file);

    my $tmp_file_name = $out_file . '.tmp';

    my ( $q1, $q3 ) =
    _peak_score_interquartile_range( $in_file, $$in_columns{score}, $verbose )
      if ( defined $rescale_scores && $rescale_scores eq 'iqr' );

    _transform_input_file(
      $in_file,     $tmp_file_name, $rescale_scores, $add_names,
      $name_prefix, $q1,            $q3,             $min_score,
      $max_score,   $verbose,       \@columns_in,    \@columns_out 
     );

    _convert_to_big_bed($bedToBigBedPath, $tmp_file_name, $out_file, $chr_file, $verbose );
    delete_file($tmp_file_name);
  }
  $self->output_param( 'bigbed', \@output_files );
}

sub _check_column_mapping_indices {
    my ($col_mappings) = @_;
    for my $index ( values %$col_mappings ) {
        throw("Column mapping: $index not a valid index")
          if ( defined $index && !( isint($index) && $index > 0 ) );
    }
}

sub _check_output_mapping_possible {
    my ( $in_columns, $out_columns, $add_names ) = @_;

    for my $out_key ( keys %$out_columns ) {
        throw("Cannot find matching input column for output column $out_key")
          unless (
            $in_columns->{$out_key}
            || ( $add_names && $out_key eq 'name' ) # we intend to name features
            || ( $out_key eq 'strand' )
          );                                        # strand can be given as '.'
    }

}

sub _peak_score_interquartile_range {
    my ( $bed_file, $column, $verbose ) = @_;
    my $score_index = $column - 1;

    my $stats = Statistics::Descriptive::Full->new();
    my $fh    = get_open_file_handle($bed_file);
    print STDERR "Reading $bed_file for peak scores$/" if ($verbose);
    while (<$fh>) {
        next if (m/^#/);
        chomp;
        my @vals  = split /\t/;
        my $score = $vals[$score_index];
        $stats->add_data($score);
    }
    close $fh;

    my $q1 = $stats->percentile(25);
    my $q3 = $stats->percentile(75);

    if ($verbose) {
        my $min    = $stats->min();
        my $median = $stats->percentile(50);
        my $max    = $stats->max();
        print STDERR Dumper(
            { peak_score_5_num_summary => [ $min, $q1, $median, $q3, $max ] } );
    }

    return ( $q1, $q3 );
}

sub _column_mappings {
    my ($col_mappings) = @_;
    my @cols;
    map { $cols[ $col_mappings->{$_} - 1 ] = $_ if $col_mappings->{$_} }
      keys %$col_mappings;
    for ( my $i = 0 ; $i < @cols ; $i++ ) {
        $cols[$i] = 'ignore' unless ( defined $cols[$i] );
    }
    return @cols;
}

sub _transform_input_file {
    my (
        $in_file,     $tmp_file_name, $rescale_scores, $add_names,
        $name_prefix, $q1,            $q3,             $min_score,
        $max_score,   $verbose,       $columns_in,     $columns_out
    ) = @_;
    my $in_fh = get_open_file_handle($in_file);
    open( my $out_fh, '>', $tmp_file_name )
      or throw("Could not open $tmp_file_name: $!");

    my $entry_count = 0;

    while (<$in_fh>) {
        if (m/^#/) {
            print $out_fh $_;
            next;
        }
        chomp;
        $entry_count++;

        my @vals = split /\t/;
        my %entry;
        @entry{ @{$columns_in} } = @vals;

        if ( defined $rescale_scores ) {
            if ( $rescale_scores eq 'iqr' ) {
                $entry{score} =
                  rescale_score_iqr( $entry{score}, $q1, $q3, $min_score,
                    $max_score );
            }
            if ( isnum($rescale_scores) ) {
                $entry{score} = $entry{score} * $rescale_scores;
            }
        }
        $entry{score} = int( $entry{score} + 0.5 ) if ( defined $entry{score} );
        if ( !$entry{name} && $add_names ) {
            $entry{name} = $name_prefix . $entry_count;
        }
        if ( !$entry{strand} ) {
            $entry{strand} = '.';
        }
        @vals = @entry{ @{$columns_out} };
        print $out_fh join( "\t", @vals ) . $/;
    }
    close $in_fh;
    close $out_fh;
    print STDERR "Transformed $entry_count rows from $in_file to $tmp_file_name.$/"
      if $verbose;
}

sub rescale_score_iqr {
    my ( $original_score, $q1, $q3, $min_score, $max_score ) = @_;

    if ( $original_score <= $q1 ) {
        return $min_score;
    }
    if ( $original_score >= $q3 ) {
        return $max_score;
    }
    my $range     = $max_score - $min_score;
    my $new_score = $min_score + (
        $range * (
            ( asinh($original_score) - asinh($q1) ) /
              ( asinh($q3) - asinh($q1) )
        )
    );
    return int( $new_score + 0.5 );
}

sub _convert_to_big_bed {
    my ( $bedToBigBedPath, $tmp_file_name, $out_file, $chr_file, $verbose ) = @_;
    print STDERR
      "Converting $tmp_file_name to big bed file $out_file using $chr_file $/"
      if ($verbose);
    execute_system_command("$bedToBigBedPath $tmp_file_name $chr_file $out_file");
}


1;
