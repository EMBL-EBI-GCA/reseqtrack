#this code is shamelessly stolen from Bio::EnsEMBL::Utils::Argument

=head1 NAME

ReseqTrack::Tools::Argument

=head1 SYNOPSIS

    use ReseqTrack::Tools::Argument qw(rearrange);

    sub new{
      my ($class, @args) = @_;
      my $self = bless {}, $class;

      my ($list) = rearrange([qw(LIST)], @args);

      $self->{'_ID_list'} = [];

      if($list){
	if(ref($list) eq 'ARRAY'){
	  $self->{'_ID_list'} = $list;
	}else{
	  confess("LIST need to be a listref not a $list");
	}
      }
      return $self;
    }


=head1 DESCRIPTION

This is a utility class providing argument parsing methods for object creation

=head1 CONTACT

 laura@ebi.ac.uk or zi@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details exported static class methods. 

=cut


use strict;
use warnings;

package ReseqTrack::Tools::Argument;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(rearrange);



=head2 rearrange

 Usage     : rearrange( array_ref, list_of_arguments)
 Purpose   : Rearranges named parameters to requested order.
 Example   : use Bio::EnsEMBL::Utils::Argument qw(rearrange);
           : rearrange([qw(SEQUENCE ID DESC)],@param);
           : Where @param = (-sequence => $s, 
	         :                 -id       => $i, 
	         :	               -desc     => $d);
 Returns   : @params - an array of parameters in the requested order.
           : The above example would return ($s, $i, $d)
 Argument  : $order : a reference to an array which describes the desired
           :          order of the named parameters.
           : @param : an array of parameters, either as a list (in
           :          which case the function simply returns the list),
           :          or as an associative array with hyphenated tags
           :          (in which case the function sorts the values 
           :          according to @{$order} and returns that new array.)
	   :	      The tags can be upper, lower, or mixed case
           :          but they must start with a hyphen (at least the
           :          first one should be hyphenated.)
 Source    : This function was taken from CGI.pm, written by Dr. Lincoln
           : Stein, and adapted for use in Bio::Seq by Richard Resnick and
           : then adapted for use in Bio::Root::Object.pm by Steve A. Chervitz.
           : This has since been adapted as an exported static method in this 
             class Bio::EnsEMBL::Utils::Argument 
 Comments  : (SAC)
           : This method may not be appropriate for method calls that are
           : within in an inner loop if efficiency is a concern.
           :
           : Parameters can be specified using any of these formats:
           :  @param = (-name=>'me', -color=>'blue');
           :  @param = (-NAME=>'me', -COLOR=>'blue');
           :  @param = (-Name=>'me', -Color=>'blue');
           : A leading hyphenated argument is used by this function to 
           : indicate that named parameters are being used.
           : Therefore, a ('me', 'blue') list will be returned as-is.
           :
           : Note that Perl will confuse unquoted, hyphenated tags as 
           : function calls if there is a function of the same name 
           : in the current namespace:
           :    -name => 'foo' is interpreted as -&name => 'foo'
           :
           : For ultimate safety, put single quotes around the tag:
	         :    ('-name'=>'me', '-color' =>'blue');
           : This can be a bit cumbersome and I find not as readable
           : as using all uppercase, which is also fairly safe:
           :    (-NAME=>'me', -COLOR =>'blue');
	         :
           : Personal note (SAC): I have found all uppercase tags to
           : be more managable: it involves less single-quoting,
           : the code is more readable, and there are no method naming 
           : conlicts.
           : Regardless of the style, it greatly helps to line
	         : the parameters up vertically for long/complex lists.

=cut


sub rearrange {
  my $order = shift;
  $order = shift if($order eq "Bio::EnsEMBL::Utils::Argument"); #skip object if one provided

  # If we've got parameters, we need to check to see whether
  # they are named or simply listed. If they are listed, we
  # can just return them.
  return @_ unless (@_ && $_[0] && substr($_[0], 0,1) eq '-');

  # Convert all of the parameter names to uppercase, and create a
  # hash with parameter names as keys, and parameter values as values
  my $i = 0;
  my (%param) = map {if($i) { $i--; $_; } else { $i++; uc($_); }} @_;
  

  # What we intend to do is loop through the @{$order} variable,
  # and for each value, we use that as a key into our associative
  # array, pushing the value at that key onto our return array.
  return map {$param{uc("-$_")}} @$order;
}
