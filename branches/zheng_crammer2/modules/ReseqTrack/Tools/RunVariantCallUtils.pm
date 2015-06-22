package ReseqTrack::Tools::RunVariantCallUtils;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

use File::Basename;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	validate_input_string 
);

sub validate_input_string {
	my ($string, $algorithm) = @_;
	my $hash = standard_parameter_hash();
	
	if (!$hash->{$string}) {
		if ($algorithm =~ /samtools/i) {
			throw("parameter has to be one or more of the follows:
					mpileup=>, 
					bcfview=>, 
					vcfutils=>, 
					bcftools_path=>, 
					vcfutils_path=>");
		}	
		elsif ($algorithm =~ /gatk/i ) {
			throw("parameter has to be one or more of the follows:
	  				 		dcov=>
							stand_call_conf=>
							stand_emit_conf=>
							glm=>
							computeSLOD=>
					  		debug_file=>
					  		genotyping_mode=> 
					  		pcr_error_rate=>
					  		p_nonref_model=>
					  		minIndelCnt=>
					  		output_mode=>
					  		min_mapping_quality_score=>
					  		min_base_quality_score=>
					  		metrics_file=>
					  		max_deletion_fraction=>
					  		indel_heterozygosity=>
					  		heterozygosity=>
					  		group=>
					  		annotation=>
					  		alleles=>" );
		}
		elsif ($algorithm =~ /umake/i) {
			throw("parameter has to be one or more of the follows:
							FILTER_MAX_SAMPLE_DP=>, 
							FILTER_MIN_SAMPLE_DP=>, 
							offset_off_target=>, 
							dbSNP=>, 
							hm3_prefix=>,
							indel_prefix=>");
		}	
	}	
	return 1;	
}	

sub standard_parameter_hash{  
  		my %hash;
  		#gatk
  		$hash{"-dcov"} = 1;
        $hash{"-stand_call_conf"} = 1;
        $hash{"-stand_emit_conf"} = 1;
        $hash{"-glm"} = 1;
        $hash{"-computeSLOD"} = 1;
  		$hash{"-debug_file"} = 1;
  		$hash{"-genotyping_mode"} = 1; 
  		$hash{"-pcr_error_rate"} = 1;
  		$hash{"-p_nonref_model"} = 1;
  		$hash{"-minIndelCnt"} = 1;
  		$hash{"-output_mode"} = 1;
  		$hash{"-min_mapping_quality_score"} = 1;
  		$hash{"-min_base_quality_score"} = 1;
  		$hash{"-metrics_file"} = 1;
  		$hash{"-max_deletion_fraction"} = 1;
  		$hash{"-indel_heterozygosity"} = 1;
  		$hash{"-heterozygosity"} = 1;
  		$hash{"-group"} = 1;
  		$hash{"-annotation"} = 1;
  		$hash{"-alleles"} = 1;

		#samtools
  		$hash{"-mpileup"} =1;
  		$hash{"-bcfview"} = 1;
  		$hash{"-vcfutils"} = 1;
  		$hash{"-bcftools_path"} = 1;
  		$hash{"-vcfutils_path"} = 1;
		
		#umake
		$hash{"-FILTER_MIN_SAMPLE_DP"} = 1;
        $hash{"-FILTER_MAX_SAMPLE_DP"} = 1; 
        $hash{"-offset_off_target"} = 1;  
        $hash{"-dbSNP"} = 1; 
		$hash{"-hm3_prefix"} = 1;
		$hash{"-indel_prefix"} = 1; 

  		return \%hash;
  		## FIXME: add more options when other algorithms are implemented 
}  

1;
