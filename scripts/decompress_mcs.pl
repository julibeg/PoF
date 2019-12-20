#! /usr/bin/perl
################################################################################use strict;
################################################################################
# Author:  Christian Jungreuthmayer
# Date:    Mon Mar  2 10:59:05 CET 2015
# Company: Austrian Centre of Industrial Biotechnology (ACIB)
################################################################################

use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_c $opt_p $opt_o $opt_h);
use Data::Dumper;

use constant LINEAR  => 1;
use constant ONEMANY => 2;

my $reac_separator = '';
my ($cfile,$pfile,$ofile);
my $comp_events;
my $all_deco_mcs = [];
my $fho;

# handle command line arguments
read_arguments();

# read in and process change protocol
my $protocol = read_change_protocol($pfile);

# process/decompress cutsets
process_mcs_file($cfile, $protocol);

# open output file
print "INFO: going to open output file '$ofile' for writing\n";
open $fho, ">$ofile" or die "ERROR: couldn't open file '$ofile' for writing: $!\n";
select((select($fho),$|=1)[0]);

# do superset test -> remove all supersets
remove_supersets();

# close output file
close $fho;
print "INFO: decompressed minimal cutsets were written to file '$ofile'\n";
################################################################################


################################################################################
################################################################################
sub remove_supersets
{
   my $hash_ref_arr = [];
   # sort by increase number of knockouts
   print "INFO: remove_supersets(): sorting ...\n";
   @$all_deco_mcs = sort{ scalar(@$a) <=> scalar(@$b) } @$all_deco_mcs;

   print "INFO: remove_supersets(): creating hashes for superset tests ...\n";
   for( my $i = 0; $i < @$all_deco_mcs; $i++ )
   {
      my $hash = {};
      %$hash = map{ $_ => 1 } @{$all_deco_mcs->[$i]};
      push @$hash_ref_arr, $hash;
   }

   print "INFO: remove_supersets(): checking for supersets ...\n";
   for( my $i = 0; $i < @$all_deco_mcs; $i++ )
   {
      print_mcs($all_deco_mcs->[$i]);
      print "INFO: remove_supersets(): processing cutset $i/",scalar(@$all_deco_mcs),"\n";
      for( my $j = $i + 1; $j < @$all_deco_mcs; $j++ )
      {
         my @seen = grep{ exists $hash_ref_arr->[$j]{$_} } keys(%{$hash_ref_arr->[$i]});
         if( scalar(@seen) == scalar(keys(%{$hash_ref_arr->[$i]})) )
         {
            # print "Found superset!!!\n";
            splice @$hash_ref_arr, $j, 1;
            splice @$all_deco_mcs, $j, 1;
            $j--;
         }
      }
   }
}
################################################################################


################################################################################
################################################################################
sub process_mcs_file
{
   my $fn    = shift;
   my $proto = shift;
   my $fh;
   my $num_processed = 0;
   my $tot_num_deco  = 0;

   ############################################################################
   # process file containing compress cutsets
   ############################################################################
   print "INFO: going to open cutset file '$fn'\n";
   open $fh, "$fn" or die "ERROR: couldn't open file '$fn' for reading: $!\n";

   while(<$fh>)
   {
      my $reacs = [];
      my $deco_mcs;
      chomp;
      next if /^\s*#/;
      next if /^\s*$/;
      s/^\s*//;
      s/\s*$//;
      @$reacs = split;
      $deco_mcs = decompress($reacs, $proto);
      my $num_deco = @$deco_mcs;
      $tot_num_deco += $num_deco;

      push @$all_deco_mcs, @$deco_mcs;
      # print_deco_mcs($deco_mcs);

      # print "INFO: number of processed modes: $num_processed\n" if $num_processed%1000 == 0;
      print "INFO: number of processed modes: $num_processed, decompression of last compressed MCS resulted in $num_deco uncompressed MCSs\n";
      $num_processed++;
   }
   close $fh;
   ############################################################################
   
   print "INFO: total number of uncompressed MCSs: $tot_num_deco\n";
}
################################################################################


################################################################################
################################################################################
sub decompress
{
   my $reacs = shift;
   my $proto = shift;
   my $deco_mcs = [];

   # copy compressed cutset to deco_mcs
   @{$deco_mcs->[0]} = @$reacs;

   for( my $i = 0; $i < @$proto; $i++ )
   {
      my $new_deco_mcs = [];
      if( $proto->[$i]->{type} == LINEAR )
      {
         ######################################################################
         # decompress a LINEAR event
         ######################################################################
         for( my $c = 0; $c < @$deco_mcs; $c++ )
         {
            my $found_deco_react = 0;
            my $r;
            for( $r = 0; $r < @{$deco_mcs->[$c]}; $r++ )
            {
               if( $deco_mcs->[$c][$r] eq $proto->[$i]->{new} )
               {
                  $found_deco_react = 1;
                  last;
               }
            }
            if( $found_deco_react )
            {
               # print "LINEAR: orig mode:    @{$deco_mcs->[$c]}\n";
               my $mcs_reacs;
               @$mcs_reacs = @{$deco_mcs->[$c]};
               $deco_mcs->[$c][$r] = $proto->[$i]->{in};
               $mcs_reacs->[$r]    = $proto->[$i]->{out};
               push @$new_deco_mcs, $mcs_reacs;

               # print "LINEAR: changed mode: @{$deco_mcs->[$c]}\n";
               # print "LINEAR: new mode:     @{$mcs_reacs}\n";
            }
         }
         ######################################################################
      }
      elsif( $proto->[$i]->{type} == ONEMANY )
      {
         ######################################################################
         # decompress a ONEMANY event
         ######################################################################
         for( my $c = 0; $c < @$deco_mcs; $c++ )
         {
            my $pos_in_mcs = [];
            my $pos_in_proto_event = [];
            for( my $n = 0; $n < @{$proto->[$i]->{new}}; $n++ )
            {
               for( my $r = 0; $r < @{$deco_mcs->[$c]}; $r++ )
               {
                  if( $deco_mcs->[$c][$r] eq $proto->[$i]->{new}[$n] )
                  {
                     push @$pos_in_mcs, $r;
                     push @$pos_in_proto_event, $n;
                  }
               }
            }

            if( @$pos_in_mcs )
            {
               # print "pos_in_mcs: @$pos_in_mcs\n";
               # print "pos_in_proto_event: @$pos_in_proto_event\n";

               # sort - largest index first, smallest index last
               @$pos_in_mcs = sort{ $b <=> $a } @$pos_in_mcs;
               # print "pos_in_mcs after sorting: @$pos_in_mcs\n";
               # print "cutsets before removing: @{$deco_mcs->[$c]}\n";
               # now remove elements from MCS
               foreach my $r_index (@$pos_in_mcs)
               {
                  splice @{$deco_mcs->[$c]}, $r_index, 1;
               }
               # print "cutsets before removing: @{$deco_mcs->[$c]}\n";
               # print "ONEMANY: orig mode:    @{$deco_mcs->[$c]}\n";
               my $mcs_reacs;
               @$mcs_reacs = @{$deco_mcs->[$c]};

               push @{$deco_mcs->[$c]}, $proto->[$i]->{one}[0];

               foreach my $p_index (@$pos_in_proto_event)
               {
                  my $r_many = $proto->[$i]->{many}[$p_index];
                  push @$mcs_reacs, $r_many;
               }
               push @$new_deco_mcs, $mcs_reacs;

               # print "ONEMANY: changed mode: @{$deco_mcs->[$c]}\n";
               # print "ONEMANY: new mode:     @{$mcs_reacs}\n";
            }
         }
         ######################################################################
      }
      else
      {
         die "ERROR: unknown decompression type '$proto->[$i]->{type}'\n";
      }
      push @$deco_mcs, @$new_deco_mcs;
   }

   # check if there is still a separator character in a reaction name
   # if so, then decompression did not work and exit with error message
   for( my $c = 0; $c < @$deco_mcs; $c++ )
   {
      for( my $r = 0; $r < @{$deco_mcs->[$c]}; $r++ )
      {
         if( $deco_mcs->[$c][$r] =~ /$reac_separator/ )
         {
            die "ERROR: decompressed reaction '$deco_mcs->[$c][$r]' contains separator character '$reac_separator'\n";
         }
      }
   }

   return $deco_mcs;
}
################################################################################


################################################################################
################################################################################
sub print_mcs
{
   my $mcs = shift;

   print $fho "@$mcs\n";
}
################################################################################


################################################################################
################################################################################
sub print_deco_mcs
{
   my $deco_mcs = shift;

   for( my $c = 0; $c < @$deco_mcs; $c++ )
   {
      print $fho "@{$deco_mcs->[$c]}\n";
   }
}
################################################################################


################################################################################
# read in program options
################################################################################
sub read_arguments
{
   getopts('c:p:o:h');

   if( $opt_h )
   {
      usage();
   }

   if( $opt_p )
   {
      $pfile = $opt_p;
   }
   else
   {
      usage('ERROR: name of input file containing change protocol not provided',-1);
   }

   if( $opt_c )
   {
      $cfile = $opt_c;
   }
   else
   {
      usage('ERROR: name of input file containing minimal cutsets not provided ',-1);
   }

   if( $opt_o )
   {
      $ofile = $opt_o;
   }
   else
   {
      usage('ERROR: name of output file decompressed minimal cutsets are written to not provided ',-1);
   }
}
################################################################################


################################################################################
################################################################################
sub read_change_protocol
{
   my $fn = shift;
   my $fh;
   my $proto = [];

   ############################################################################
   # print change protocol to file
   ############################################################################
   print "INFO: going to open change protocol file '$fn'\n";
   open $fh, "$fn" or die "ERROR: couldn't open file '$fn' for reading: $!\n";
   my @chg_proto = <$fh>;
   close $fh;
   ############################################################################


   ############################################################################
   # get reaction separator
   ############################################################################
   my $sep_line = shift @chg_proto;
   chomp $sep_line;

   if( $sep_line =~ /^REACTION_SEPARATOR: (.)$/ )
   {
      $reac_separator = $1;
      print "INFO: reaction separator: '$reac_separator'\n";
   }
   else
   {
      die "ERROR: couldn't get reaction separator from line '$sep_line' from file '$fn'\n";
   }
   ############################################################################


   ############################################################################
   # analyse change protocol
   ############################################################################
   for( my $l = @chg_proto - 1; $l >= 0; $l-- )
   {
      chomp $chg_proto[$l];

      # print "DEBUG: l=$l: change protocol line: '$chg_proto[$l]'\n";

      if( $chg_proto[$l] =~ /^EMPTY_REAC: / )
      {
         # nothing to do -> ignore
         print "INFO: found event 'EMPTY_REAC' -> nothing to do\n";
      }
      elsif( $chg_proto[$l] =~ /^MERGED_IN_IRREV_OUT_REV: / )
      {
         # MERGED_IN_IRREV_OUT_REV: 18 1/1 R_CS 19 -1/1 R_ACONTa=R_ACONTb => 18 R_CS=R_ACONTa=R_ACONTb meta_name=M_cit_c meta_idx=9
         my $event_elem = {};
         ($event_elem->{in}, $event_elem->{out}, $event_elem->{new}) = (split(/\s+/, $chg_proto[$l]))[3,6,9];
         $event_elem->{type} = LINEAR;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_IN_REV_OUT_IRREV: / )
      {
         my $event_elem = {};
         ($event_elem->{in}, $event_elem->{out}, $event_elem->{new}) = (split(/\s+/, $chg_proto[$l]))[3,6,9];
         $event_elem->{type} = LINEAR;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_MANY_IN_ONE_OUT: / )
      {
         my $event_elem = {};
         my $idx;
         $event_elem->{many} = [];
         $event_elem->{one} = [];
         $event_elem->{new} = [];
         # MERGED_MANY_IN_ONE_OUT: i=1: 40 1/1 R_ACKr 59 -1/1 R_ACt2r => 40 R_ACKr=R_ACt2r meta_name=M_ac_c meta_idx=3
         # MERGED_MANY_IN_ONE_OUT: i=0: 34 1/1 R_FEM2 59 -1/1 R_ACt2r => 34 R_FEM2=R_ACt2r meta_name=M_ac_c meta_idx=3
         do
         {
            my ($idxstr, $rmany, $rone, $rnew) = (split(/\s+/, $chg_proto[$l]))[1,4,7,10];
            push @{$event_elem->{many}}, $rmany;
            push @{$event_elem->{one}},  $rone;
            push @{$event_elem->{new}},  $rnew;
            $idx = $idxstr;
            $idx =~ s/i=//;
            $idx =~ s/://;
            $l--;
         }while( $idx > 0 );
         # we went one index to far -> increment it
         $l++;
         $event_elem->{type} = ONEMANY;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_MANY_IN_ONE_REV: / )
      {
         my $event_elem = {};
         my $idx;
         $event_elem->{many} = [];
         $event_elem->{one} = [];
         $event_elem->{new} = [];
         do
         {
            my ($idxstr, $rmany, $rone, $rnew) = (split(/\s+/, $chg_proto[$l]))[1,4,7,10];
            push @{$event_elem->{many}}, $rmany;
            push @{$event_elem->{one}},  $rone;
            push @{$event_elem->{new}},  $rnew;
            $idx = $idxstr;
            $idx =~ s/i=//;
            $idx =~ s/://;
            $l--;
         }while( $idx > 0 );
         # we went one index to far -> increment it
         $l++;
         $event_elem->{type} = ONEMANY;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_MANY_OUT_ONE_REV: / )
      {
         my $event_elem = {};
         my $idx;
         $event_elem->{many} = [];
         $event_elem->{one} = [];
         $event_elem->{new} = [];
         do
         {
            my ($idxstr, $rone, $rmany, $rnew) = (split(/\s+/, $chg_proto[$l]))[1,4,7,10];
            push @{$event_elem->{many}}, $rmany;
            push @{$event_elem->{one}},  $rone;
            push @{$event_elem->{new}},  $rnew;
            $idx = $idxstr;
            $idx =~ s/i=//;
            $idx =~ s/://;
            $l--;
         }while( $idx > 0 );
         # we went one index to far -> increment it
         $l++;
         $event_elem->{type} = ONEMANY;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_ONE_IN_MANY_OUT: / )
      {
         my $event_elem = {};
         my $idx;
         $event_elem->{many} = [];
         $event_elem->{one} = [];
         $event_elem->{new} = [];
         # MERGED_ONE_IN_MANY_OUT: i=1: 47 1/1 R_TRA11 35 -1/1 R_OPM2 => 35 R_TRA11=R_OPM2 meta_name=M_o2_c meta_idx=19
         # MERGED_ONE_IN_MANY_OUT: i=0: 47 1/1 R_TRA11 34 -1/1 R_OPM1 => 34 R_TRA11=R_OPM1 meta_name=M_o2_c meta_idx=19
         do
         {
            my ($idxstr, $rone, $rmany, $rnew) = (split(/\s+/, $chg_proto[$l]))[1,4,7,10];
            push @{$event_elem->{many}}, $rmany;
            push @{$event_elem->{one}},  $rone;
            push @{$event_elem->{new}},  $rnew;
            $idx = $idxstr;
            $idx =~ s/i=//;
            $idx =~ s/://;
            $l--;
         }while( $idx > 0 );
         # we went one index to far -> increment it
         $l++;
         $event_elem->{type} = ONEMANY;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_PROP_REACS: / )
      {
         die "ERROR: decompression of 'MERGED_PROP_REACS' is not supported, yet\n";
      }
      elsif( $chg_proto[$l] =~ /^MERGED_TWO_IRREV: / )
      {
         # MERGED_TWO_IRREV: 64 1/1 R_TRA10 48 -1/1 R_ARA1 => 63 R_TRA10^R_ARA1 meta_name=M_ara_L_c meta_idx=8
         my $event_elem = {};
         ($event_elem->{in}, $event_elem->{out}, $event_elem->{new}) = (split(/\s+/, $chg_proto[$l]))[3,6,9];
         $event_elem->{type} = LINEAR;
         push @$proto, $event_elem;
      }
      elsif( $chg_proto[$l] =~ /^MERGED_TWO_REV: / )
      {
         # MERGED_TWO_REV: 19 1/1 R_ACONTa 20 -1/1 R_ACONTb => 19 R_ACONTa=R_ACONTb meta_name=M_acon_C_c meta_id=5
         my $event_elem = {};
         ($event_elem->{in}, $event_elem->{out}, $event_elem->{new}) = (split(/\s+/, $chg_proto[$l]))[3,6,9];
         $event_elem->{type} = LINEAR;
         push @$proto, $event_elem;
         print "HUGO: found 'MERGED_TWO_REV: $event_elem->{new}'\n"
      }
      elsif( $chg_proto[$l] =~ /^REMOVE_DEADEND_META: / )
      {
         # nothing to do -> ignore
         print "INFO: found event 'REMOVE_DEADEND_META' -> nothing to do\n";
      }
      elsif( $chg_proto[$l] =~ /^REMOVE_UNUSED_META: / )
      {
         # nothing to do -> ignore
         print "INFO: found event 'REMOVE_UNUSED_META' -> nothing to do\n";
      }
      elsif( $chg_proto[$l] =~ /^REVERSED_REVERSIBLE_REACTION: / )
      {
         # nothing to do -> ignore
         print "INFO: found event 'REVERSED_REVERSIBLE_REACTION' -> nothing to do\n";
      }
      else
      {
         die "ERROR: unknown compression event in protocol line '$chg_proto[$l]' in file '$fn'\n";
      }
   }
   ############################################################################

   print_proto_array($proto);
   return $proto;
}
################################################################################


################################################################################
################################################################################
sub print_proto_array
{
   my $proto = shift;

   for( my $e = 0; $e < @$proto; $e++ )
   {
      print "compression event $e:\n";
      my $elem = $proto->[$e];

      if( $elem->{type} == LINEAR )
      {
         print "   type: LINEAR\n";
         print "   in:  $elem->{in}\n";
         print "   out: $elem->{out}\n";
         print "   new: $elem->{new}\n";
      }
      elsif( $elem->{type} == ONEMANY )
      {
         print "   type: ONEMANY\n";
         for( my $i = 0; $i < @{$elem->{one}}; $i++ )
         {
            print "   $i: one:  $elem->{one}[$i]\n";
            print "   $i: many: $elem->{many}[$i]\n";
            print "   $i: new:  $elem->{new}[$i]\n";
         }
      }
      else
      {
         die "ERROR: invalid compression event type '$elem->{type}'\n";
      }
   }

   print "INFO: number of real compression events: ", scalar(@$proto), "\n";
}
################################################################################


################################################################################
################################################################################
sub _whoami
{
   ( caller(1) )[3]
}
################################################################################


################################################################################
################################################################################
sub usage
{
   my $message   = shift || '';
   my $exit_code = shift || 0;

   print "$message\n" if $message;

   print "decompress_mcs.pl -c cfile -p pfile -o ofile [-h]\n";
   print "\n";
   print "-c ..... name of file containing compressed cutsets (input)\n";
   print "-p ..... name of file containing change protocol (input)\n";
   print "-o ..... name of file decompressed cutsets are written to (output)\n";
   print "-h ..... print this message\n";
   print "\n";
   print "decompress_mcs.pl decompresses minimal cutsets of a sytem that was compress with compress_network.pl\n";

   exit($exit_code);
}
################################################################################
