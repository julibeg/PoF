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
use Math::Fraction;
use Math::MatrixFraction;
use vars qw($opt_s $opt_m $opt_r $opt_n $opt_v $opt_p $opt_o $opt_h $opt_l $opt_k $opt_i);
use Data::Dumper;
use constant NUM_INS        => 0;
use constant NUM_OUTS       => 1;
use constant NUM_REVS       => 2;
use constant POTENTIAL_META => 3;
use constant POTENTIAL_REAC => 4;
use constant IN_IDX         => 5;
use constant OUT_IDX        => 6;
use constant REV_IDX        => 7;

# my @no_onemany_for = qw(R_Ec_biomass_iJO1366_core_53p95M);
# my @no_onemany_for = qw(R_Ec_biomass_iJO1366_core_53p95M R_Ec_biomass_iJO1366_WT_53p95M);
my @no_onemany_for     = ();
my @no_compression_for = ();

my @potential_reac_separartors = ('%', '=', '&', '!', '@');
my $reac_separator;
my $ofile_postfix;
my ($sfile,$mfile,$rfile,$rvfile,$pfile,$nfile,$skipfile);
my $orig_network = {};
my ($num_reacs_orig, $num_metas_orig);
my ($num_reacs_comp, $num_metas_comp);
my @change_protocol;
my $do_only_linear_compression    = 0;
my $skip_kernel_based_compression = 1;
my ($reac_new_names, $reac_new_stoic);

# handle command line arguments
read_arguments();

my $reacs = read_reactions($rfile);
$num_reacs_orig = @$reacs;
print "reacs: @$reacs\n";
%$reac_new_names = map{ $_ => [$_] } @$reacs;
%$reac_new_stoic = map{ $_ => Math::Fraction->new(1,1) } @$reacs;

my $metas = read_metabolites($mfile);
$num_metas_orig = @$metas;
print "metas: @$metas\n";

my $rever = read_reversibility($rvfile);
print "rever: @$rever\n";

if( $nfile )
{
   @no_onemany_for = read_no_one2many_comp($nfile);
   print "no_onemany_for: @no_onemany_for\n";
}

if( $skipfile )
{
   @no_compression_for = read_no_one2many_comp($skipfile);
   print "no_compression_for: @no_compression_for\n";
}

# my $stoim = read_stoichiomat($sfile);
my $stoimDouble = read_stoichiomat($sfile);
my $stoimObj = Math::MatrixFraction->new($stoimDouble);
my $stoim = $stoimObj->get_matrix_fraction();
# print "stoich:\n";
# print_matrix($stoim);
# print_matrix_for_octave($stoimDouble);
# print_matrix_for_perl($stoimDouble);

$orig_network->{stoim} = $stoim;
$orig_network->{reacs} = $reacs;
$orig_network->{metas} = $metas;
$orig_network->{rever} = $rever;


my $compressed_network = analyze_network( $orig_network );
$num_reacs_comp = @{$compressed_network->{reacs}};
$num_metas_comp = @{$compressed_network->{metas}};
print "original network:   number of reactions: $num_reacs_orig, number of metabolites: $num_metas_orig\n";
print "compressed network: number of reactions: $num_reacs_comp, number of metabolites: $num_metas_comp\n";
write_protocol_to_file(@change_protocol);
write_compressed_network( $compressed_network );
write_new_reac_names( $reac_new_names ) if $do_only_linear_compression;
write_new_reac_stoic( $reac_new_stoic ) if $do_only_linear_compression;

################################################################################


################################################################################
################################################################################
sub analyze_network
{
   my $orig_nw = shift;
   my $comp_nws = [];
   my $changes  = 0;

   push @$comp_nws, $orig_nw;

   my $new_nw = remove_empty_reactions( $comp_nws->[-1], \$changes );
   push @$comp_nws, $new_nw;

   do
   {
      $changes = 0;
      my $new_nw1 = remove_irrelevant_metabolites( $comp_nws->[-1], \$changes );
      push @$comp_nws, $new_nw1;

      my $new_nw2 = merge_reactions( $comp_nws->[-1], \$changes );
      push @$comp_nws, $new_nw2;

      if( $skip_kernel_based_compression == 0 )
      {
         my $new_nw3 = kernel_checks( $comp_nws->[-1], \$changes );
         push @$comp_nws, $new_nw3;
      }
   }
   while( $changes );

   return( $comp_nws->[-1] );
}
################################################################################

################################################################################
################################################################################
sub remove_empty_reactions
{
   my $network = shift;
   my $changes = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $new_network;
   my $new_stm = [];
   my $new_rcs = [];
   my $new_mts = [];
   my $new_rvs = [];
   my $keep_reacs = [];

   ############################################################################
   # analyze  metabolites
   ############################################################################
   for( my $r = 0; $r < @$rcs; $r++ )
   {
      my $all_zero = 1;
      for( my $m = 0; $m < @$mts; $m++ )
      {
         if( abs($stm->[$m][$r]) != 0 )
         {
            $all_zero = 0;
            last;
         }
      }

      if( $all_zero == 0 )
      {
         push @$keep_reacs, $r;
      }
      else
      {
         print "EMPTY_REAC: $r $rcs->[$r]\n";
         push @change_protocol, "EMPTY_REAC: $r $rcs->[$r]\n";

         update_new_names($rcs->[$r], 'FOREVER_GONE') if $do_only_linear_compression;
         update_new_stoic($rcs->[$r], 0.0)            if $do_only_linear_compression;

         $$changes++;
      }
   }
   ############################################################################


   ############################################################################
   # create new metabolites array and stoichiometric matrix
   ############################################################################
   for( my $i = 0; $i < @$mts; $i++ )
   {
      $new_mts->[$i] = $mts->[$i];

      my $c_worker = 0;
      for( my $c = 0; $c < @$keep_reacs; $c++ )
      {
         $new_stm->[$i][$c_worker] = $stm->[$i][$keep_reacs->[$c]];
         $c_worker++;
      }
   }
   ############################################################################

   ############################################################################
   ############################################################################
   my $i_worker = 0;
   for( my $i = 0; $i < @$keep_reacs; $i++ )
   {
      $new_rcs->[$i_worker] = $rcs->[$keep_reacs->[$i]];
      $new_rvs->[$i_worker] = $rvs->[$keep_reacs->[$i]];
      $i_worker++;
   }
   ############################################################################

   $new_network->{stoim} = $new_stm;
   $new_network->{reacs} = $new_rcs;
   $new_network->{metas} = $new_mts;
   $new_network->{rever} = $new_rvs;

   
   return $new_network;
}
################################################################################


################################################################################
################################################################################
sub remove_irrelevant_metabolites
{
   my $network = shift;
   my $changes = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $new_network;
   my $new_stm = [];
   my $new_rcs = [];
   my $new_mts = [];
   my $new_rvs = [];
   my $keep_metas = [];
   my $remo_reacs = {};

   ############################################################################
   # analyze  metabolites
   ############################################################################
   for( my $m = 0; $m < @$mts; $m++ )
   {
      my $ins  = 0;
      my $outs = 0;
      my $revs = 0;
      my $in_idx = -1;
      my $out_idx = -1;
      my $rev_idx = -1;
      my $t_remo_reacs = {};
      for( my $r = 0; $r < @$rcs; $r++ )
      {
         if( abs($stm->[$m][$r]) > 0 )
         {
            if( $rvs->[$r] )
            {
               # we are dealing with a reversible reaction
               $revs++;
               $t_remo_reacs->{$r} = 1;
            }
            else
            {
               if( $stm->[$m][$r] < 0 )
               {
                  $outs++;
                  $t_remo_reacs->{$r} = 1;
               }
               else
               {
                  $ins++;
                  $t_remo_reacs->{$r} = 1;
               }
            }
         }
      }

      if( $revs == 0  && $outs == 0 && $ins == 0 )
      {
         print "m=$m ($mts->[$m]): this is an unused metabolite\n";
         push @change_protocol, "REMOVE_UNUSED_META: m=$m, $mts->[$m]\n";
         $$changes++;
      }
      elsif( ($revs == 1 && $ins == 0  && $outs == 0) ||
             ($revs == 0 && $ins >= 1  && $outs == 0) ||
             ($revs == 0 && $ins == 0  && $outs >= 1) )
      {
         print "m=$m ($mts->[$m]), revs=$revs, ins=$ins, outs=$outs: this metabolite is a dead end: ", join(' ',keys(%$t_remo_reacs)), "\n";
         foreach my $tid (keys(%$t_remo_reacs))
         {
            print "    $tid: $rcs->[$tid]\n";
         }
         my @reacs_ids = keys(%$t_remo_reacs);
         my %reacs_names = map{ $rcs->[$_] => $_ } @reacs_ids;
         push @change_protocol, "REMOVE_DEADEND_META: m=$m, $mts->[$m], removed_reactions: ", join(' ',%reacs_names), "\n";

         foreach my $key (keys %$t_remo_reacs)
         {
            $remo_reacs->{$key}++;
         }
         $$changes++;
      }
      else
      {
         # we want to keep this metabolite
         push @$keep_metas, $m;
      }
   }
   ############################################################################


   ############################################################################
   # create new metabolites array and stoichiometric matrix
   ############################################################################
   for( my $i = 0; $i < @$keep_metas; $i++ )
   {
      $new_mts->[$i] = $mts->[$keep_metas->[$i]];

      my $c_worker = 0;
      for( my $c = 0; $c < @{$stm->[0]}; $c++ )
      {
         next if exists $remo_reacs->{$c};
         $new_stm->[$i][$c_worker] = $stm->[$keep_metas->[$i]][$c];
         $c_worker++;
      }
   }
   ############################################################################

   
   my @reacs_removed = map{ $rcs->[$_]} keys %$remo_reacs;
   foreach my $removed_reac (@reacs_removed)
   {
      # print "YYYY: $removed_reac\n";
      update_new_names($removed_reac, 'FOREVER_GONE') if $do_only_linear_compression;
      update_new_stoic($removed_reac, 0.0 )           if $do_only_linear_compression;
   }

   ############################################################################
   ############################################################################
   my $i_worker = 0;
   for( my $i = 0; $i < @$rcs; $i++ )
   {
      next if exists $remo_reacs->{$i};
      $new_rcs->[$i_worker] = $rcs->[$i];
      $new_rvs->[$i_worker] = $rvs->[$i];
      $i_worker++;
   }
   ############################################################################

   $new_network->{stoim} = $new_stm;
   $new_network->{reacs} = $new_rcs;
   $new_network->{metas} = $new_mts;
   $new_network->{rever} = $new_rvs;

   
   return $new_network;
}
################################################################################


################################################################################
################################################################################
sub merge_reactions
{
   my $network = shift;
   my $changes = shift;
   my $t_change;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $new_network = {};
   my $new_stm = [];
   my $new_rcs = [];
   my $new_mts = [];
   my $new_rvs = [];
   my $whoami = _whoami();

   ############################################################################
   # create new metabolites array and stoichiometric matrix
   ############################################################################
   for( my $i = 0; $i < @$stm; $i++ )
   {
      $new_mts->[$i] = $mts->[$i];

      for( my $c = 0; $c < @{$stm->[$i]}; $c++ )
      {
         $new_stm->[$i][$c] = $stm->[$i][$c];
      }
   }
   ############################################################################

   ############################################################################
   ############################################################################
   for( my $i = 0; $i < @$rcs; $i++ )
   {
      $new_rcs->[$i] = $rcs->[$i];
      $new_rvs->[$i] = $rvs->[$i];
   }
   ############################################################################

   $new_network->{stoim} = $new_stm;
   $new_network->{reacs} = $new_rcs;
   $new_network->{metas} = $new_mts;
   $new_network->{rever} = $new_rvs;

   ############################################################################
   # analyze metabolites
   ############################################################################
   do
   {
      $t_change = 0;

      for( my $m = 0; $m < @$new_mts && @$new_mts > 0; $m++ )
      # for( my $m = 0; $m < @$new_mts; $m++ )
      {
         my $ins  = 0;
         my $outs = 0;
         my $revs = 0;
         my $in_idx = [];
         my $out_idx = [];
         my $rev_idx = [];
         for( my $r = 0; $r < @$new_rcs; $r++ )
         {
            if( abs($new_stm->[$m][$r]) > 0  )
            {
               if( $new_rvs->[$r] )
               {
                  # we are dealing with a reversible reaction
                  $revs++;
                  push @$rev_idx, $r;
               }
               else
               {
                  if( $new_stm->[$m][$r] < 0 )
                  {
                     $outs++;
                     push @$out_idx, $r;
                  }
                  else
                  {
                     $ins++;
                     push @$in_idx, $r;
                  }
               }
            }
         }

         # print "revs=$revs outs=$outs ins=$ins\n";
         if( $revs == 0 && $outs == 1 && $ins > 1 && $do_only_linear_compression == 0 &&
             contains_no_bad_reac($in_idx, $out_idx, $new_network) && contains_no_skip_reac($in_idx, $out_idx, $new_network) )
         {
            print "m=$m ($mts->[$m]): 0 reversible, more than 1 input (@$in_idx), 1 output (@$out_idx)\n";
            merge_manyIn_oneOuts_irrev($m, $in_idx, $out_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }
         if( $revs == 0 && $outs >  1 && $ins ==  1 && $do_only_linear_compression == 0 &&
             contains_no_bad_reac($in_idx, $out_idx, $new_network) && contains_no_skip_reac($in_idx, $out_idx, $new_network) )
         {
            print "m=$m ($mts->[$m]): 0 reversible, 1 input (@$in_idx), more than 1 output (@$out_idx)\n";
            merge_oneIn_manyOuts_irrev($m, $in_idx, $out_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }
         if( $revs == 1  && $ins == 0 && $outs > 1 && $do_only_linear_compression == 0 &&
             contains_no_bad_reac($rev_idx, $out_idx, $new_network) && contains_no_skip_reac($rev_idx, $out_idx, $new_network))
         {
            print "m=$m ($mts->[$m]): 1 reversible, 0 input, $outs outputs (@$out_idx)\n";
            merge_manyOuts_one_rev($m, $rev_idx, $out_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }
         if( $revs == 1  && $outs == 0 && $ins > 1 && $do_only_linear_compression == 0 &&
             contains_no_bad_reac($rev_idx, $in_idx, $new_network) && contains_no_skip_reac($rev_idx, $in_idx, $new_network))
         {
            print "m=$m ($mts->[$m]): 1 reversible (@$rev_idx), 0 output, $ins inputs (@$in_idx)\n";
            merge_manyIn_one_rrev($m, $in_idx, $rev_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }

         if( $revs == 1  && $outs == 0 && $ins == 1 && contains_no_skip_reac($rev_idx, $in_idx, $new_network) )
         {
            print "m=$m ($mts->[$m]): 1 reversible (@$rev_idx), 0 output, $ins inputs (@$in_idx)\n";
            merge_manyIn_one_rrev($m, $in_idx, $rev_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }
         if( $revs == 1  && $ins == 0 && $outs == 1 && contains_no_skip_reac($rev_idx, $out_idx, $new_network) )
         {
            merge_manyOuts_one_rev($m, $rev_idx, $out_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }

         if( $revs == 2  && $outs == 0 && $ins == 0 && contains_no_skip_reac($rev_idx, $rev_idx, $new_network) )
         {
            print "MERRG_TWO_REV: m=$m ($new_mts->[$m]): 2 reversible (@$rev_idx), 0 input, 0 outputs\n";
            merge_two_rev_reacs($m, $rev_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }

         if( $revs == 0  && $outs == 1 && $ins == 1 && contains_no_skip_reac($in_idx, $out_idx, $new_network) )
         {
            print "MERGE_TWO_IRREV: m=$m ($new_mts->[$m]): 0 reversible, 1 input (@$in_idx), 1 output (@$out_idx)\n";
            merge_two_irrev_reacs($m, $in_idx, $out_idx, $new_network);
            $t_change++;
            $$changes++;
            # last;
            next;
         }
      }
   }
   while( $t_change);
   ############################################################################

   return $new_network;
}
################################################################################


################################################################################
################################################################################
sub contains_no_skip_reac
{
   my $ins     = shift;
   my $outs    = shift;
   my $network = shift;

   my $rcs = $network->{reacs};
   print "DEBUG: contains_no_skip_reac(): ins=@$ins outs=@$outs\n";

   for( my $i = 0; $i < @$ins; $i++ )
   {
      my $comp_reacs = $reac_separator . $rcs->[$ins->[$i]] . $reac_separator;
      for( my $j = 0; $j < @no_compression_for; $j++ )
      {
         my $this_reac = $reac_separator . $no_compression_for[$j] . $reac_separator;
         print "INFO: contains_no_bad_reac(): ins: comp_reacs=$comp_reacs this_reac=$this_reac\n";
         if( $comp_reacs =~ /$this_reac/ )
         {
            print "INFO: contains_no_bad_reac: found bad reaction '$no_compression_for[$j]' -> do not compress\n";
            return 0;
         }
      }
   }

   for( my $i = 0; $i < @$outs; $i++ )
   {
      my $comp_reacs = $reac_separator . $rcs->[$outs->[$i]] . $reac_separator;
      for( my $j = 0; $j < @no_compression_for; $j++ )
      {
         my $this_reac = $reac_separator . $no_compression_for[$j] . $reac_separator;
         print "INFO: contains_no_bad_reac(): outs: comp_reacs=$comp_reacs this_reac=$this_reac\n";
         if( $comp_reacs =~ /$this_reac/ )
         {
            print "INFO: contains_no_bad_reac: found bad reaction '$no_compression_for[$j]' -> do not compress\n";
            return 0;
         }
      }
   }

   return 1;
}
################################################################################


################################################################################
################################################################################
sub contains_no_bad_reac
{
   my $ins     = shift;
   my $outs    = shift;
   my $network = shift;

   my $rcs = $network->{reacs};

   for( my $i = 0; $i < @$ins; $i++ )
   {
      my $comp_reacs = $reac_separator . $rcs->[$ins->[$i]] . $reac_separator;
      for( my $j = 0; $j < @no_onemany_for; $j++ )
      {
         my $this_reac = $reac_separator . $no_onemany_for[$j] . $reac_separator;
         print "INFO: contains_no_bad_reac(): ins: comp_reacs=$comp_reacs this_reac=$this_reac\n";
         if( $comp_reacs =~ /$this_reac/ )
         {
            print "INFO: contains_no_bad_reac: found bad reaction '$no_onemany_for[$j]' -> do not compress\n";
            return 0;
         }
      }
   }

   for( my $i = 0; $i < @$outs; $i++ )
   {
      my $comp_reacs = $reac_separator . $rcs->[$outs->[$i]] . $reac_separator;
      for( my $j = 0; $j < @no_onemany_for; $j++ )
      {
         my $this_reac = $reac_separator . $no_onemany_for[$j] . $reac_separator;
         print "INFO: contains_no_bad_reac(): outs: comp_reacs=$comp_reacs this_reac=$this_reac\n";
         if( $comp_reacs =~ /$this_reac/ )
         {
            print "INFO: contains_no_bad_reac: found bad reaction '$no_onemany_for[$j]' -> do not compress\n";
            return 0;
         }
      }
   }

   return 1;
}
################################################################################


################################################################################
################################################################################
sub kernel_checks
{
   my $network = shift;
   my $changes = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $new_network = {};
   my $new_stm = [];
   my $new_rcs = [];
   my $new_mts = [];
   my $new_rvs = [];
   my $whoami = _whoami();


   ############################################################################
   # create new metabolites array and stoichiometric matrix
   ############################################################################
   for( my $i = 0; $i < @$stm; $i++ )
   {
      $new_mts->[$i] = $mts->[$i];

      for( my $c = 0; $c < @{$stm->[$i]}; $c++ )
      {
         # $new_stm->[$i][$c] = Math::Fraction->new($stm->[$i][$c], $Math::Fraction::REDUCE);
         $new_stm->[$i][$c] = Math::Fraction->new($stm->[$i][$c]);
         # $new_stm->[$i][$c] = $stm->[$i][$c];
      }
   }
   ############################################################################

   ############################################################################
   ############################################################################
   for( my $i = 0; $i < @$rcs; $i++ )
   {
      $new_rcs->[$i] = $rcs->[$i];
      $new_rvs->[$i] = $rvs->[$i];
   }
   ############################################################################

   $new_network->{stoim} = $new_stm;
   $new_network->{reacs} = $new_rcs;
   $new_network->{metas} = $new_mts;
   $new_network->{rever} = $new_rvs;
   ############################################################################

   # we do not need to create the kernel if the stoichiometric matrix
   # is already reduced to a system without metabolites!
   return $new_network unless @$stm;

   ############################################################################
   ############################################################################
   print "INFO: $whoami: going to create Math::MatrixFraction object ...\n";
   my $matrixObj = Math::MatrixFraction->new_from_fraction($new_stm);

   print "INFO: $whoami: going to compute kernel of stoichiometric matrix ...\n";
   my $kernel = $matrixObj->compute_kernel();
   print "INFO: $whoami: kernel:\n";
   print_matrix($kernel);
   # debugging: begin
   # my $kernelObj = Math::MatrixFraction->new_from_fraction($kernel);
   # my $resMultiObj = $matrixObj->multiply($kernel);
   # print "stoichiometric matrix x kernel (fraction):\n";
   # $resMultiObj->print_matrix_fraction();
   # debugging: end
   ############################################################################

   ############################################################################
   ############################################################################
   my $num_reacs_removed = 0;
   my $r_worker = 0;
   for( my $r = 0; $r < @$kernel; $r++ )
   {
      my $all_zero = 1;
      for( my $c = 0; $c < @{$kernel->[$r]}; $c++ )
      {
         if( $kernel->[$r][$c] != 0 )
         {
            $all_zero = 0;
            last;
         }
      }

      if( $all_zero )
      {
         print "Reaction $new_rcs->[$r_worker] (worker_id=$r_worker) has a zero-row in kernel -> no flux will ever be carried -> can be removed\n";
         # remove reaction from reversibility array
         splice @$new_rvs, $r_worker, 1;
         # remove reaction from reaction name array
         splice @$new_rcs, $r_worker, 1;
         # remove reaction from stoichiometric matrix
         extract_reac( $new_stm, $r_worker);

         $num_reacs_removed++;
         $r_worker--;
         $$changes++;
      }
      $r_worker++;
   }
   print "INFO: $whoami: removed $num_reacs_removed reactions from network\n";
   ############################################################################


   ############################################################################
   # get proportional set of reactions
   ############################################################################
   if( @$kernel )
   {
      my ($prop_sets, $prop_facs) = get_prop_sets($kernel);

      for( my $i = 0; $i < @$prop_sets; $i++ )
      {
         print "DEBUG; $whoami: set i=$i:";
         for( my $j = 0; $j < @{$prop_sets->[$i]}; $j++ )
         {
            print " $prop_sets->[$i][$j] ($new_rcs->[$prop_sets->[$i][$j]]): prop_fac: $prop_facs->[$i][$j]";
         }
         print "\n";
      }

      if( @$prop_sets )
      {
         # merge reactions which are proportional in nullspace
         merge_coupled_reacs($new_network, $prop_sets, $prop_facs);
         $$changes++;
      }
   }
   ############################################################################

   return $new_network;
}
################################################################################


################################################################################
################################################################################
sub merge_coupled_reacs
{
   my $network   = shift;
   my $prop_sets = shift;
   my $prop_facs = shift;
   my $whoami    = _whoami();

   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};

   my $to_extract;

   # copy reaction name for output in change protocol
   my $rcs_save;
   for( my $i = 0; $i < @$rcs; $i++ )
   {
      $rcs_save->[$i] = $rcs->[$i];
   }

   my $revers_hash;
   my $new_name;
   my $new_ids;
   for( my $i = 0; $i < @$prop_sets; $i++ )
   {
      print "DEBUG; $whoami: set i=$i:";
      for( my $j = 0; $j < @{$prop_sets->[$i]}; $j++ )
      {
         print " $prop_sets->[$i][$j] ($rcs->[$prop_sets->[$i][$j]])";
      }
      print "\n";

      my $dest_id     = $prop_sets->[$i][0];
      $new_ids->[$i]  = $dest_id;
      $new_name->[$i] = $rcs->[$dest_id]; 
      $revers_hash->{$dest_id} = $rvs->[$dest_id];
      for(  my $j = 1; $j < @{$prop_sets->[$i]}; $j++ )
      {
         my $src_id = $prop_sets->[$i][$j];
         my $factor = $prop_facs->[$i][$j];

         $revers_hash->{$dest_id} = $revers_hash->{$dest_id} && $rvs->[$src_id];
         $new_name->[$i] .= $reac_separator . $rcs->[$src_id];
         push @$to_extract, $src_id;

         # multiply_matrix_row_by( $stm, $src_id, $factor );
         divide_matrix_row_by( $stm, $src_id, $factor );
         add_col2_to_col1_of_matrix( $stm, $dest_id, $src_id );
      }

      print "DEBUG: $whoami: revers_hash->{$dest_id}=$revers_hash->{$dest_id}\n";
      # set new reversibility of reaction
      $rvs->[$dest_id] = $revers_hash->{$dest_id};
      $rcs->[$dest_id] = $new_name->[$i];
   }

   @$to_extract = sort { $b <=> $a } @$to_extract;
   print "DEBUG: $whoami: to_extract: @$to_extract\n";

   for( my $i = 0; $i < @$to_extract; $i++ )
   {
      my $id = $to_extract->[$i];
      print "DEBUG: $whoami: i=$i: extracting reaction $id: $rcs->[$id]\n";
      # update stoichiometric matrix
      extract_reac( $stm, $id );

      # update array of reaction names
      # splice( $rcs, $id, 1 );
      splice( @$rcs, $id, 1 );

      # update array of reversibilities
      # splice( $rvs, $id, 1 );
      splice( @$rvs, $id, 1 );

      for( my $j = 0; $j < @$new_ids; $j++ )
      {
         $new_ids->[$j]-- if $id < $prop_sets->[$j][0];
      }
   }

   # create new ids
   

   for( my $i = 0; $i < @$prop_sets; $i++ )
   {
      my $dest_id  = $prop_sets->[$i][0];
      my $chg_str = "MERGED_PROP_REACS: $i: ";
      for(  my $j = 1; $j < @{$prop_sets->[$i]}; $j++ )
      {
         my $src_id = $prop_sets->[$i][$j];
         my $factor = $prop_facs->[$i][$j];
         
         $chg_str .= " $src_id $rcs_save->[$src_id] $factor";
      }
      $chg_str .= ", $dest_id $rcs_save->[$dest_id] => ";
      $chg_str .= " $new_ids->[$i] $new_name->[$i]\n";
      push @change_protocol, $chg_str;
   }
}
################################################################################


################################################################################
################################################################################
sub get_prop_sets
{
   my $kernel = shift;
   my $whoami = _whoami();

   my $kernel_idx;

   for( my $r = 0; $r < @$kernel; $r++ )
   {
      my $row;
      for( my $c = 0; $c < @{$kernel->[$r]}; $c++ )
      {
         $row->[$c] = Math::Fraction->new($kernel->[$r][$c]);
      }
      $kernel_idx->[$r]->[0] = $r;
      $kernel_idx->[$r]->[1] = $row;
   }

   @$kernel_idx = sort kernel_sorter @$kernel_idx;
   # print "INFO: $whoami: kernel (sorted):\n";
   # print_matrix_idx($kernel_idx);
   my $proportionals = [];
   my $proportio_fac = [];
   my $num_proportionals = -1;
   my $last_two_were_proportional = 0;
   my $idx_1 = 0;
   my $idx_2 = 1;
   do
   {
      my $row_1 = get_row( $kernel_idx, $idx_1);
      my $row_2 = get_row( $kernel_idx, $idx_2);

      if( my $prop_fac = are_proportional($row_1, $row_2) )
      {
         if( $last_two_were_proportional == 0 )
         {
            $num_proportionals++;
            # print "num_proportionals=$num_proportionals\n";
            # print "idx_1=$idx_1, index_1=",$kernel_idx->[$idx_1]->[0],"\n";
            # print "idx_2=$idx_2, index_2=",$kernel_idx->[$idx_2]->[0],"\n";
            push @{$proportionals->[$num_proportionals]}, $kernel_idx->[$idx_1]->[0];
            push @{$proportio_fac->[$num_proportionals]}, Math::Fraction->new(1,1);
            push @{$proportionals->[$num_proportionals]}, $kernel_idx->[$idx_2]->[0];
            push @{$proportio_fac->[$num_proportionals]}, $prop_fac;
         }
         else
         {
            # print "num_proportionals=$num_proportionals\n";
            # print "idx_n=$idx_2, index_n=",$kernel_idx->[$idx_2]->[0],"\n";
            push @{$proportionals->[$num_proportionals]}, $kernel_idx->[$idx_2]->[0];
            push @{$proportio_fac->[$num_proportionals]}, $prop_fac;
         }
         $last_two_were_proportional = 1;
      }
      else
      {
         $last_two_were_proportional = 0;
      }
      $idx_1++;
      $idx_2++;
   } while( $idx_2 < @$kernel_idx );

   print "INFO: we found ", $num_proportionals + 1, " sets of proportional reactions\n";

   # for( my $i = 0; $i < @$proportionals; $i++ )
   # {
   #    print "DEBUG: $whoami: set i=$i:";
   #    for( my $j = 0; $j < @{$proportionals->[$i]}; $j++ )
   #    {
   #       print " $proportionals->[$i][$j]";
   #    }
   #    print "\n";
   # }

   return $proportionals, $proportio_fac;
}
################################################################################


################################################################################
################################################################################
sub get_row
{
   my $mat = shift;
   my $idx = shift;
   my $row;
   # my $sum = Math::Fraction->new(0,1);

   my $whoami = _whoami();

   print "DEBUG: $whoami: getting row $idx\n";

   for( my $c = 0; $c < @{$mat->[$idx]->[1]}; $c++ )
   {
      $row->[$c] = Math::Fraction->new($mat->[$idx]->[1]->[$c]);
   }

   print "DEBUG: $whoami: row: @$row\n";

   return $row;
}
################################################################################


################################################################################
################################################################################
sub are_proportional
{
   my $vec1 = shift;
   my $vec2 = shift;
   my $g_prop = undef;

   my $whoami = _whoami();

   my $len1 = @$vec1;
   my $len2 = @$vec2;

   print "DEBUG: $whoami: comparing the following rows:\n";
   print "DEBUG: $whoami: row1: @$vec1\n";
   print "DEBUG: $whoami: row2: @$vec2\n";

   die "ERROR: $whoami: length of vector #1 ($len1) is not equalt to length of vector #2 ($len2)" unless $len1 == $len2;

   for( my $c = 0; $c < $len1; $c++ )
   {
      next if $vec1->[$c] == 0 && $vec2->[$c] == 0;

      return 0 if $vec1->[$c] == 0 && $vec2->[$c] != 0;
      return 0 if $vec1->[$c] != 0 && $vec2->[$c] == 0;

      my $t_prop = $vec1->[$c]/$vec2->[$c];

      if( defined $g_prop )
      {
         return 0 unless $t_prop == $g_prop;
      }
      else
      {
         $g_prop = $t_prop;
      }
   }

   print "INFO: $whoami: we found two proportional rows:\n";
   print "INFO: $whoami: row #1: @$vec1\n";
   print "INFO: $whoami: row #2: @$vec2\n";

   return $g_prop;
}
################################################################################


################################################################################
################################################################################
sub kernel_sorter
{
   for( my $i = 0; $i < @{$a->[1]}; $i++ )
   {
      if( $a->[1]->[$i] != 0 && $b->[1]->[$i] == 0 )
      {
         return 1;
      }
      if( $a->[1]->[$i] == 0 && $b->[1]->[$i] != 0 )
      {
         return -1;
      }
   }
   return 0;
}
################################################################################


################################################################################
################################################################################
sub merge_two_rev_reacs
{
   my $m_idx = shift;
   my $rev_idx = shift;
   my $network = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $fac1;
   my $fac2;

   my $whoami = _whoami();

   print "reac1: "; print_reaction($network, $rev_idx->[0]);
   print "reac2: "; print_reaction($network, $rev_idx->[1]);

   my $reac_name_1 = $rcs->[$rev_idx->[0]];
   my $reac_name_2 = $rcs->[$rev_idx->[1]];

   my $stoim_val_1 = $stm->[$m_idx][$rev_idx->[0]];
   my $stoim_val_2 = $stm->[$m_idx][$rev_idx->[1]];
   die "ERROR: $whoami: stoim_val_1 is equal to zero\n" if $stoim_val_1 == 0;
   die "ERROR: $whoami: stoim_val_2 is equal to zero\n" if $stoim_val_2 == 0;

   if( $stoim_val_1 > 0.0 && $stoim_val_2 > 0.0 )
   {
      $fac1 =  1.0;
      $fac2 = -1.0;
   }
   elsif( $stoim_val_1 < 0.0 && $stoim_val_2 < 0.0 )
   {
      $fac1 = -1.0;
      $fac2 =  1.0;
   }
   elsif( $stoim_val_1 > 0.0 && $stoim_val_2 < 0.0 )
   {
      $fac1 =  1.0;
      $fac2 =  1.0;
   }
   elsif( $stoim_val_1 < 0.0 && $stoim_val_2 > 0.0 )
   {
      $fac1 = -1.0;
      $fac2 = -1.0;
   }

   divide_matrix_row_by( $stm, $rev_idx->[0], $stoim_val_1);
   divide_matrix_row_by( $stm, $rev_idx->[1], $stoim_val_2);

   # negate fac2, as reac2 is subtracted from reac1
   sub_col2_from_col1_of_matrix( $stm, $rev_idx->[0], $rev_idx->[1] );

   # set new reaction name
   my $reac_name_new = $reac_name_1 . $reac_separator . $reac_name_2;
   $rcs->[$rev_idx->[0]] = $reac_name_new;

   my $reac_name_extract = splice @$rcs, $rev_idx->[1], 1;

   my $reac_extract      = extract_reac($stm, $rev_idx->[1]);

   # remove reversibility entry from array
   # splice @$rvs, $rev_idx->[1], 1;
   splice @$rvs, $rev_idx->[1], 1;

   # first check if really all entries for metabolite are zero
   print_metabolite($network, $m_idx);
   for( my $r =  0; $r < @{$stm->[$m_idx]}; $r++ )
   {
      die "ERROR: $whoami: r=$r element for metabolite $m_idx is not zero. This should never happen!\n" if $stm->[$m_idx][$r] != 0;
   }

   # remove involved metabolite from name array
   my $meta_name = splice @$mts, $m_idx, 1;

   # remove involved metabolite from stoichiometric matrix
   splice @$stm, $m_idx, 1;

   my $new_idx = $rev_idx->[0];
   $new_idx-- if $rev_idx->[1] < $rev_idx->[0];

   push @change_protocol, "MERGED_TWO_REV: $rev_idx->[0] $stoim_val_1 $reac_name_1 $rev_idx->[1] $stoim_val_2 $reac_name_2 => $new_idx $reac_name_new meta_name=$meta_name meta_id=$m_idx\n";
   update_new_names($reac_name_new, $reac_name_new)    if $do_only_linear_compression;
   update_new_stoic($reac_name_1, $fac1/$stoim_val_1)  if $do_only_linear_compression;
   update_new_stoic($reac_name_2, $fac2/$stoim_val_2)  if $do_only_linear_compression;

   print "merged: "; print_reaction($network, $new_idx);
}
################################################################################


################################################################################
################################################################################
sub merge_manyIn_one_rrev
{
   my $m_idx = shift;
   my $ins_idx = shift;
   my $out_idx = shift;
   my $network = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $fac1 =  1.0;
   my $fac2 = -1.0;

   my $whoami = _whoami();

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      print "i=$i: reac1: "; print_reaction($network, $ins_idx->[$i]);
   }
   print "reac2: "; print_reaction($network, $out_idx->[0]);

   my $reac_name_out = $rcs->[$out_idx->[0]];
   my $reac_name_in  = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      push @$reac_name_in, $rcs->[$ins_idx->[$i]];
   }

   my $stoim_val_out = $stm->[$m_idx][$out_idx->[0]];
   my $stoim_val_in  = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      push @$stoim_val_in, $stm->[$m_idx][$ins_idx->[$i]];
   }
   die "ERROR: $whoami: stoim_val_out is equal to zero\n"  if $stoim_val_out == 0;
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      die "ERROR: $whoami: i=$i: stoim_val_in is equal to zero\n" if $stoim_val_in->[$i] == 0;
   }

   if( $stoim_val_out > 0 )
   {
      # we are dealing with a reversible reaction that has a positive contribution
      # hence, we reverse the direction by multiplying/dividing the column by -1
      divide_matrix_row_by( $stm, $out_idx->[0], Math::Fraction->new(-1,1));
      $fac2 *= -1;
      push @change_protocol, "REVERSED_REVERSIBLE_REACTION: $out_idx->[0] $stoim_val_out $reac_name_out\n";
   }

   divide_matrix_row_by( $stm, $out_idx->[0], abs($stoim_val_out));
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      divide_matrix_row_by( $stm, $ins_idx->[$i], abs($stoim_val_in->[$i]));
   }

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      add_col2_to_col1_of_matrix( $stm, $ins_idx->[$i], $out_idx->[0] );
   }

   # set new reaction name
   my $reac_name_new = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      $reac_name_new->[$i]   = $reac_name_in->[$i] . $reac_separator . $reac_name_out;
      $rcs->[$ins_idx->[$i]] = $reac_name_new->[$i];
   }

   my $reac_name_extract = splice @$rcs, $out_idx->[0], 1;

   my $reac_extract      = extract_reac($stm, $out_idx->[0]);

   # remove reversibility entry from array
   splice @$rvs, $out_idx->[0], 1;

   print_metabolite($network, $m_idx);
   for( my $r =  0; $r < @{$stm->[$m_idx]}; $r++ )
   {
      die "ERROR: $whoami: r=$r element for metabolite $m_idx is not zero. This should never happen!\n" if $stm->[$m_idx][$r] != 0;
   }

   # remove involved metabolite from name array
   my $meta_name = splice @$mts, $m_idx, 1;

   # remove involved metabolite from stoichiometric matrix
   splice @$stm, $m_idx, 1;

   my $new_idx = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      $new_idx->[$i] = $ins_idx->[$i];
      $new_idx->[$i]-- if $out_idx->[0] < $ins_idx->[$i];
   }

   if( @{$ins_idx} == 1 )
   {
      push @change_protocol, "MERGED_IN_IRREV_OUT_REV: $ins_idx->[0] $stoim_val_in->[0] $reac_name_in->[0] $out_idx->[0] $stoim_val_out $reac_name_out => $new_idx->[0] $reac_name_new->[0] meta_name=$meta_name meta_idx=$m_idx\n";
      update_new_names($reac_name_new->[0], $reac_name_new->[0])     if $do_only_linear_compression;
      update_new_stoic($reac_name_in->[0], $fac1/$stoim_val_in->[0]) if $do_only_linear_compression;
      update_new_stoic($reac_name_out,     $fac2/$stoim_val_out)     if $do_only_linear_compression;
   }
   else
   {
      for( my $i = 0; $i < @{$ins_idx}; $i++ )
      {
         push @change_protocol, "MERGED_MANY_IN_ONE_REV: i=$i: $ins_idx->[$i] $stoim_val_in->[$i] $reac_name_in->[$i] $out_idx->[0] $stoim_val_out $reac_name_out => $new_idx->[$i] $reac_name_new->[$i] meta_name=$meta_name meta_idx=$m_idx\n";
      }
   }

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      print "i=$i: merged: "; print_reaction($network, $new_idx->[$i]);
   }
}
################################################################################


################################################################################
################################################################################
sub merge_manyIn_oneOuts_irrev
{
   my $m_idx = shift;
   my $ins_idx = shift;
   my $out_idx = shift;
   my $network = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};

   my $whoami = _whoami();

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      print "i=$i: reac1: "; print_reaction($network, $ins_idx->[$i]);
   }
   print "reac2: "; print_reaction($network, $out_idx->[0]);

   my $reac_name_out = $rcs->[$out_idx->[0]];
   my $reac_name_in  = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      push @$reac_name_in, $rcs->[$ins_idx->[$i]];
   }

   my $stoim_val_out = $stm->[$m_idx][$out_idx->[0]];
   my $stoim_val_in  = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      push @$stoim_val_in, $stm->[$m_idx][$ins_idx->[$i]];
   }
   die "ERROR: $whoami: stoim_val_out is equal to zero\n"  if $stoim_val_out == 0;
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      die "ERROR: $whoami: i=$i: stoim_val_in is equal to zero\n" if $stoim_val_in->[$i] == 0;
   }

   divide_matrix_row_by( $stm, $out_idx->[0], abs($stoim_val_out));
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      divide_matrix_row_by( $stm, $ins_idx->[$i], abs($stoim_val_in->[$i]));
   }

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      add_col2_to_col1_of_matrix( $stm, $ins_idx->[$i], $out_idx->[0] );
   }

   # set new reaction name
   my $reac_name_new = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      $reac_name_new->[$i]   = $reac_name_in->[$i] . $reac_separator . $reac_name_out;
      $rcs->[$ins_idx->[$i]] = $reac_name_new->[$i];
   }

   my $reac_name_extract = splice @$rcs, $out_idx->[0], 1;

   my $reac_extract      = extract_reac($stm, $out_idx->[0]);

   # remove reversibility entry from array
   splice @$rvs, $out_idx->[0], 1;

   print_metabolite($network, $m_idx);
   for( my $r =  0; $r < @{$stm->[$m_idx]}; $r++ )
   {
      die "ERROR: $whoami: r=$r element for metabolite $m_idx is not zero. This should never happen!\n" if $stm->[$m_idx][$r] != 0;
   }

   # remove involved metabolite from name array
   my $meta_name = splice @$mts, $m_idx, 1;

   # remove involved metabolite from stoichiometric matrix
   splice @$stm, $m_idx, 1;

   my $new_idx = [];
   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      $new_idx->[$i] = $ins_idx->[$i];
      $new_idx->[$i]-- if $out_idx->[0] < $ins_idx->[$i];
   }

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      push @change_protocol, "MERGED_MANY_IN_ONE_OUT: i=$i: $ins_idx->[$i] $stoim_val_in->[$i] $reac_name_in->[$i] $out_idx->[0] $stoim_val_out $reac_name_out => $new_idx->[$i] $reac_name_new->[$i] meta_name=$meta_name meta_idx=$m_idx\n";
   }

   for( my $i = 0; $i < @{$ins_idx}; $i++ )
   {
      print "i=$i: merged: "; print_reaction($network, $new_idx->[$i]);
   }
}
################################################################################


################################################################################
################################################################################
sub merge_manyOuts_one_rev
{
   my $m_idx = shift;
   my $ins_idx = shift;
   my $out_idx = shift;
   my $network = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $fac1 = 1.0;
   my $fac2 = -1.0;

   my $whoami = _whoami();

   print "reac1: "; print_reaction($network, $ins_idx->[0]);
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      print "i=$i: reac2: "; print_reaction($network, $out_idx->[$i]);
   }

   my $reac_name_in  = $rcs->[$ins_idx->[0]];
   my $reac_name_out = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      push @$reac_name_out, $rcs->[$out_idx->[$i]];
   }

   my $stoim_val_in  = $stm->[$m_idx][$ins_idx->[0]];
   my $stoim_val_out = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      push @$stoim_val_out, $stm->[$m_idx][$out_idx->[$i]];
   }
   die "ERROR: $whoami: stoim_val_in is equal to zero\n"  if $stoim_val_in  == 0;
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      die "ERROR: $whoami: i=$i: stoim_val_out is equal to zero\n" if $stoim_val_out->[$i] == 0;
   }

   if( $stoim_val_in < 0 )
   {
      # we are dealing with a reversible reaction that has a negative contribution
      # hence, we reverse the direction by multiplying/dividing the column by -1
      divide_matrix_row_by( $stm, $ins_idx->[0], Math::Fraction->new(-1,1));
      $fac1 *= -1;
      push @change_protocol, "REVERSED_REVERSIBLE_REACTION: $ins_idx->[0] $stoim_val_in $reac_name_in\n";
   }

   divide_matrix_row_by( $stm, $ins_idx->[0], abs($stoim_val_in));
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      divide_matrix_row_by( $stm, $out_idx->[$i], abs($stoim_val_out->[$i]));
   }

   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      add_col2_to_col1_of_matrix( $stm, $out_idx->[$i], $ins_idx->[0] );
   }

   # set new reaction name
   my $reac_name_new = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      $reac_name_new->[$i] = $reac_name_in . $reac_separator . $reac_name_out->[$i];
      $rcs->[$out_idx->[$i]] = $reac_name_new->[$i];
   }

   my $reac_name_extract = splice @$rcs, $ins_idx->[0], 1;

   my $reac_extract      = extract_reac($stm, $ins_idx->[0]);

   # remove reversibility entry from array
   splice @$rvs, $ins_idx->[0], 1;

   print_metabolite($network, $m_idx);
   for( my $r =  0; $r < @{$stm->[$m_idx]}; $r++ )
   {
      die "ERROR: $whoami: r=$r element for metabolite $m_idx is not zero. This should never happen!\n" if $stm->[$m_idx][$r] != 0;
   }

   # remove involved metabolite from name array
   my $meta_name = splice @$mts, $m_idx, 1;

   # remove involved metabolite from stoichiometric matrix
   splice @$stm, $m_idx, 1;

   my $new_idx = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      $new_idx->[$i] = $out_idx->[$i];
      $new_idx->[$i]-- if $ins_idx->[0] < $out_idx->[$i];
   }

   if( @{$out_idx} == 1 )
   {
      push @change_protocol, "MERGED_IN_REV_OUT_IRREV: $ins_idx->[0] $stoim_val_in $reac_name_in $out_idx->[0] $stoim_val_out->[0] $reac_name_out->[0] => $new_idx->[0] $reac_name_new->[0] meta_name=$meta_name meta_idx=$m_idx\n";
      update_new_names($reac_name_new->[0], $reac_name_new->[0])       if $do_only_linear_compression;
      update_new_stoic($reac_name_in, $fac1/$stoim_val_in)             if $do_only_linear_compression;
      update_new_stoic($reac_name_out->[0], $fac2/$stoim_val_out->[0]) if $do_only_linear_compression;
   }
   else
   {
      for( my $i = 0; $i < @{$out_idx}; $i++ )
      {
         push @change_protocol, "MERGED_MANY_OUT_ONE_REV: i=$i: $ins_idx->[0] $stoim_val_in $reac_name_in $out_idx->[$i] $stoim_val_out->[$i] $reac_name_out->[$i] => $new_idx->[$i] $reac_name_new->[$i] meta_name=$meta_name meta_idx=$m_idx\n";
      }
   }

   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      print "i=$i: merged: "; print_reaction($network, $new_idx->[$i]);
   }
}
################################################################################


################################################################################
################################################################################
sub merge_oneIn_manyOuts_irrev
{
   my $m_idx = shift;
   my $ins_idx = shift;
   my $out_idx = shift;
   my $network = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};

   my $whoami = _whoami();

   print "reac1: "; print_reaction($network, $ins_idx->[0]);
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      print "i=$i: reac2: "; print_reaction($network, $out_idx->[$i]);
   }

   my $reac_name_in  = $rcs->[$ins_idx->[0]];
   my $reac_name_out = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      push @$reac_name_out, $rcs->[$out_idx->[$i]];
   }

   my $stoim_val_in  = $stm->[$m_idx][$ins_idx->[0]];
   my $stoim_val_out = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      push @$stoim_val_out, $stm->[$m_idx][$out_idx->[$i]];
   }
   die "ERROR: $whoami: stoim_val_in is equal to zero\n"  if $stoim_val_in  == 0;
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      die "ERROR: $whoami: i=$i: stoim_val_out is equal to zero\n" if $stoim_val_out->[$i] == 0;
   }

   divide_matrix_row_by( $stm, $ins_idx->[0], abs($stoim_val_in));
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      divide_matrix_row_by( $stm, $out_idx->[$i], abs($stoim_val_out->[$i]));
   }

   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      add_col2_to_col1_of_matrix( $stm, $out_idx->[$i], $ins_idx->[0] );
   }

   # set new reaction name
   my $reac_name_new = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      $reac_name_new->[$i] = $reac_name_in . $reac_separator . $reac_name_out->[$i];
      $rcs->[$out_idx->[$i]] = $reac_name_new->[$i];
   }

   my $reac_name_extract = splice @$rcs, $ins_idx->[0], 1;

   my $reac_extract      = extract_reac($stm, $ins_idx->[0]);

   # remove reversibility entry from array
   splice @$rvs, $ins_idx->[0], 1;

   print_metabolite($network, $m_idx);
   for( my $r =  0; $r < @{$stm->[$m_idx]}; $r++ )
   {
      die "ERROR: $whoami: r=$r element for metabolite $m_idx is not zero. This should never happen!\n" if $stm->[$m_idx][$r] != 0;
   }

   # remove involved metabolite from name array
   my $meta_name = splice @$mts, $m_idx, 1;

   # remove involved metabolite from stoichiometric matrix
   splice @$stm, $m_idx, 1;

   my $new_idx = [];
   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      $new_idx->[$i] = $out_idx->[$i];
      $new_idx->[$i]-- if $ins_idx->[0] < $out_idx->[$i];
   }

   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      push @change_protocol, "MERGED_ONE_IN_MANY_OUT: i=$i: $ins_idx->[0] $stoim_val_in $reac_name_in $out_idx->[$i] $stoim_val_out->[$i] $reac_name_out->[$i] => $new_idx->[$i] $reac_name_new->[$i] meta_name=$meta_name meta_idx=$m_idx\n";
   }

   for( my $i = 0; $i < @{$out_idx}; $i++ )
   {
      print "i=$i: merged: "; print_reaction($network, $new_idx->[$i]);
   }
}
################################################################################


################################################################################
################################################################################
sub merge_two_irrev_reacs
{
   my $m_idx = shift;
   my $ins_idx = shift;
   my $out_idx = shift;
   my $network = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};
   my $fac1 =  1.0;
   my $fac2 = -1.0;

   my $whoami = _whoami();

   print "reac1: "; print_reaction($network, $ins_idx->[0]);
   print "reac2: "; print_reaction($network, $out_idx->[0]);

   my $reac_name_in  = $rcs->[$ins_idx->[0]];
   my $reac_name_out = $rcs->[$out_idx->[0]];

   my $stoim_val_in  = $stm->[$m_idx][$ins_idx->[0]];
   my $stoim_val_out = $stm->[$m_idx][$out_idx->[0]];
   die "ERROR: $whoami: stoim_val_in is equal to zero\n"  if $stoim_val_in  == 0;
   die "ERROR: $whoami: stoim_val_out is equal to zero\n" if $stoim_val_out == 0;

   divide_matrix_row_by( $stm, $ins_idx->[0], abs($stoim_val_in));
   divide_matrix_row_by( $stm, $out_idx->[0], abs($stoim_val_out));

   add_col2_to_col1_of_matrix( $stm, $ins_idx->[0], $out_idx->[0] );

   # set new reaction name
   my $reac_name_new = $reac_name_in . $reac_separator . $reac_name_out;
   $rcs->[$ins_idx->[0]] = $reac_name_new;

   my $reac_name_extract = splice @$rcs, $out_idx->[0], 1;

   my $reac_extract      = extract_reac($stm, $out_idx->[0]);

   # remove reversibility entry from array
   splice @$rvs, $out_idx->[0], 1;

   # remove involved metabolite from stoichiometric matrix
   print_metabolite($network, $m_idx);

   # remove involved metabolite from name array
   my $meta_name = splice @$mts, $m_idx, 1;

   for( my $r =  0; $r < @{$stm->[$m_idx]}; $r++ )
   {
      die "ERROR: $whoami: r=$r element for metabolite $m_idx is not zero. This should never happen!\n" if $stm->[$m_idx][$r] != 0;
   }
   splice @$stm, $m_idx, 1;

   my $new_idx = $ins_idx->[0];
   $new_idx-- if $out_idx->[0] < $ins_idx->[0];

   push @change_protocol, "MERGED_TWO_IRREV: $ins_idx->[0] $stoim_val_in $reac_name_in $out_idx->[0] $stoim_val_out $reac_name_out => $new_idx $reac_name_new meta_name=$meta_name meta_idx=$m_idx\n";
   update_new_names($reac_name_new, $reac_name_new)   if $do_only_linear_compression;
   update_new_stoic($reac_name_in, $fac1/$stoim_val_in)   if $do_only_linear_compression;
   update_new_stoic($reac_name_out, $fac2/$stoim_val_out) if $do_only_linear_compression;

   print "merged: "; print_reaction($network, $new_idx);
}
################################################################################


################################################################################
################################################################################
sub sub_col2_from_col1_of_matrix
{
   my $stm = shift;
   my $idx1 = shift;
   my $idx2 = shift;

   for( my $i = 0; $i < @$stm; $i++ )
   {
      $stm->[$i][$idx1] -= $stm->[$i][$idx2];
   }

}
################################################################################


################################################################################
################################################################################
sub add_col2_to_col1_of_matrix
{
   my $stm = shift;
   my $idx1 = shift;
   my $idx2 = shift;

   for( my $i = 0; $i < @$stm; $i++ )
   {
      $stm->[$i][$idx1] += $stm->[$i][$idx2];
   }

}
################################################################################


################################################################################
################################################################################
sub divide_matrix_row_by
{
   my $stm = shift;
   my $idx = shift;
   my $div = shift;

   for( my $i = 0; $i < @$stm; $i++ )
   {
      $stm->[$i][$idx] /= $div;
   }

}
################################################################################


################################################################################
################################################################################
sub multiply_matrix_row_by
{
   my $stm   = shift;
   my $idx   = shift;
   my $multi = shift;

   for( my $i = 0; $i < @$stm; $i++ )
   {
      $stm->[$i][$idx] *= $multi;
   }

}
################################################################################


################################################################################
################################################################################
sub extract_reac
{
   my $stm = shift;
   my $idx = shift;
   my $reac;

   for( my $i = 0; $i < @$stm; $i++ )
   {
      my $arr_ref;
      push @$arr_ref, splice @{$stm->[$i]}, $idx, 1;
      push @$reac, $arr_ref;
   }

   return $reac;
}
################################################################################


################################################################################
################################################################################
sub extract_reac_name
{
   my $rcs = shift;
   my $idx = shift;

   my $name = splice @$rcs, $idx, 1;

   return $name;
}
################################################################################


################################################################################
################################################################################
sub read_stoichiomat
{
   my $file = shift;
   my $sm;
   my $num_line     = 0;
   my $last_num_col = -1;

   my $whoami = _whoami();

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";

   while( <$fh> )
   {
      my @st_facs = split;
      my $num_cols = @st_facs;

      push @$sm, \@st_facs;

      $num_line++;
      if( $num_line > 1 && $num_cols != $last_num_col )
      {
         die  "ERROR: $whoami: number of columns ($num_cols) in line $num_line in file '$file' is not equal to number of columns ($last_num_col) in previous line\n";
      }
      
      $last_num_col = @st_facs;
   }

   close $fh;

   return $sm;
}
################################################################################


################################################################################
################################################################################
sub read_reversibility
{
   my $file = shift;
   my $rvs;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
   $_ = <$fh>;
   @$rvs = split;
   close $fh;

   return $rvs;
}
################################################################################


################################################################################
################################################################################
sub read_metabolites
{
   my $file = shift;
   my $mbs;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
   $_ = <$fh>;
   s/"//g;
   # s/_//g;
   s/#//g;
   @$mbs = split;
   close $fh;

   return $mbs;
}
################################################################################


################################################################################
################################################################################
sub read_no_one2many_comp
{
   my $file = shift;
   my @rcs;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
   $_ = <$fh>;
   s/"//g;
   s/>//g;
   # s/_//g;
   s/#//g;
   @rcs = split;
   close $fh;

   return @rcs;
}
################################################################################


################################################################################
################################################################################
sub read_reactions
{
   my $file = shift;
   my $rcs;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
   $_ = <$fh>;
   s/"//g;
   s/>//g;
   # s/_//g;
   s/#//g;
   @$rcs = split;
   close $fh;

   ############################################################################
   # find a suitable record separator for merged reaction names
   ############################################################################
   $reac_separator = '';
   for( my $s = 0; $s < @potential_reac_separartors; $s++ )
   {
      my $sep = $potential_reac_separartors[$s];
      print "INFO: testing reaction separator '$sep'\n";
      my $found_sep = 0;
      for( my $r = 0; $r < @$rcs; $r++ )
      {
         for( my $i = 0; $i < length $rcs->[$r]; $i++ )
         {
            if( $sep eq substr $rcs->[$r], $i, 1 )
            {
               print "INFO: found potential separator '$sep' in reaction '$rcs->[$r]' at position $i -> this separator is not suitable\n";
               $found_sep = 1;
               last ;
            }
         }
         last if $found_sep;
      }
      if( $found_sep == 0 )
      {
         print "INFO: we found a suitable record separator: '$sep'\n";
         $reac_separator = $sep;
         last;
      }
   }

   if( $reac_separator eq '' )
   {
      warn "ERROR: list of potential reactor separators (@potential_reac_separartors) did not contain a suitable candidate\n";
      die  "       update list and restart program\n";
   }
   ############################################################################

   push @change_protocol, "REACTION_SEPARATOR: $reac_separator\n";

   return $rcs;
}
################################################################################


################################################################################
################################################################################
sub update_new_names
{
   my $all_chg  = shift;
   my $new_name = shift;

   my @changed_reacs = split $reac_separator, $all_chg;

   # print "XXXXXXXXXXXXXXXXXXXXXXX\n";
   # print "@change_protocol";
   # print "XXXXXXXXXXXXXXXXXXXXXXX\n";

   foreach my $changed_reac (@changed_reacs)
   {
      die "ERROR: update_new_names(): original reaction '$changed_reac' not found in hash\n" unless exists $reac_new_names->{$changed_reac};

      push @{$reac_new_names->{$changed_reac}}, $new_name;
   }
}
################################################################################


################################################################################
################################################################################
sub update_new_stoic
{
   my $all_chg  = shift;
   my $factor = shift;

   my @changed_reacs = split $reac_separator, $all_chg;

   foreach my $changed_reac (@changed_reacs)
   {
      die "ERROR: update_new_stoic(): original reaction '$changed_reac' not found in hash\n" unless exists $reac_new_names->{$changed_reac};

      $reac_new_stoic->{$changed_reac} *= $factor;
   }
}
################################################################################


################################################################################
# read in program options
################################################################################
sub read_arguments
{
   getopts('s:m:r:v:p:o:n:i:hlk');

   if( $opt_h )
   {
      usage();
   }

   if( $opt_l )
   {
      $do_only_linear_compression = 1;
      print "INFO: we only do linear compression and do not compress many-in-one-out or one-in-many-out metabolites\n";
   }

   if( $opt_k )
   {
      $skip_kernel_based_compression = 1;
      print "INFO: we skip all compressed methods based on the kernel/nullspace of the stoichiometric matrix\n";
   }

   if( $opt_p )
   {
      $pfile = $opt_p;
   }
   else
   {
      usage('ERROR: name of output file change protocol is written to',-1);
   }

   if( $opt_n )
   {
      $nfile = $opt_n;
   }
   else
   {
      warn "INFO: name of file containing reactions which may only be 'linearly' compressed was not provided\n";
   }

   if( $opt_i )
   {
      $skipfile = $opt_i;
   }
   else
   {
      warn "INFO: name of file containing reactions which may not be compressed at all was not provided\n";
   }

   if( $opt_s )
   {
      $sfile = $opt_s;
   }
   else
   {
      usage('ERROR: name of input file containing stoichiometric matrix not provided ',-1);
   }

   if( $opt_m )
   {
      $mfile = $opt_m;
   }
   else
   {
      usage('ERROR: name of input file containing metabolite names not provided ',-1);
   }

   if( $opt_r )
   {
      $rfile = $opt_r;
   }
   else
   {
      usage('ERROR: name of input file containing reaction names not provided ',-1);
   }

   if( $opt_v )
   {
      $rvfile = $opt_v;
   }
   else
   {
      usage('ERROR: name of input file containing reaction reversibility information not provided ',-1);
   }

   if( $opt_o )
   {
      $ofile_postfix = $opt_o;
   }
   else
   {
      usage('ERROR: postfix for output files not provided',-2);
   }
}
################################################################################


################################################################################
################################################################################
sub write_new_reac_names
{
   my ($new_name_hash) = @_;

   my $fn = $rfile . '.new_reac_names';
   my $fh;

   ############################################################################
   # print change protocol to file
   ############################################################################
   print "INFO: going to write new reaction names to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   foreach my $orig_name (sort keys %$new_name_hash)
   {
      print $fh "$orig_name: $new_name_hash->{$orig_name}->[-1]\n";
   }
   close $fh;
   ############################################################################
}
################################################################################


################################################################################
################################################################################
sub write_new_reac_stoic
{
   my ($new_stoic_hash) = @_;

   my $fn = $rfile . '.new_reac_stoich';
   my $fh;

   ############################################################################
   # print change protocol to file
   ############################################################################
   print "INFO: going to write new reaction stoichiometric coefficient to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   foreach my $orig_name (sort keys %$new_stoic_hash)
   {
      print $fh "$orig_name: ", $new_stoic_hash->{$orig_name}->decimal(), "\n";
   }
   close $fh;
   ############################################################################
}
################################################################################


################################################################################
################################################################################
sub write_protocol_to_file
{
   my @chg_proto = @_;

   my $fn = $pfile;
   my $fh;

   ############################################################################
   # print change protocol to file
   ############################################################################
   print "INFO: going to write change protocol to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   print $fh @chg_proto;
   close $fh;
   ############################################################################
}
################################################################################



################################################################################
################################################################################
sub write_compressed_network
{
   my $comp_nw = shift;

   my $fn;
   my $fh;

   ############################################################################
   # print reactions to file
   ############################################################################
   $fn = $rfile . $ofile_postfix;
   print "INFO: going to write reactions to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   for( my $i = 0; $i < @{$comp_nw->{reacs}}; $i++ )
   {
      print $fh '"' . "$comp_nw->{reacs}[$i]" . '"';
      print $fh ' ' unless $i == @{$comp_nw->{reacs}} - 1;
   }
   close $fh;
   ############################################################################

   ############################################################################
   # print metabolites to file
   ############################################################################
   $fn = $mfile . $ofile_postfix;
   print "INFO: going to write metabolites to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   for( my $i = 0; $i < @{$comp_nw->{metas}}; $i++ )
   {
      print $fh '"' . "$comp_nw->{metas}[$i]" . '"';
      print $fh ' ' unless $i == @{$comp_nw->{metas}} - 1;
   }
   close $fh;
   ############################################################################

   ############################################################################
   # print reversiblities to file
   ############################################################################
   $fn = $rvfile . $ofile_postfix;
   print "INFO: going to write reversiblities to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   print $fh "@{$comp_nw->{rever}}";
   close $fh;
   ############################################################################

   ############################################################################
   # print stoichiometric matrix to file
   ############################################################################
   my $num_tests = 0;
   my $abs_min = Math::Fraction->new(10000000);
   $fn = $sfile . $ofile_postfix;
   print "INFO: going to write stoichiometric matrix to file '$fn'\n";
   open $fh, ">$fn" or die "ERROR: couldn't open file '$fn' for writing: $!\n";
   select((select($fh), $| = 1)[0]);
   for( my $i = 0; $i < @{$comp_nw->{stoim}}; $i++ )
   {
      for( my $j = 0; $j < @{$comp_nw->{stoim}->[$i]}; $j++ )
      {
         die "i=$i j=$j: value of stoichiometric matrix is not defined\n" unless defined $comp_nw->{stoim}->[$i][$j];
         $comp_nw->{stoim}->[$i][$j]->modify_reduce();
         print $fh $comp_nw->{stoim}->[$i][$j]->num();
         print $fh ' ' if $j < @{$comp_nw->{stoim}->[$i]} - 1;
      }
      print $fh "\n";
      
   }
   close $fh;
   # print "Minium not-zero element in stoichiometric matrix: $abs_min=" . $abs_min->num() . "\n" if $num_tests > 0;
   ############################################################################
}
################################################################################


################################################################################
################################################################################
sub print_metabolite
{
   my $network  = shift;
   my $meta_idx = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};

   print "$mts->[$meta_idx] (id=$meta_idx):";

   for( my $r = 0; $r < @{$stm->[$meta_idx]}; $r++ )
   {
      print " $stm->[$meta_idx][$r] $rcs->[$r] (id=$r, rev=$rvs->[$r])" if $stm->[$meta_idx][$r] != 0;
   }
   print "\n";
}
################################################################################


################################################################################
################################################################################
sub print_reaction
{
   my $network  = shift;
   my $reac_idx = shift;
   my $stm = $network->{stoim};
   my $rcs = $network->{reacs};
   my $mts = $network->{metas};
   my $rvs = $network->{rever};

   print "$rcs->[$reac_idx] (id=$reac_idx, rev=$rvs->[$reac_idx]):";

   # input
   for( my $m = 0; $m < @$mts; $m++ )
   {
      if( $stm->[$m][$reac_idx] < 0 )
      {
         print ' + ' . abs($stm->[$m][$reac_idx]) . "(id=$m) " . $mts->[$m];
      }
   }

   if( $rvs->[$reac_idx] )
   {
      print " <==>";
   }
   else
   {
      print " ==>";
   }


   # output
   for( my $m = 0; $m < @$mts; $m++ )
   {
      if( $stm->[$m][$reac_idx] > 0 )
      {
         print ' + ' . $stm->[$m][$reac_idx] . "(id=$m) " . $mts->[$m];
      }
   }

   print "\n";
}
################################################################################


################################################################################
################################################################################
sub print_matrix_idx
{
   my $mat = shift;
   for( my $m = 0; $m < @$mat; $m++ )
   {
      for( my $r = 0; $r < @{$mat->[$m]->[1]}; $r++ )
      {
         print "$mat->[$m]->[1]->[$r] ";
      }
      print "\n";
   }
}
################################################################################


################################################################################
################################################################################
sub print_matrix
{
   my $mat = shift;
   for( my $m = 0; $m < @$mat; $m++ )
   {
      for( my $r = 0; $r < @{$mat->[$m]}; $r++ )
      {
         print "$mat->[$m][$r] ";
      }
      print "\n";
   }
}
################################################################################


################################################################################
################################################################################
sub print_matrix_for_octave
{
   my $mat = shift;
   print "[";
   for( my $m = 0; $m < @$mat; $m++ )
   {
      for( my $r = 0; $r < @{$mat->[$m]}; $r++ )
      {
         print "$mat->[$m][$r] ";
      }
      print ";";
   }
   print "]\n";
}
################################################################################


################################################################################
################################################################################
sub print_matrix_for_perl
{
   my $mat = shift;
   print "[";
   for( my $m = 0; $m < @$mat; $m++ )
   {
      print "[";
      for( my $r = 0; $r < @{$mat->[$m]}; $r++ )
      {
         print "$mat->[$m][$r], ";
      }
      print "],\n";
   }
   print "];\n";
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

   print "compress_network.pl -s sfile -m file -r rfile -v rvfile -p chg_proto_file -o postfix_filename [-h -l -k -n noOnToMany_file]\n";
   print "\n";
   print "-s ..... name of file containing stoichiometric matrix (input)\n";
   print "-m ..... name of file containing metabolites (input)\n";
   print "-r ..... name of file containing eactions (input)\n";
   print "-v ..... name of file containing reversibility information (input)\n";
   print "-p ..... name of file change protocol is written to (output)\n";
   print "-o ..... postfix for output files\n";
   print "-l ..... do only linear compression and no many-in-one-out or one-in-many-out\n";
   print "-k ..... skip kernel/nullspace based compression methods\n";
   print "-n ..... name of file containing reactions which may only be 'linearly' compressed (input)\n";
   print "-i ..... name of file containing reactions which may not be compressed at all (input)\n";
   print "-h ..... print this message\n";
   print "\n";
   print "compress_network.pl compresses provided metabolic network\n";

   exit($exit_code);
}
################################################################################
