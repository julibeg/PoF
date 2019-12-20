#! /usr/bin/perl
################################################################################
################################################################################
# Author:  Christian Jungreuthmayer
# Date:    Tue Apr  1 09:05:38 CEST 2014
# Company: Austrian Centre of Industrial Biotechnology (ACIB)
################################################################################

use strict;
use warnings;

use constant CONSIDER_ZERO => 1e-09;
use Data::Dumper;
use Getopt::Std;
use vars qw($opt_c $opt_w $opt_o $opt_h $opt_s);

my $verbose = 0;
my $string_len = 4;

################################################################################
################################################################################
my ($cfile,$hfile,$ofile_basis);
my ($num_cancer_cells, $num_healthy_cells);
my ($cancer_data,$healthy_data);
my ($hash_rever, $update_reversibility);
################################################################################


################################################################################
# handle command line arguments
################################################################################
read_arguments();
################################################################################

################################################################################
# read cancer cells input data
################################################################################
my $cancer_inp  = read_meta_inp_file($cfile);
# print_meta_data($cancer_inp);
$num_cancer_cells = @$cancer_inp;
################################################################################

################################################################################
################################################################################
for( my $cc = 0; $cc < $num_cancer_cells; $cc++ )
{
   my $struct;

   my $reacs = read_reactions($cancer_inp->[$cc]{rfile});
   $struct->{reacs} = $reacs;
   $struct->{num_reacs} = scalar @$reacs;
   warn "reacs: @$reacs\n" if $verbose;

   my $reac_map;
   %$reac_map = map{ $_ => 1 } @$reacs;
   $struct->{reac_map} = $reac_map;

   my $metas = read_metabolites($cancer_inp->[$cc]{mfile});
   warn "metas: @$metas\n" if $verbose;
   $struct->{metas} = $metas;
   $struct->{num_metas} = scalar @$metas;

   my $rever = read_reversibility($cancer_inp->[$cc]{rvfile});
   warn "rever: @$rever\n" if $verbose;
   $struct->{rever} = $rever;
   $struct->{num_rever} = scalar @$rever;

   if( $struct->{num_rever} != $struct->{num_reacs} )
   {
      die "reading model $cc: number of reversibilities ($struct->{num_rever})  and number of reactions ($struct->{num_reacs}) do not match\n";
   }
   # print "reading model $cc: number of reversibilities ($struct->{num_rever})  and number of reactions ($struct->{num_reacs})\n";

   my ($biomass_reacs_cancer,$bm_cc_idx,$bm_cc_arr) = read_biomass_reacs_cancer_cells($cancer_inp->[$cc]{Rbm}, $reacs, $reac_map, $rever);
   warn "biomass_reac_cancer: @$biomass_reacs_cancer\n" if $verbose;
   warn "bm_cc_idx=@$bm_cc_idx\n" if $verbose;
   warn "bm_cc_arr=@$bm_cc_arr\n" if $verbose;
   $struct->{bm_reac}   = $biomass_reacs_cancer;
   $struct->{bm_idx} = $bm_cc_idx;
   $struct->{bm_arr} = $bm_cc_arr;

   my $stoim = read_stoichiomat($cancer_inp->[$cc]{sfile});

   if( $struct->{num_metas} != scalar(@$stoim) )
   {
      die "reading cancer model $cc: number of metabolites ($struct->{num_metas}) and number of lines in stoichiometric matrix (",scalar(@$stoim),") do not match\n";
   }
   print "reading cancer model $cc: number of metabolites ($struct->{num_metas}) and number of lines in stoichiometric matrix (",scalar(@$stoim),")\n";

   for( my $r = 0; $r < @{$stoim}; $r++ )
   {
      print "INFO: XXX testing line $r of stoichiometric matrix read from file '$cancer_inp->[$cc]{sfile}'\n";
      for( my $c = 0; $c < @{$stoim->[$r]}; $c++ )
      {
         die "ERROR: read_stoichiomat(): undefined value at line $r and element $c\n" unless defined $stoim->[$r][$c];
      }
   }
   print_matrix("stoich:\n", $stoim) if $verbose;
   $struct->{stoim} = $stoim;

   push @$cancer_data, $struct;
}
################################################################################

################################################################################
my ($n_stoich, $n_reacs, $n_metas, $n_rever) = create_combined_system($cancer_data, $healthy_data);
print "combined reacs: @$n_reacs\n" if $verbose;
print "combined metas: @$n_metas\n" if $verbose;
print "combined rever: @$n_rever\n" if $verbose;
print_matrix("combined stoich:\n", $n_stoich) if $verbose;

my $coupled = find_coupled_reactions($cancer_data, $healthy_data);

write_efmtool_files($n_stoich, $n_reacs, $n_metas, $n_rever,$ofile_basis);

my $name;
$name = $ofile_basis . '_dual.cfile';
print_coupled_reacs_to_file($name, $coupled);

my $tmp1 = transpose($n_stoich);
print_matrix("transposed:\n", $tmp1) if $verbose;

my $tmp2 = append_I_right($tmp1);
print_matrix("I appended:\n", $tmp2) if $verbose;

my $tmp3 = append_Iirrev_right($tmp2, $n_rever);
print_matrix("Iirrev appended:\n", $tmp3) if $verbose;

my $tmp4 = append_target_col($tmp3, $cancer_data);
print_matrix("targets appended:\n", $tmp4) if $verbose;

$name = $ofile_basis . '_dual.sfile';
print_matrix_to_file($name,$tmp4);

my $total_metas;
my $total_reacs;
my $total_revers;
my $total_revers_names;
my $colum_values;

for( my $i = 0; $i < @$n_reacs; $i++ )
{
   push @$total_metas, "\"$n_reacs->[$i]\"";
}

for( my $i = 0; $i < @$n_metas; $i++ )
{
   push @$total_reacs, "\"$n_metas->[$i]\"";
   push @$total_revers, 1;
   push @$total_revers_names, "\"$n_metas->[$i]\"";
}
push @$colum_values, scalar(@$n_metas);

for( my $i = 0; $i < @$n_reacs; $i++ )
{
   push @$total_reacs, "\"v_$n_reacs->[$i]\"";
   push @$total_revers, 1;
   push @$total_revers_names, "\"v_$n_reacs->[$i]\"";
}
push @$colum_values, scalar(@$n_reacs); # number of reaction

my $num_revers = 0;
for( my $i = 0; $i < @$n_reacs; $i++ )
{
   if( $n_rever->[$i] == 0 )
   {
      push @$total_reacs, "\"z_$n_reacs->[$i]\"";
      push @$total_revers, 0;
      $num_revers++;
   }
}
push @$total_reacs, "\"target\"";
push @$total_revers, 0;
push @$colum_values, $num_revers;  # number of irreversible reactions
push @$colum_values, 1;            # one target reaction

$name = $ofile_basis . '_dual.mfile';
print_array_to_file($name, $total_metas);

$name = $ofile_basis . '_dual.rvfile';
print_array_to_file($name, $total_revers);

$name = $ofile_basis . '_dual.rfile';
print_array_to_file($name, $total_reacs);

$name = $ofile_basis . '_dual.vfile';
print_array_to_file($name, $total_revers_names);

$name = $ofile_basis . '_dual.xfile';
print_array_to_file($name, $colum_values);
################################################################################




################################################################################
# look for coupled reaction
# if the have the same name, we assume it is the same reactions and they are coupled
# this subroutine creates an array that contains
################################################################################
sub find_coupled_reactions
{
   my ($c_data,$h_data) = @_;

   my $seen;

   for( my $cc = 0; $cc < @$c_data; $cc++ )
   {
      for( my $r = 0; $r < @{$c_data->[$cc]{reacs}}; $r++ )
      {
         # my $rev_marker = 'r' . $c_data->[$cc]{rever}[$r];
         $hash_rever->{$c_data->[$cc]{reacs}[$r]} = $c_data->[$cc]{rever}[$r];
         push @{$seen->{$c_data->[$cc]{reacs}[$r]}}, $c_data->[$cc]{reacs}[$r];
      }
   }

   return $seen;
}
################################################################################


################################################################################
################################################################################
sub write_efmtool_files
{
   my $stoich = shift;
   my $reacs  = shift;
   my $metas  = shift;
   my $rever  = shift;
   my $fbasis = shift;
   my $fh;


   if( $update_reversibility )
   {
      print "INFO: reversibility IS set to 1 if any of the coupled reaction is reversible\n";
      my $rvfile = $fbasis . "_efmtool.rvfile";
      print "INFO: going to write file '$rvfile'\n";
      open $fh, ">$rvfile" or die "ERROR: couldn't open file '$rvfile' for writing: $!\n";
      my $num_elems = @$rever;
      my $num_switch_reversibility = 0;
      for( my $i = 0; $i < $num_elems; $i++ )
      {
         my $reversibility = $hash_rever->{$reacs->[$i]};
         (my $re = $reacs->[$i]) =~ s/_[ch]\d\d\d\d//;
         if( exists $coupled->{$re} && @{$coupled->{$re}} > 1 )
         {
            for( my $j = 0; $j < @{$coupled->{$re}}; $j++ )
            {
               my $reac_name = $coupled->{$re}[$j];
               next if $reacs->[$i] eq $reac_name;

               if( $hash_rever->{$reacs->[$i]} == 0 && $hash_rever->{$reac_name} == 1 )
               {
                  $reversibility = 1;
                  $num_switch_reversibility++;
                  last;
               }
            }
         }
         print $fh $reversibility;
         print $fh " " if $i != $num_elems - 1;
      }
      print "num_switch_reversibility=$num_switch_reversibility\n";
      close $fh;
   }
   else
   {
      print "INFO: reversibility is NOT set to 1 if any of the coupled reaction is reversible\n";
      my $rvfile = $fbasis . "_efmtool.rvfile";
      print "INFO: going to write file '$rvfile'\n";
      open $fh, ">$rvfile" or die "ERROR: couldn't open file '$rvfile' for writing: $!\n";
      print $fh "@$rever";
      close $fh;
   }


   my $rfile = $fbasis . "_efmtool.rfile";
   print "INFO: going to write file '$rfile'\n";
   open $fh, ">$rfile" or die "ERROR: couldn't open file '$rfile' for writing: $!\n";
   print $fh join(" ",map{'"' . $_ . '"'} @$reacs);
   close $fh;

   my $mfile = $fbasis . "_efmtool.mfile";
   print "INFO: going to write file '$mfile'\n";
   open $fh, ">$mfile" or die "ERROR: couldn't open file '$mfile' for writing: $!\n";
   # print $fh "@$metas";
   print $fh join(" ",map{'"' . $_ . '"'} @$metas);
   close $fh;

   my $sfile = $fbasis . "_efmtool.sfile";
   print "INFO: going to write file '$sfile'\n";
   open $fh, ">$sfile" or die "ERROR: couldn't open file '$sfile' for writing: $!\n";
   for( my $r = 0; $r < @$n_stoich; $r++ )
   {
      for( my $c = 0; $c < @{$n_stoich->[$r]}; $c++ )
      {
         if( !defined($n_stoich->[$r][$c]) )
         {
            die "undefined value n_stoich for c=$c and r=$r; (number of columns=",scalar(@{$n_stoich->[$r]}),")\n";
         }
         print $fh " " if $c != 0;
         print $fh $n_stoich->[$r][$c];
      }
      print $fh "\n";
   }
   close $fh;
}
################################################################################


################################################################################
################################################################################
sub create_combined_system
{
   my ($c_data) = @_;

   my ($stoich, $reacs, $metas, $rever);

   my $row = 0;
   my $col = 0;

   #############################################################################
   # cancer cells
   #############################################################################
   for( my $cc = 0; $cc < @$c_data; $cc++ )
   {
     for( my $r = 0; $r < @{$c_data->[$cc]{reacs}}; $r++ )
     {
        push @$reacs, $c_data->[$cc]{reacs}[$r];
     }
     for( my $m = 0; $m < @{$c_data->[$cc]{metas}}; $m++ )
     {
        push @$metas, $c_data->[$cc]{metas}[$m];
     }
     for( my $v = 0; $v < @{$c_data->[$cc]{rever}}; $v++ )
     {
        push @$rever, $c_data->[$cc]{rever}[$v]
     }

     my $num_new_rows = @{$c_data->[$cc]{metas}};
     my $num_new_cols = @{$c_data->[$cc]{reacs}};
     print "num_new_rows=$num_new_rows\n"; # if $verbose;
     print "num_new_cols=$num_new_cols\n"; # if $verbose;

     for( my $r = $row; $r < $row + $num_new_rows; $r++ )
     {
        for( my $c = 0; $c < $col; $c++ )
        {
           # print "r=$r c=$c 0.0\n";
           $stoich->[$r][$c] = 0.0;
        }
     }

     for( my $r = $row, my $rid = 0; $r < ($row + $num_new_rows); $r++, $rid++ )
     {
        # print "r=$r row=$row col=$col num_new_rows=$num_new_rows num_new_cols=$num_new_cols\n";
        for( my $c = $col, my $cid = 0; $c < ($col + $num_new_cols); $c++, $cid++ )
        {
           die "AAA: undefined value cc=$cc rid=$rid cid=$cid r=$r c=$c row=$row col=$col num_new_rows=$num_new_rows num_new_cols=$num_new_cols\n" unless defined $c_data->[$cc]{stoim}[$rid][$cid];
           $stoich->[$r][$c] = $c_data->[$cc]{stoim}[$rid][$cid];
        }
     }

     for( my $r = 0; $r < $row; $r++ )
     {
        for( my $c = $col; $c < $col + $num_new_cols; $c++ )
        {
           # print "r=$r c=$c 0.0\n" if $verbose;
           $stoich->[$r][$c] = 0.0;
        }
     }

     $col += $num_new_cols;
     $row += $num_new_rows;
     print "row: $row\n"; # if $verbose;
     print "col: $col\n"; # if $verbose;
   }
   #############################################################################

   return $stoich, $reacs, $metas, $rever;
}
################################################################################


################################################################################
################################################################################
sub print_meta_data
{
   my $data = shift;

   for( my $i = 0; $i < @$data; $i++ )
   {
      print "i: $i\n";
      foreach my $key (keys %{$data->[$i]})
      {
         print "   $key -> $data->[$i]{$key}\n";
      }
   }
}
################################################################################


################################################################################
################################################################################
sub read_meta_inp_file
{
   my $filename = shift;
   my $meta_data;

   open my $fh, $filename or die "ERROR: couldn't openfile '$filename' for reading: $!\n";

   while(<$fh>)
   {
      my $struct;
      # print "Read line form '$filename':", $_, "\n";
      s/^\s*//;
      s/\s*$//;

      next if /^$/;
      next if /^#/;
      # print "Read line form '$filename':", $_, "\n";

      my @fields = split ',';
      die "ERROR: invalid number of elements in line '$_'\n" unless @fields == 5;
      # print "read fields: @fields\n";

      $struct->{rfile}  = compose_filename($filename, $fields[0]);
      $struct->{mfile}  = compose_filename($filename, $fields[1]);
      $struct->{sfile}  = compose_filename($filename, $fields[2]);
      $struct->{rvfile} = compose_filename($filename, $fields[3]);
      $struct->{Rbm}    = compose_filename($filename, $fields[4]);

      push @$meta_data, $struct;
   }

   close $fh;

   # print Data::Dumper->Dump($meta_data);
   # print_meta_data($meta_data);

   return $meta_data;
}
################################################################################


################################################################################
################################################################################
sub compose_filename
{
   my $meta_filename = shift;
   my $field_filename = shift;

   my $new_filename;

   $field_filename =~ s/^\s*//;
   $field_filename =~ s/\s*$//;

   if( $field_filename =~ /^\// )
   {
      $new_filename = $field_filename;
   }
   else
   {
      if( $meta_filename =~ /(.*\/)/ )
      {
         $new_filename = $1 . $field_filename;
      }
      else
      {
         $new_filename = $field_filename;
      }
   }

   # print "DEBUG: compose_filename(): meta_filename=$meta_filename, field_filename=$field_filename: new_filename=$new_filename\n";
   return $new_filename;
}
################################################################################


################################################################################
################################################################################
sub append_target_col
{
   my $in  = shift;
   my $c_data = shift;
   my $out;

   # copy matrix to new variable
   for( my $m = 0; $m < @$in; $m++ )
   {
      for( my $r = 0; $r < @{$in->[$m]}; $r++ )
      {
         $out->[$m][$r] = $in->[$m][$r];
      }
   }

   # add a column at last position of matrix
   my $c = @{$in->[0]};

   for( my $r = 0; $r < @$in; $r++ )
   {
      $out->[$r][$c] =  0;
   }

   # set cells at target rows to -1
   my $row = 0;
   for( my $cc = 0; $cc < @$c_data; $cc++ )
   {
      for( my $r = 0; $r < @{$c_data->[$cc]{bm_arr}}; $r++ )
      {
         if( $c_data->[$cc]{bm_arr}[$r] )
         {
            $out->[$row][$c] = -1;
         }
         $row++;
      }
   }

   return $out;
}
################################################################################


################################################################################
################################################################################
sub append_I_right
{
   my $in  = shift;
   my $out = [];

   # copy upper part
   my $m = 0;
   my $r = 0;
   for( ; $m < @$in; $m++ )
   {
      # print "m=$m\n";
      for( $r = 0; $r < @{$in->[$m]}; $r++ )
      {
         # print "   r=$r\n";
         $out->[$m][$r] = $in->[$m][$r];
      }
   }
   my $last = $r;

   for( $m = 0; $m < @$in; $m++ )
   {
      for( my $r = 0; $r < @$in; $r++ )
      {
         $out->[$m][$r+$last] = 0.0;
         $out->[$m][$r+$last] = 1.0 if $m == $r;
      }
   }

   return $out;
}
################################################################################


################################################################################
################################################################################
sub append_Iirrev_right
{
   my $in     = shift;
   my $revers = shift;
   my $out    = [];

   # copy upper part
   my $m = 0;
   my $r = 0;
   for( ; $m < @$in; $m++ )
   {
      # print "m=$m\n";
      for( $r = 0; $r < @{$in->[$m]}; $r++ )
      {
         # print "   r=$r\n";
         $out->[$m][$r] = $in->[$m][$r];
      }
   }
   my $last = $r;

   my $r_cnt = $last;
   for( my $r = 0; $r < @$in; $r++ )
   {
      if( $revers->[$r] == 0 )
      {
         # we dealing with an irrversible reaction -> add column
         for( $m = 0; $m < @$in; $m++ )
         {
            $out->[$m][$r_cnt] = 0.0;
            $out->[$m][$r_cnt] = -1.0 if $m == $r;
         }
         $r_cnt++;
      }
   }

   return $out;
}
################################################################################


################################################################################
################################################################################
sub setLength
{
   my $inp = shift;
   my $len = shift;

   while( length($inp) < $len )
   {
      $inp = '0' . $inp;
   }

   return $inp;
}
################################################################################


################################################################################
################################################################################
sub transpose
{
   my $mat = shift;
   my $trp = [];;

   # number of rows
   my $m = @$mat;
   # number of colums
   my $n = @{$mat->[0]};

   for( my $r = 0; $r < $m; $r++ )
   {
      for( my $c = 0; $c < $n; $c++ )
      {
         $trp->[$c][$r] = $mat->[$r][$c];
      }
   }

   return $trp;
}
################################################################################


################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
################################################################################
sub read_stoichiomat
{
   my $file = shift;
   my $sm;
   my $line_cnt = 0;
   my $num_elems_per_line;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";

   while( <$fh> )
   {
      s/^\s+//;
      s/\s+$//;
      my @st_facs = split;
      push @$sm, \@st_facs;

      $line_cnt++;
      if( $line_cnt == 1 )
      {
         $num_elems_per_line = scalar @st_facs;
      }
      elsif( $num_elems_per_line != scalar @st_facs )
      {
         die "ERROR: read_stoichiomat(): file '$file' in line $line_cnt: invalid number of elements (",scalar(@st_facs),"), expected $num_elems_per_line elements\n";
      }

      for( my $i = 0; $i < @st_facs; $i++ )
      {
         die "$i element in line $line_cnt is not defined!\n" unless defined $st_facs[$i];
      }
   }

   close $fh;

   print "read $line_cnt lines from file '$file' each containing $num_elems_per_line\n";

   for( my $r = 0; $r < $line_cnt; $r++ )
   {
      print "INFO: testing line $r of stoichiometric matrix read from file '$file'\n";
      for( my $c = 0; $c < @{$sm->[$r]}; $c++ )
      {
         die "ERROR: read_stoichiomat(): undefined value at line $r and element $c\n" unless defined $sm->[$r][$c];
      }
   }

   return $sm;
}
################################################################################


################################################################################
################################################################################
sub print_matrix
{
   my $message = shift;
   my $matrix  = shift;

   warn $message;

   for( my $m = 0; $m < @$matrix; $m++ )
   {
      for( my $r = 0; $r < @{$matrix->[$m]}; $r++ )
      {
         print STDERR " ";
         print STDERR " " if $matrix->[$m][$r] >= 0;
         print STDERR $matrix->[$m][$r];
      }
      warn "\n";
   }
}
################################################################################


################################################################################
################################################################################
sub print_coupled_reacs_to_file
{
   my $filename = shift;
   my $coupled   = shift;

   open my $fh, ">$filename" or die "ERROR: couldn't open file '$filename' for writing: $!\n";

   print "INFO: going to write dual cfile to '$filename'\n";
   foreach my $reacs (sort keys %{$coupled})
   {
      if( @{$coupled->{$reacs}} > 1 )
      {
         for( my $i = 0; $i < @{$coupled->{$reacs}}; $i++ )
         {
            print $fh " " if $i != 0;
            print $fh "\"$coupled->{$reacs}[$i]\"";
         }
         print $fh "\n";
         # print $fh join(':',@{$coupled->{$reacs}}),"\n";
      }
   }

   close $fh;
}
###############################################################################


################################################################################
################################################################################
sub print_matrix_to_file
{
   my $filename = shift;
   my $matrix   = shift;

   open my $fh, ">$filename" or die "ERROR: couldn't open file '$filename' for writing: $!\n";

   print "INFO: going to write dual sfile to '$filename'\n";
   for( my $m = 0; $m < @$matrix; $m++ )
   {
      for( my $r = 0; $r < @{$matrix->[$m]}; $r++ )
      {
         print $fh " ";
         print $fh " " if $matrix->[$m][$r] >= 0;
         print $fh $matrix->[$m][$r];
      }
      print $fh "\n";
   }

   close $fh;
}
###############################################################################


################################################################################
################################################################################
sub print_array_to_file
{
   my $filename = shift;
   my $array    = shift;

   open my $fh, ">$filename" or die "ERROR: couldn't open file '$filename' for writing: $!\n";
   print "INFO: going to write dual mfile to '$filename'\n";

   print $fh "@$array";

   close $fh;
}
###############################################################################


################################################################################
################################################################################
sub print_vector
{
   my $mes = shift;
   my $vec = shift;

   print $mes;

   print join(",",@$vec), "\n"
}
################################################################################


################################################################################
################################################################################
sub read_biomass_reacs_cancer_cells
{
   my $file     = shift;
   my $reacs    = shift;
   my $reac_map = shift;
   my $rever    = shift;
   my $bm_cc;
   my $idx;
   my $arr;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
   $_ = <$fh>;
   s/"//g;
   s/>//g;
   # s/_//g;
   s/#//g;
   @$bm_cc = split;
   close $fh;

   foreach my $reac (@$bm_cc)
   {
      die "ERROR: '$reac' not found in list of reactions @$reacs\n" unless defined $reac_map->{$reac};

      for( my $i = 0; $i < @$reacs; $i++ )
      {
         if( $reac eq $reacs->[$i] )
         {
            push @$idx, $i;
            $arr->[$i] = 1;

            # check if this biomass reaction is (by accident) a reversible reaction
            die "ERROR biomass reaction '$reac' is not allowed to be reversible!\n" if $rever->[$i] != 0;
         }
         else
         {
            $arr->[$i] = 0;
         }
      }
   }

   return($bm_cc, $idx, $arr);
}
################################################################################


################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
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

   return $rcs;
}
################################################################################


################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
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
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
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
# read in program options
################################################################################
sub read_arguments
{
   getopts('c:o:h');

   if( $opt_h )
   {
      usage();
   }

   if( $opt_c )
   {
      $cfile = $opt_c;
   }
   else
   {
      usage('ERROR: name of input file containing file names for cancer cells data not provided ',-1);
   }

   if( $opt_o )
   {
      $ofile_basis = $opt_o;
   }
   else
   {
      usage('ERROR: base name for output files not provided',-3);
   }
}
################################################################################


################################################################################
################################################################################
sub usage
{
   my $message   = shift || '';
   my $exit_code = shift || 0;

   print "$message\n" if $message;

   print "create_ccds_files.pl -c cfile -w wfile -o ofile [-s -h]\n";
   print "\n";
   print "-c ..... name of file containing file names for meta input data files (input)\n";
   print "-o ..... base file name for output files to be written (output)\n";
   print "-h ..... print this message\n";
   print "\n";
   print "\n";
   print "create_ccds_files.pl creates files which are then used by defigueiredo.pl\n";

   exit($exit_code);
}
################################################################################
