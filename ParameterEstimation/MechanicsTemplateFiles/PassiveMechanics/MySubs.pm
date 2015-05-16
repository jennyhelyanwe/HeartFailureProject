
package MySubs;
require 5.6.0;


=head1 MyUtils

Plotting functions 

=cut

BEGIN {
}

sub findCalciumGivenT
{
	# This perl script is designed to find the corresponding calcium level given T
	my ($index)=@_;
	
	
	# Define the time for one calcium transient
	my @t=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,53,60,70,75,85,90,93,95,100,105,110,112,125,140,150,165,180,200,220,260,280,300,330,350,370,390,410,430,445,465,480,500);
	
	my $maxIndex=$#t;
	
	print "Size of the array Time is: $maxIndex \n";
	
	# Define the calcium concentration at a given time
	my @Calcium=(0,12.4,24.8,37.2,49.6,62,74.4,86.8,99.2,111.6,124,136.4,148.8,161.2,173.6,186,198.4,210.8,223.2,235.6,248,260.4,272.8,285.2,297.6,310,322.4,334.8,347.2,359.6,372,384.4,396.8,409.2,421.6,434,446.4,458.8,471.2,483.6,496,508.4,520.8,533.2,545.6,558,570.4,582.8,595.2,607.6,620,610,600,580,572,565,558,550,535,515,500,483,467,450,433,417,400,383,367,350,333,317,300,283,267,250,233,217,200,183,167,150,133,117,100,80,65,50,40,35,30,25,20,15,10,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	
	my $ca_t=0; 
	my $current_t=0;
	
	$current_t=@t[$index];
	
	$ca_t=@Calcium[$index];
	
	if ($index>$maxIndex)
	{
		print " Maximum index of the calcium array has been reached \n"
	} else
	{  
		print "Current Time is: $current_t, Current Calium is: $ca_t\n";
	}
	
	return ($ca_t);
}

sub findCalciumGivent3
{
	# This perl script is designed to find the corresponding calcium level given T
	my ($index)=@_;
	
	
	# Define the time values for a single calcium transient
	my @t=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,53,60,70,75,85,90,93,95,100,105,110,112,125,140,150,165,180,200,220,260,280,300,330,350,370,390,410,430,445,465,480,500);

	#my $maxIndex=$#t+1;
	
	#print "Size of the array T is: $maxIndex \n";
	

	my @Calcium=(0,12.4,24.8,37.2,49.6,62,74.4,86.8,99.2,111.6,124,136.4,148.8,161.2,173.6,186,198.4,210.8,223.2,235.6,248,260.4,272.8,285.2,297.6,310,322.4,334.8,347.2,359.6,372,384.4,396.8,409.2,421.6,434,446.4,458.8,471.2,483.6,496,508.4,520.8,533.2,545.6,558,570.4,582.8,595.2,607.6,620,610,600,580,572,565,558,550,535,515,500,483,467,450,433,417,400,383,367,350,333,317,300,283,267,250,233,217,200,183,167,150,133,117,100,80,65,50,40,35,30,25,20,15,10,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	
	my $maxIndex=$#Calcium+1;
	
	print "Size of the array Calcium is: $maxIndex \n";
	

	my $CA_t=0; 
	my $current_t=0;
	
	$current_t=@t[$index-1];
	
	$CA_t=@Calcium[$index-1];
	#$CA_t=$CA_t-100; #Offset calcium by -100
	
	
	if ($index>$maxIndex)
	{
		print " Maximum index of the calcium array has been reached \n"
		
		#$CA_t=0;
		
	} else 
	
	{ 
		print "Current Time is: ${current_t}, Current Calcium is: ${CA_t}\n";	
	}
	
	return ($CA_t);
}


	
sub findVolumeGivenPEspen
{
  # findVolumeGivenP finds (interpolates linearly) volume for a given pressure.
  my ($p0,$V_ED,$EF_at_peak_P,$P_start_ejection) = @_;
  
  my @Peject=(10.0,10.2,10.4,10.6,10.8,11.0,11.2,11.4,11.6,11.8,
           12.0,12.2,12.4,12.6,12.8,13.0,13.2,13.4,13.6,13.8,
           14.0,14.2,14.4,14.6,14.8,15.0);
  my @Veject=(1.0000,0.9800,0.9601,0.9401,0.9202,0.8992,0.8748,
           0.8504,0.8262,0.8070,0.7878,0.7662,0.7441,0.7219,
           0.6844,0.6301,0.6008,0.5548,0.5067,0.4697,0.4299,
           0.3670,0.3058,0.2367,0.1458,0.0000);
  my $maxIndex = $#Peject;
  
  my $V_at_peak_P = $V_ED*(100-$EF_at_peak_P)/100;
  my $index=0;
  my $v=0; my $p1=0; my $p2=0; my $v1=0; my $v2=0;my $a=0;
  
  foreach my $v (@Veject)   #find volumes based on EDV and ejection fraction
  {
    $Veject[$index]=$v*($V_ED-$V_at_peak_P) + $V_at_peak_P;
    $index++;
  }
  if($P_start_ejection != 10)
  {
    $index=0;
    foreach my $p (@Peject)
    {
      $Peject[$index]=($p-10)/5*(15-$P_start_ejection)+$P_start_ejection;
      $index++;
    }
  }
  
  $index=0;
  foreach my $p (@Peject) #find volume at current press by linear interpolation
  {
    if($p>=$p0)
    {
      $p2=$p;
      $v2=$Veject[$index];
      if($index>0)
      {
        $p1=$Peject[$index-1];
        $v1=$Veject[$index-1];
      } else
      {
        $p1=$p;
        $v1=$Veject[$index];
      }
      last;
    }
    $index++;
  }
  if($p2==$p1)
  {
    $a = 0;
  } else
  {
    $a = ($v2-$v1)/($p2-$p1);
  }
  $v = $a*($p0-$p1) + $v1;
 
  if($p0 > $Peject[$maxIndex])
  {
    print "Warning: Subroutine findVolumeGivenP: input pressure exceeds max pressure!\n";
    $v=$Veject[$maxIndex];
  }
  print "p1: $p1, p2: $p2, v1: $v1, v2: $v2, v: $v\n";
  return($v); 
}

sub findVolumeGivenP
{
	# This perl script is designed to find the corresponding LV cavity volume given pressure
	my ($index)=@_;
	
	# Define the pressure extracted from the pressure recordings
	my @p=(11.69,11.75,11.76,11.72,11.73,11.75,11.80,11.82,11.85,11.86,11.85,11.81,11.72,11.61,11.45,11.31,11.14,11.00,10.82,10.57);
	
	my $maxIndex=$#p;
	
	print "Size of the array Pressure is: $maxIndex \n";
	print "Current index is:  $index \n";
	
	
	# Define the LV cavity volume for a given pressure
	my @LV_Volume=(19.22,18.51,17.45,17.20,16.07,14.97,13.27,11.65,11.50,10.79,10.27,9.14,8.02,7.61,7.04,6.41,5.91,5.20,4.80,4.20);
	
	# Correct LV End-Diastolic Volume
	#my @LV_Volume=(21.26,18.51,17.45,17.20,16.07,14.97,13.27,11.65,11.50,10.79,10.27,9.14,8.02,7.61,7.04,6.41,5.91,5.20,4.80,4.20);
	
	if ($index > $maxIndex)
	{
		print " Maximum index of the pressure array has been reached \n";
	} 
	
	my $p_current=0; 
	my $p_next=0;
	my $current_volume=0;
	
	if ($index == 1)
	{
		$p_next=@p[$index-1];
		$current_volume=@LV_Volume[$index-1];
		print "Pressure is $p_next, Cavity Volume for this pressure is $current_volume \n";
		return ($current_volume,$p_next);
	} else
	{
	$p_current=@p[$index-2];
	$p_next=@p[$index-1];
	$current_volume=@LV_Volume[$index-1];
	print "Current Pressure is $p_current, Next Presure is $p_next, Cavity Volume for this pressure is $current_volume \n";
	return ($current_volume,$p_current,$p_next);
	}
	
}

# changeCalcium
# Function creates/edit the file calcium0xx.ipacti and enters the received parameter as the "Enter initial calcium level [Ca]i [0]: xxx"
sub changeCalcium
{
  my ($Cactn) = @_;

  my $outputFile = "calcium0xx.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $Cactn\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.10000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   1\n";
  print OUTFID sprintf(" Enter initial calcium level [Ca]i [0]:%f\n", $Cactn);
  CORE::close(OUTFID);
  
  return(1);
} 

sub changeCalciumBaseOff
{
  my ($Cactn) = @_;

  my $outputFile = "calcium0xx.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $Cactn\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.10000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.5000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   2\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..4\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*0);
  print OUTFID " Enter element #s/name [EXIT]: 5..16\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn);
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  
  return(1);
} 

sub changeCalciumEjection
{
  my ($Cactn) = @_;

  my $outputFile = "calcium0xxEjection.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $Cactn\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.80000D+02\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   1\n";
  print OUTFID sprintf(" Enter initial calcium level [Ca]i [0]:%f\n", $Cactn);
  CORE::close(OUTFID);
  
  return(1);
} 

sub changeCalciumIan
{
  my ($CA) = @_;

  my $outputFile = "calcium0xx.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.11500D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.31000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.62000D+00\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   1\n";
  #print OUTFID " Enter element #s/name [EXIT]: 1,2,3,4\n";
  #print OUTFID " The value is [0]: 0\n";
  #print OUTFID " Enter element #s/name [EXIT]: 5..16\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $CA);
  #print OUTFID " Enter element #s/name [EXIT]: 0\n";  
  print OUTFID sprintf(" Enter initial calcium level [Ca]i [0]:%f\n", $CA);
  CORE::close(OUTFID);
  
  return(1);
} 
sub changeCalciumFixApex
{
  my ($CA) = @_;

  my $outputFile = "calcium0xx.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.10000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..4\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..27\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $CA*0);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 5..16\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..27\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $CA);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeCalciumEspen
{
  my ($Cactn,$ca_index) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "ActivationStep/calcium0xx_".$ca_index.".ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.1000000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  #print OUTFID " Enter element #s/name [EXIT]: 1..4\n";
  #print OUTFID " Gauss Point #s[EXIT]: 1..4,17..20,33..36,49..52\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_1);
  #print OUTFID " Gauss Point #s[EXIT]: 5..8,21..24,37..40,53..56\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_2);
  #print OUTFID " Gauss Point #s[EXIT]: 9..12,25..28,41..44,57..60\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_3);
  #print OUTFID " Gauss Point #s[EXIT]: 13..16,29..32,45..48,61..64\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_4);
  #print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..16\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.0);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeCalciumEspen_TCa
{
  my ($Cactn,$ca_index) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "output_debug/TCa_".$ca_index.".ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.1000000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..16\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeCalciumEspen_ESTCa
{
  my ($Cactn) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "calcium0xx_OptTCa_Initial.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.1000000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..16\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeCalciumEspen_Steps
{
  my ($Cactn,$Tref_tmp,$Beta_tmp) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "calcium0xx_Tref_tmp.ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID sprintf(" Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: %f\n",$Tref_tmp);
  print OUTFID sprintf(" Enter non-dimensional slope parameter (beta) [1.45]: %f\n",$Beta_tmp);
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  #print OUTFID " Enter element #s/name [EXIT]: 1..4\n";
  #print OUTFID " Gauss Point #s[EXIT]: 1..4,17..20,33..36,49..52\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_1);
  #print OUTFID " Gauss Point #s[EXIT]: 5..8,21..24,37..40,53..56\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_2);
  #print OUTFID " Gauss Point #s[EXIT]: 9..12,25..28,41..44,57..60\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_3);
  #print OUTFID " Gauss Point #s[EXIT]: 13..16,29..32,45..48,61..64\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_4);
  #print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..16\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.0);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeTref
{
  my ($Tref,$Cactn,$ca_index) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "calcium0xx_Tref_".$ca_index.".ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: $Tref\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..7\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..4,17..20,33..36,49..52\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_1);
  print OUTFID " Gauss Point #s[EXIT]: 5..8,21..24,37..40,53..56\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_2);
  print OUTFID " Gauss Point #s[EXIT]: 9..12,25..28,41..44,57..60\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_3);
  print OUTFID " Gauss Point #s[EXIT]: 13..16,29..32,45..48,61..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_4);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 8..49\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.0);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeBeta
{
  my ($Beta,$Cactn,$ca_index) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "calcium0xx_Beta_".$ca_index.".ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.10000000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: $Beta\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..7\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..4,17..20,33..36,49..52\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_1);
  print OUTFID " Gauss Point #s[EXIT]: 5..8,21..24,37..40,53..56\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_2);
  print OUTFID " Gauss Point #s[EXIT]: 9..12,25..28,41..44,57..60\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_3);
  print OUTFID " Gauss Point #s[EXIT]: 13..16,29..32,45..48,61..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_4);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 8..49\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.0);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 

sub changeCalciumEspen_Ejection
{
  my ($Cactn,$ca_index) = @_;
  my $Xi2_1=0.69432E-01;
  my $Xi2_2=0.33001E+00;
  my $Xi2_3=0.66999E+00;
  my $Xi2_4=0.93057E+00;
  my $outputFile = "calcium0xx_".$ca_index.".ipacti";
  CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
  print OUTFID " CMISS Version 1.21 ipacti File Version 2\n";
  print OUTFID " Heading: ipacti file. Calcium level: $CA\n";
  print OUTFID "\n";
  print OUTFID " Specify type of contraction model [1]:\n";
  print OUTFID "   (1) SS tension-length-Ca relation\n";
  print OUTFID "  *(2) Hill\n";
  print OUTFID "   (3) Fading memory\n";
  print OUTFID "   (4) Outer hair cell stiffness\n";
  print OUTFID "    1\n";
  print OUTFID " Enter max isometric tension at ext.ratio=1 (Tref) [100kPa]: 0.1000000D+03\n";
  print OUTFID " Enter non-dimensional slope parameter (beta) [1.45]: 0.14500D+01\n";
  print OUTFID " Enter c50 for [Ca]i saturation curve (0<c<1) [0.5]: 0.50000D+00\n";
  print OUTFID " Enter Hill coeff. for [Ca]i saturation curve (h) [3.0]: 0.30000D+01\n";
  print OUTFID " Enter max [Ca]i (Ca_max) [1]: 0.10000D+01\n";
  print OUTFID " Specify whether the initial calcium level [Ca]i is [1]:\n";
  print OUTFID "  (1) Constant spatially\n";
  print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  print OUTFID "  (3) Defined by Gauss points\n";
  print OUTFID "   3\n";
  print OUTFID " Enter element #s/name [EXIT]: 1..7\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..4,17..20,33..36,49..52\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_1);
  print OUTFID " Gauss Point #s[EXIT]: 5..8,21..24,37..40,53..56\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_2);
  print OUTFID " Gauss Point #s[EXIT]: 9..12,25..28,41..44,57..60\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_3);
  print OUTFID " Gauss Point #s[EXIT]: 13..16,29..32,45..48,61..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*$Xi2_4);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  #print OUTFID " Enter element #s/name [EXIT]: 50..59\n";
  #print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.1);
  #print OUTFID " Gauss Point #s[EXIT]: \n";
  print OUTFID " Enter element #s/name [EXIT]: 8..49\n";
  print OUTFID " Gauss Point #s[EXIT]: 1..64\n";
  print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.0);
  print OUTFID " Gauss Point #s[EXIT]: 0\n";
  #print OUTFID " Enter element #s/name [EXIT]: 13..16\n";
  #print OUTFID " Gauss Point #s[EXIT]: 1..4,17..20,33..36,49..52\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*1.0);
  #print OUTFID " Gauss Point #s[EXIT]: 5..8,21..24,37..40,53..56\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*0.9);
  #print OUTFID " Gauss Point #s[EXIT]: 9..12,25..28,41..44,57..60\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*0.8);
  #print OUTFID " Gauss Point #s[EXIT]: 13..16,29..32,45..48,61..64\n";
  #print OUTFID sprintf(" The value is [0]: %f\n", $Cactn*0.7);
  #print OUTFID " Gauss Point #s[EXIT]: 0\n";
  print OUTFID " Enter element #s/name [EXIT]: 0\n";
  CORE::close(OUTFID);
  print $CA, "\n";
  return(1);
} 


sub writeCalciumIVC
{
	my ($index,$ca, $p,$volume) = @_;
	my $outputFile = "RecordingOutputIVC.txt";
  	CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  	print OUTFID sprintf(" Current output index is %d, Current Calcium is %f Current Pressure is %f Current Cavity Volume is %f\n",$index,$ca,$p,$volume);
  	CORE::close(OUTFID);
	return(1);	
}

sub writeCalciumEjection
{
	my ($ca, $p,$index,$v,$ejection_fraction) = @_;
	my $outputFile = "RecordingOutputEjection.txt";
  	CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  	print OUTFID sprintf(" Current output index is %d, Current Calcium is %f Current Pressure is %f Current Volume is %f Current Ejection Fraction is %f\n",$index,$ca,$p,$v,$ejection_fraction);
  	CORE::close(OUTFID);
	return(1);	
}

sub writePVdiastole
{
	my ($p,$v) = @_;
	my $outputFile = "RecordingOutputDiastole.txt";
  	CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  	print OUTFID sprintf("%f  %f\n",$p,$v);
  	CORE::close(OUTFID);
	return(1);	
}

sub writeCalciumSystole
{
	my ($TCa, $p,$v) = @_;
	my $outputFile = "RecordingOutputSystole.txt";
  	CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  	print OUTFID sprintf("%f  %f  %f\n",$TCa,$p,$v);
  	CORE::close(OUTFID);
	return(1);	
}


sub writeconvergedvolume
{
   my ($LV_Press,$LV_V,$frame_number,$LV_V_of_P) = @_;
   my $outputFile = "Ejection_Pressure_Volume.txt";
  CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  print "\n Current Frame is ${frame_number}, Pressure is ${LV_Press}, Converged Cavity Volume is ${LV_V} Desired Cavity Volume is ${LV_V_of_P}\n";
  print OUTFID sprintf(" Current Frame is %d, Pressure is %f, Converged Cavity Volume is %f \n",$frame_number,$LV_Press,$LV_V);
  CORE::close(OUTFID);
  return(1);
}


sub writeCalciumLevel
{
  my ($ca, $frame) = @_;

  my $outputFile = "Ca_Vs_frame.txt";
  CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  print OUTFID sprintf(" Current Frame is %d, Converged Calcium is %f\n",$frame,$ca);
  CORE::close(OUTFID);
  return(1);
}

sub writeSolving1
{
	# Write the calcium changes for frame 1
  my ($ca_old,$ca,$LV_V,$frame_number) = @_;

  my $outputFile = "SolvingProgress.txt";
  CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  print OUTFID sprintf(" Current Frame is %d\n",$frame_number);
  print OUTFID sprintf(" Current Cavity Volume is %f, Increasing calcium from %f to %f to reduce volume \n",$LV_V,$ca_old,$ca);
  CORE::close(OUTFID);
  return(1);
}

sub writeSolving2
{
  my ($ca_old, $ca,$LV_V) = @_;

  my $outputFile = "SolvingProgress.txt";
  CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  print OUTFID sprintf(" Current Cavity Volume is %f, Increasing calcium from %f to %f to reduce volume \n",$LV_V,$ca_old,$ca);
  CORE::close(OUTFID);
  return(1);
}

sub writeSolving3
{
  my ($ca_old, $ca,$LV_V) = @_;

  my $outputFile = "SolvingProgress.txt";
  CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  print OUTFID sprintf(" Current Cavity Volume is %f, Decreasing calcium from %f to %f to increase volume \n",$LV_V,$ca_old,$ca);
  CORE::close(OUTFID);
  return(1);
}

sub writeSolving4
{
  my ($ca,$LV_V) = @_;

  my $outputFile = "SolvingProgress.txt";
  CORE::open(OUTFID, ">> $outputFile") || die "Error: Couldn't open $outputFile\n";
  print OUTFID sprintf(" Current Cavity Volume is %f, Converged Calcium Concentration was found from %f \n",$LV_V,$ca);
  CORE::close(OUTFID);
  return(1);
}

sub createNewIVCInitial
{
	my $filenameIn="EndOfInf_BeforeIVC_tmp.ipinit";
	my $filenameIn2="output_FM_GAUSS/EndOfInfla_FixAX.ipinit";
	my $filenameOut="EndOfInf_BeforeIVC.ipinit";
	my $line;
	my $flag=0;
	CORE::open (INPUT_FILE,"<$filenameIn")||(die "Could not open file $filenameIn");
	CORE::open (INPUT_FILE2,"<$filenameIn2")||(die "Could not open file $filenameIn2");
	CORE::open (OUTPUT_FILE,">$filenameOut")||(die "Could not open file $filenameOut");
	
	while (defined($line=<INPUT_FILE>))
	{		
		print OUTPUT_FILE $line;
		
	}
	
	while (defined($line=<INPUT_FILE2>))
	{
		if($line =~ m/ Force boundary conditions/) 
		{
			$flag=1;
		}
		if ($flag==1)
		{
			print OUTPUT_FILE $line;
		}
	}
}

sub createNewEjectionInitial
{
	my $filenameIn="EndOfIVC_BeforeEjection_tmp.ipinit";
	my $filenameIn2="output_FM_GAUSS/EndOfIVC.ipinit";
	my $filenameOut="EndOfIVC_BeforeEjection.ipinit";
	my $line;
	my $flag=0;
	CORE::open (INPUT_FILE,"<$filenameIn")||(die "Could not open file $filenameIn");
	CORE::open (INPUT_FILE2,"<$filenameIn2")||(die "Could not open file $filenameIn2");
	CORE::open (OUTPUT_FILE,">$filenameOut")||(die "Could not open file $filenameOut");
	
	while (defined($line=<INPUT_FILE>))
	{		
		print OUTPUT_FILE $line;
		
	}
	
	while (defined($line=<INPUT_FILE2>))
	{
		if($line =~ m/ Force boundary conditions/) 
		{
			$flag=1;
		}
		if ($flag==1)
		{
			print OUTPUT_FILE $line;
		}
	}
}
	
sub createIpmate
{
	my ($C1,$C3,$C4,$C5) = @_;
 	my $Xi2_1=0.69432E-01;
  	my $Xi2_2=0.33001E+00;
  	my $Xi2_3=0.66999E+00;
  	my $Xi2_4=0.93057E+00;
  	my $outputFile = "LV_CubicOpt.ipmate";
  	CORE::open(OUTFID, ">$outputFile") || die "Couldn't open $outputFile";
	print OUTFID " CMISS Version 2.1  ipmate File Version 2\n";
 	print OUTFID " Heading: Elements created in Perl\n";
 	print OUTFID "                                  \n";
 	print OUTFID " Stresses in constitutive law are referred to [2]:\n";
   	print OUTFID "   (1) Reference (theta) coordinates\n";
   	print OUTFID "   (2) Body (fibre/transverse) coordinates\n";
   	print OUTFID "   (3) Body coordinates with active fibre stress\n";
   	print OUTFID "    3\n";
 	print OUTFID " Specify whether the constitutive law is defined by [1]:\n";
   	print OUTFID "   (1) a Green strain energy function (hyperelastic)\n";
   	print OUTFID "   (2) a stress/strain-rate relation\n";
   	print OUTFID "   (3) Gauss point stresses (grid coupling)\n";
    	print OUTFID "    1\n";
	print OUTFID " Specify whether the strain energy W is given as a function of [3]:\n";
   	print OUTFID "   (1) the principal strain invariants\n";
   	print OUTFID "   (2) the principal extension ratios\n";
   	print OUTFID "   (3) the fibre & transverse strains\n";
    	print OUTFID "    3\n";
 	print OUTFID " Specify the form of strain energy function W [3]:\n";
   	print OUTFID "   (1) Polynomial  in fibre and transverse normal strains\n";
   	print OUTFID "   (2) Exponential in fibre and transverse strains\n";
   	print OUTFID "   (3) Pole-zero in fibre and transverse strains\n";
   	print OUTFID "   (4) Fibre in fluid\n";
   	print OUTFID "   (5) Defined in subroutine USER53\n";
   	print OUTFID "    2\n";
 	print OUTFID " Specify the form of exponential funtion[1]:\n";
   	print OUTFID "   (1) W=C1*exp(Q) where Q=2*C2*(Ef+Es+En)+C3*(Ef^2)+C4*(Es^2+En^2+2*Esn^2)+2*C5*(Efs^2+Enf^2)\n";
   	print OUTFID "   (2) Dr J.W. Holmes distributed fibre formulation\n";
   	print OUTFID "   (3) Tong & Fung skin function\n";
   	print OUTFID "   (4) W=C1*exp(Q) where Q=C2*Ef^2+C3*Es^2+C4*En^2+2*C5*Efs*Esf+2*C6*Efn*Enf+2*C7*Esn*Ens\n";
    	print OUTFID "    1\n";
 	print OUTFID " Enter the number of constitutive law parameters [1]: 5\n";
 	print OUTFID "                                                   \n";
 	print OUTFID " Specify whether the exponential  material parameter 1 is [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
   	print OUTFID "   1\n";
 	print OUTFID sprintf(" The value is [0]: %f\n", $C1);
 	print OUTFID "                                \n";
 	print OUTFID " Specify whether the exponential  material parameter 2 is [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
   	print OUTFID "   1\n";
 	print OUTFID " The value is [0.220E+00]: 0\n";
 	print OUTFID "                            \n";
 	print OUTFID " Specify whether the exponential  material parameter 3 is [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
  	print OUTFID "   1\n";
 	print OUTFID sprintf(" The value is [0]: %f\n", $C3);
 	print OUTFID "                                \n";
 	print OUTFID " Specify whether the exponential  material parameter 4 is [1]:\n";
 	print OUTFID "  (1) Constant spatially\n";
 	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
 	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
 	print OUTFID "  (4) Defined by Gauss points\n";
 	print OUTFID "  (5) Defined by Grid points\n";
  	print OUTFID "   1\n";
 	print OUTFID sprintf(" The value is [0]: %f\n", $C4);
 	print OUTFID "                              \n";
 	print OUTFID " Specify whether the exponential  material parameter 5 is [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
   	print OUTFID "   1\n";
 	print OUTFID sprintf(" The value is [0]: %f\n", $C5);
 	print OUTFID "                           \n";
 	print OUTFID " Specify whether stress-free sarcomere extension ratios are [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
    	print OUTFID "   1\n";
	print OUTFID " The value is [0.100E+01]: 0.10000E+01\n";
 	print OUTFID "                                       \n";
 	print OUTFID " Specify whether time delay (in secs) to activation is [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
    	print OUTFID "   1\n";
 	print OUTFID " The value is [0.000E+00]: 0.00000E+00\n";
 	print OUTFID "                                     \n";
 	print OUTFID " Specify whether the density is [1]:\n";
  	print OUTFID "  (1) Constant spatially\n";
  	print OUTFID "  (2) Piecewise constant (defined by elements)\n";
  	print OUTFID "  (3) Piecewise linear   (defined by nodes)\n";
  	print OUTFID "  (4) Defined by Gauss points\n";
  	print OUTFID "  (5) Defined by Grid points\n";
   	print OUTFID "   1\n";
 	print OUTFID " The value is [0.000E+00]: 0.00000E+00\n";
 	CORE::close(OUTFID);
  	return($C1,$C3,$C4,$C5);
 }

sub findCalciumGivenFrame
{
	# This perl script is designed to find the corresponding calcium level given T
	my ($index)=@_;
	
	
	# Define the time for one calcium transient
	my @frame=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24);
	# corresponding frame number
	#          1,2,3,4,5,6,7,8,11,14,17,20,23,26,29
	
	my $maxIndex=$#frame+1;
	
	print "Size of the array Frame is: $maxIndex \n";
	
	# Define the calcium concentration at a given frame
	my @Calcium=(0,0.05,0.1,0.15,0.18,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.39,0.4,0.405,0.41,0.412,0.42,0.434,0.4536,0.48,0.516,0.5616,0.60,0.68,0.70,0.76,0.8072,0.90,1.0724,1.3376);
	#my @Calcium=(0.05,0.1,0.15,0.18,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.434,0.4536,0.516,0.48,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0);
	
	my $ca_frame=0; 
	my $current_frame=0;
	
	$current_frame=@frame[$index-1];
	
	$ca_frame=@Calcium[$index-1];
	
	if ($index>$maxIndex)
	{
		print " Maximum index of the calcium array has been reached \n"
	} else
	{  
		print "Current Frame is: $current_frame, Current Calium is: $ca_frame\n";
	}
	
	return ($ca_frame);
}

sub findPressureGivenFrameIVC
{
	# This perl script is designed to find the corresponding calcium level given T
	my ($index)=@_;
	
	
	# Define the time for one calcium transient
	my @frame=(1,2,3,4,5,6,7);
	# corresponding frame number
	#          1,2,3,4,5,6,7,8,11,14,17,20,23,26,29
	
	my $maxIndex=$#frame+1;
	
	print "Size of the array Frame is: $maxIndex \n";
	
	
	# Define the pressure increment at a given frame
	my @Pressure=(1.18,1.96,3.29,5.05,7.19,9.29,11.07);
	
	my $Pressure_frame=0; 
	my $current_frmae=0;
	
	$current_frame=@frame[$index-1];
	
	$Pressure_frame=@Pressure[$index-1];
	if ($index == 1)
	{
		$p_next=@Pressure[$index-1];
		print "Pressure is $p_next \n";
		return ($p_next);
	} else
	{
	$p_current=@Pressure[$index-2];
	$p_next=@Pressure[$index-1];
	print "Current Pressure is $p_current, Next Presure is $p_next \n";
	return ($p_current,$p_next);
	}
	
	if ($index>$maxIndex)
	{
		print " Maximum index of the calcium array has been reached \n"
	} else
	{  
		print "Current Frame is: $current_frame, Current Pressure is: $Pressure_frame\n";
	}
	
}

sub findPressureGivenFrameEjection
{
	# This perl script is designed to find the corresponding calcium level given T
	my ($index)=@_;
	
	
	# Define the time for one calcium transient
	my @frame=(1,2,3,4,5,6,7,8,9,10,11,12,13,14);
	# corresponding frame number
	#          1,2,3,4,5,6,7,8,11,14,17,20,23,26,29
	
	my $maxIndex=$#frame+1;
	
	print "Size of the array Frame is: $maxIndex \n";
	
	
	# Define the pressure increment at a given frame
	my @Pressure=(12.10,12.21,12.50,13,13.01,13.17,13.30,13.42,13.32,13.30,12.92,12.50,12.0,11.34);
	
	my $Pressure_frame=0; 
	my $current_frmae=0;
	
	$current_frame=@frame[$index-1];
	
	$Pressure_frame=@Pressure[$index-1];
	if ($index == 1)
	{
		$p_next=@Pressure[$index-1];
		print "Pressure is $p_next \n";
		return ($p_next);
	} else
	{
	$p_current=@Pressure[$index-2];
	$p_next=@Pressure[$index-1];
	print "Current Pressure is $p_current, Next Presure is $p_next \n";
	return ($p_current,$p_next);
	}
	
	if ($index>$maxIndex)
	{
		print " Maximum index of the calcium array has been reached \n"
	} else
	{  
		print "Current Frame is: $current_frame, Current Pressure is: $Pressure_frame\n";
	}
	
}



