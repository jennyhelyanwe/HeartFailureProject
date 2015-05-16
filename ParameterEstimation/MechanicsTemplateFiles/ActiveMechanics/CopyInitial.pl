	use strict;

	my $filenameIn="LV_CubicPreEpiBase_tmp.ipinit";
	my $filenameIn2="output/LV_Inflation_OptC1.ipinit";
	my $filenameOut="LV_CubicPreEpiBase.ipinit";
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


