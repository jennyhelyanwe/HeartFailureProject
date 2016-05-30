BEGIN {i=-2;}

{
if (match($0,"Node number") && $NF== node1 ) {i=0;}
	
	if (i>=0) {
		if (i==42) {$NF=1.0000000;printf(" ");} #dy/ds2 version 1	
		if (i==52) {$NF=-1.5707963;printf(" ");} #d2y/ds1ds2 version 2
		if (i==60) {$NF=-1.0000000;printf(" ");} #dy/ds2 version 3
		if (i==70) {$NF=1.5707963;printf(" ");} #d2y/ds1ds2 version 4
		if (i==80) {$NF=-1.5707963;printf(" ");} #d2z/ds1ds2 version 1
		if (i==88) {$NF=-1.0000000;printf(" ");} #dz/ds2 version 2
		if (i==98) {$NF=1.5707963;printf(" ");} #d2z/ds1ds2 version 3
		if (i==106) {$NF=1.0000000;printf(" ");} #dz/ds2 version 4
		i++;
	}
	 gsub("--"," ");
 	 printf("%s\n",$0);
}

END{}

