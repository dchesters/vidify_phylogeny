
$in0 = $ARGV[0];
$in1 = $ARGV[1];
$in2 = $ARGV[2];
$in3 = $ARGV[3];
$in4 = $ARGV[4];
unless($in3 =~ /./){die "\n\nCOMMAND ERROR\n"};

open(IN0, $in0) || die "\nerror 6\n";

while(my $line = <IN0>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /\#(.+)/)	
		{
		$command_class = $1;$command_class =~ s/\s//g;
		$commands{$command_class} .= "$line\n";
		};

	};

close IN0;

####################################################

@classes = keys %commands;@classes = sort @classes;
print "
command classes read:@classes
";
# BRANCH,1  BRANCH,10  BRANCH,11  BRANCH,12  BRANCH,13  BRANCH,14  BRANCH,15  BRANCH,16  BRANCH,17  BRANCH,18  BRANCH,19  BRANCH,2  BRANCH,20  BRANCH,3  BRANCH,4  BRANCH,5  BRANCH,6  BRANCH,7  BRANCH,8  BRANCH,9  
# FRAME_EDGE,1  FRAME_EDGE,10  FRAME_EDGE,2  FRAME_EDGE,3  FRAME_EDGE,4  FRAME_EDGE,5  FRAME_EDGE,6  FRAME_EDGE,7  FRAME_EDGE,8  FRAME_EDGE,9


open(OUT, ">$in1") || die "";

print OUT "
library(plotrix)
library(scales)

current_x_limit <- $in3
current_age <- $in4
if(current_age <= 1)
	{
	printtext <- \"Current era\"
	}else{
	printtext <- paste(current_age, \"mya\", sep = \" \")
	}

colorfunc = colorRamp(c(
	\"darkblue\" ,\"darkblue\" ,\"darkblue\" ,
	\"red\" ,\"red\" ,\"red\" , 
	\"yellow\",\"yellow\", 
	\"white\",\"white\"))


jpeg( \"$in2\" , width = 4000 , height = 4000, res=200)
par(bg= \"black\")
plot(0, type = \"n\", xlim=c(-0.5 + current_x_limit , 0.5 + current_x_limit), ylim=c(0,1), xlab = \"\",ylab = \"\",xaxt = \"n\", yaxt = \"n\", bty = \"n\")
text(0.3 + current_x_limit , 0.9, labels=printtext, col = \"white\", cex = 4)

";

#

# smaller indexes are larger points, plot these first

# ,"FRAME_EDGE,","FRAME_EDGE,","FRAME_EDGE,",

# plot large diffuse frame edges first
foreach my $class("FRAME_EDGE,1","FRAME_EDGE,2","FRAME_EDGE,3","FRAME_EDGE,4","FRAME_EDGE,5","FRAME_EDGE,6","FRAME_EDGE,7")
	{
	my $commands_for_class = $commands{$class}; print OUT "$commands_for_class";
	};

# plot branches
foreach my $class("BRANCH,1","BRANCH,2","BRANCH,3","BRANCH,4","BRANCH,5","BRANCH,6","BRANCH,7","BRANCH,8","BRANCH,9","BRANCH,10",
		"BRANCH,11","BRANCH,12","BRANCH,13","BRANCH,14","BRANCH,15","BRANCH,16","BRANCH,17","BRANCH,18","BRANCH,19","BRANCH,20")
	{
	my $commands_for_class = $commands{$class}; print OUT "$commands_for_class";
	};

# then remaining small points
foreach my $class("FRAME_EDGE,8","FRAME_EDGE,9","FRAME_EDGE,10")
	{
	my $commands_for_class = $commands{$class}; print OUT "$commands_for_class";
	};

print OUT "
dev.off()
";

close OUT;



