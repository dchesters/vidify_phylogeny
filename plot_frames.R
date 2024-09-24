
###################################################################################################################################
#
# 
#	plot_frames.R, for depicting a phylogeny as video frames.
#		  
#    	Copyright (C) 2021-2024 Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
##########################################################################################################################################
#
#
#	Dependencies:
#	Perl script sort_commands_for_plotting_frame.pl
#	R libraries plotrix, scales, tuneR
#
#
#
#
#
#
#	
#	CHANGE LOG
#	 
#	2021-SEP-11 (version 1.01)
#	2021-OCT-03 (version 1.02): 	Checks for Perl script dependency in working directory.
#	2024-JUL-06 (version 1.03): 	Terminus point size scaled by number in current frame rather than decendent count.
#	2024-JUL-14 (version 1.04): 	Works on phylograms (non-ultrametric).
#	 				Couple of ways added to speed up where rendering is very slow on large trees.
#	2024-AUG-13 (version 1.05):	Bugfix 0-1 scaling the table. 
#	2024-SEP-14 (version 1.06):	Option to alter video pan speed.
#	 				Some variables automatically optimized according to number of branches input.
#	 
#	 
#	 
#	 
#
#
#
##########################################################################################################################################


library("plotrix");
library(scales)
library(tuneR)

pitchshift 	<- 1 	# default==0, 1== pitch ascends or decends according to child branch delta Y
			# uses different input table
# track_length 	<- 120	# default is 60 seconds, for terminal-rich trees, make longer



if(file.exists("sort_commands_for_plotting_frame.pl") == TRUE)
	{
	print ("required perl script found in working directory")
	}else{
	print("")
	print("ERROR. required perl script seems NOT in working directory.")
	print("please copy sort_commands_for_plotting_frame.pl to current directory and try again. Quitting.")
	print("")
	quit(save = "no")
	}

# video and audio columns
#                    first 5 columns  VIDEO, 6-9 audio     
# print VIDIFY "$branchnumber\t$sum_branchlength\t$y1_proportion\t$sum_branchlength\t$y2_proportion\t" , 
#			"$sum_branchlength\t$y1_proportion\t$sum_branchlength\t$y1_proportion\t$count_terminalsfromnode\n";


if(pitchshift == 1)
	{
	table2 <- read.table("vidify_table2", header=F, row.names = 1, sep = "\t")
	aud.table1 <- table2[ , c(1,2,3,4,5)]
	vid.table1 <- table2[ , c(1,2,3,4,5)]
	}else{
	table2 <- read.table("vidify_table", header=F, row.names = 1, sep = "\t")
	# each have column 9 which is count of terminals derived from each brnach
	# VIDEO:
	vid.table1 <- table2[ , c(1,2,3,4,9)]
	# AUDIO:
	aud.table1 <- table2[ , c(5,6,7,8,9)]
	}
head(table2)

number_branches <- length(table2[,1]); # about 2x number of terminals.


if(number_branches < 160)
{track_length<-30;total_number_frames <- track_length*25;bin_size <- (1 / total_number_frames);pan_speed <- 0.5 # slow pan for small trees
print("tree < 160 branches")
}else{

if(number_branches < 1200)
{track_length <- 60;total_number_frames <- track_length*25;bin_size <- (1 / total_number_frames);pan_speed <- 0.2
# default is 60 seconds, for terminal-rich trees, make longer
print("tree < 1200 branches")
}else{

if(number_branches < 2400)
{track_length <- 60;total_number_frames <- track_length*25;bin_size <- (1 / total_number_frames);pan_speed <- 0.05
print("tree < 2400 branches")
}else{

if(number_branches < 3200)
{track_length <- 60;total_number_frames <- track_length*25;bin_size <- (1 / total_number_frames);pan_speed <- 0.05
print("tree < 3200 branches")
}else{

if(number_branches >= 3200)
{track_length <- 120;total_number_frames <- track_length*25;bin_size <- (1 / total_number_frames);pan_speed <- 0.05
# consider decrease branch width
print("tree > 3200 branches")
}}}}}

# stop()




colorfunc = colorRamp(c(
	"darkblue" ,"darkblue" ,"darkblue" ,  "red" ,"red" ,"red" ,  "yellow","yellow",  "white","white"
	))



# process both tables

root_age <- max(vid.table1[  , c(1,3)]); # in millions years
print(c("root_age", root_age))

max_number_terminals <- max(vid.table1[  , 5])
major_nodes_A <- quantile(vid.table1[  , 5],  probs = 0.995)
major_nodes_B <- (1:length(vid.table1[  , 5]))[ vid.table1[  , 5] >= major_nodes_A ]
vid.table1[ major_nodes_B , 5] <- major_nodes_A;vid.table1[  , 5] <- vid.table1[  , 5] / major_nodes_A
aud.table1[ major_nodes_B , 5] <- major_nodes_A;aud.table1[  , 5] <- aud.table1[  , 5] / major_nodes_A

#######################################

# either parent to child vertical
# branchnumber sum_branchlength y1_proportion sum_branchlength y2_proportion
# or parent to child horizontal
# branchnumber sum_branchlength y2_proportion assign_new_x    y2_proportion

# COLUMN 	1	2	3	4
#		x	y	x	y

# simplfy matters, scale x 0-1, and y 0-1
min.1.3 <- min(vid.table1[  , c(1,3)]); max.1.3 <- max(vid.table1[  , c(1,3)]);
min.2.4 <- min(vid.table1[  , c(2,4)]); max.2.4 <- max(vid.table1[  , c(2,4)]);

temp1 <- vid.table1[  , 1] +  - min.1.3; temp1b <- temp1 / max.1.3
temp2 <- vid.table1[  , 3] +  - min.1.3; temp2b <- temp2 / max.1.3
temp3 <- vid.table1[  , 2] +  - min.2.4; temp3b  <- temp3 / max.2.4
temp4 <- vid.table1[  , 4] +  - min.2.4; temp4b <- temp4 / max.2.4
vid.table1[  , 1] <- temp1b; vid.table1[  , 3] <- temp2b; vid.table1[  , 2] <- temp3b; vid.table1[  , 4] <- temp4b

dont_plot_short_branches <- 0 # default == 0, if large tree and too slow to render, assign 1
dont_plot_any_branches 	<- 0 	# ==1 if extremly slow to render.
dont_plot_tiny_points	<- 0; 	# another way to speed up a bit, only use for very slow renders
tiny_points_cutoff 	<- 0.1

##########################################################################################
if( dont_plot_short_branches == 1 )
	{
	remove_array <- 0
	for(row in 1:length(vid.table1[ , 1 ] )) 
		{
		x1 <- vid.table1[ row, 1 ];y1 <- vid.table1[ row, 2 ];x2 <- vid.table1[ row , 3 ];y2 <- vid.table1[ row, 4 ];
	#	branch_euclid <- dist(matrix(c(x1,x2,y1,y2), nrow=2)); # print(c("branch_euclid",branch_euclid))
	#	if(branch_euclid < 0.005){remove_array[row] <- 1}else{remove_array[row] <- 0}
		# due to print style, best remove vertical branches
		if(x1==x2 & (y2-y1)<100.0){remove_array[row] <- 1}else{remove_array[row] <- 0}
		}
	keep_indices <- (1:length(remove_array))[remove_array == 0]
	print(c("retained:",length(keep_indices)))
	print(c("total rows:",length(vid.table1[ , 1 ] )))
	vid.table1 <- vid.table1[ keep_indices , ]
#	print (rm_indices);
	}
##########################################################################################



 # Parameters 

# NOTE, 2000 is very slow, test on much less.
# one minute of tree traversal plus 10 seconds out
# 25 fps * 60 seconds is 1500
# total_number_frames <- 1500 # default 1500

# 25 * 5 seconds
finishing_frames <- 125 # default; 125 after tree terminals reached, video should not immediatly cut off, add a few frames, stationary.


postsegment <- 5 # default 5




number_rows <- length(vid.table1[ , 1 ] );print(number_rows)





# for video only

terminus_point_scale_type <- 2

# Make table for audify, number_of_frames * max_notes
max_notes <- 20

#phaseout<-0.95 # 0.99 = cex reduces one quarter per second; 0.95 = reduces more than 50 percent

alphas <- c(0.06, 0.06, 0.06, 0.06, 0.06, 0.07, 0.1, 0.2, 0.4, 0.8)
colorindices <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9)
xoffsets <- c(0.010, 0.005, 0.004, 0.003, 0.0025, 0.002, 0.0010, 0.0005, 0.0002, 0.000)

# old ones before flicker effect:
point_cexes <- c(20, 16, 12, 8, 6, 4, 3, 2, 1, 0.8)
point_cexes <- point_cexes*3
pointcex_x_diminish_factor <- 1.25 # larger value==larger points; default 1.5

branchlwd.outer <- 4;branchlwd.inner <- 0.5;

point_diminish <- 0.5 	# default==1, reduce if points too large.
			# nb this value reduces itself during end phase


##################################################################

print(total_number_frames)
system("rm video_frame_commands*");



 system("rm video_frame.*")
 for (frame_number in 1:(total_number_frames+finishing_frames))

# if plotting last few frames, hash out prev 2 commands, and do this:
# for (frame_number in 2879:(total_number_frames+finishing_frames))
	{
#	frame_number <- 50

	if(frame_number < total_number_frames)
		{
		current_frame_x_limit <- bin_size * frame_number
	#	branchlength <- x2-x1;
		}else{
		############################################################################################
		# end phase, stop incrementing x axis
		current_frame_x_limit <- bin_size * (total_number_frames-1)
	#	branchlength <- (bin_size * frame_number)-x1;

		# glow of points denoting branch edge will fade in end phase		
	#	for( i in 1:10 ){point_cexes[i] <- point_cexes[i]*phaseout}
		point_diminish <- point_diminish * 0.95;
		branchlwd.outer <- branchlwd.outer * 0.92; branchlwd.inner <- branchlwd.inner *0.92
		############################################################################################
		}
	current_filename <- paste ( "video_frame." , frame_number , ".jpeg" , sep = "" )
	system("rm commands_for_current_frame");
	print(c("frame_number" , frame_number, "current_frame_x_limit" , current_frame_x_limit))
	frame_speciations <- 0; length_after_speciation <- 0; frame_speciations_count <- 1;
	case1<-0;case2<-0;case3<-0;case4<-0;case5<-0;
	splits_in_current_bin <-0; splits_in_current_count <- 1; 
	current_frame_points_printed <- 0


	# alternative point scaling is how many terminals at each position along Y
	# vid.table1[ row, 1 ];
	# seems not exactly counting what i want, instead count
	# current_frame_terminal_count <- (1:length(vid.table1[ , 1 ]))[ vid.table1[ , 1 ] <= current_frame_x_limit & vid.table1[ , 3 ] >= current_frame_x_limit ];

	# inefficient, but couldnt figure out other way to calculate point size prior to checking them all
	for(row in 1:number_rows) # for current frame, look through entire tree and see what fits limit.
		{
		current_x1 <- vid.table1[ row, 1 ];	current_x2 <- vid.table1[ row , 3 ];
		if(current_x1 <= current_frame_x_limit & current_x2 <= current_frame_x_limit)
			{
			# start and end within current X
			}else{
			if(current_x1 > current_frame_x_limit & current_x2 > current_frame_x_limit)
				{
				}else{
				if(current_x1 <= current_frame_x_limit & current_x2 >= current_frame_x_limit)
					{
					current_frame_points_printed <-	current_frame_points_printed+1
					}
				}
			}
		}

	 current_frame_terminal_count <- current_frame_points_printed;
	 print(c("current_frame_terminal_count:",current_frame_terminal_count))


	for(row in 1:number_rows) # for current frame, look through entire tree and see what fits limit.
		{
		# current branch X left		    	current branch X right
		current_x1 <- vid.table1[ row, 1 ];	current_x2 <- vid.table1[ row , 3 ];	print_line <- 0; print_point <- 0; point_y <- NA
		

		if(current_x1 > (current_frame_x_limit-bin_size) & current_x2 <= current_frame_x_limit)
			{
			splits_in_current_bin[splits_in_current_count] <- row; splits_in_current_count <- splits_in_current_count+1
			}

		if(current_x1 <= current_frame_x_limit & current_x2 <= current_frame_x_limit)
			{
			# start and end within current X
			print_line <- 1;x1<-vid.table1[ row, 1 ];y1<-vid.table1[ row, 2 ];x2<-vid.table1[ row, 3 ];y2<-vid.table1[ row, 4 ];point_y <- y1
			case1 <- case1 +1
			
			}else{

			if(current_x1 > current_frame_x_limit & current_x2 > current_frame_x_limit)
				{
			#	print("both x's outside (RHS) frame")
				case2 <- case2 + 1
				}else{
			
				if(current_x1 <= current_frame_x_limit & current_x2 >= current_frame_x_limit)
					{
					# current_x1 on the left, before frame limit, 
					# current_x2 on the right, beyond frame limit (so plot only to frame limit)
					print_line <- 1; x1 <- current_x1; y1<-vid.table1[ row, 2 ]; x2 <- current_frame_x_limit; y2 <- vid.table1[ row, 4 ]
					print_point <- 1; point_y <- y2 
					case3 <- case3 + 1
					}else{
				
					if(current_x1 >= current_frame_x_limit & current_x2 <= current_frame_x_limit)
						{
						# this one maybe not encountered, check case4 printout
						print_line <- 1; x1 <- current_x2; x2 <- current_frame_x_limit;y1<- vid.table1[ row, 2 ];y2<- vid.table1[ row, 4 ];
						print_point <- 1; point_y <- y1 
						case4 <- case4 +1
						}else{
						print("HUH");print(c(current_x1, current_x2))
						case5 <- case5 +1
						}
					}


				}
	
			}

		# last few frames should be blank
		if(frame_number >= ((total_number_frames+finishing_frames)-3) ){print_line <- 0; print_point <- 0}
		# branches fade out nice, point persist too long, so:
		if(frame_number >= ((total_number_frames+(finishing_frames/2))) ){print_point <- 0}

	#	if(print_line == 1)
	#		{
	#		branch_euclid <- dist(matrix(c(x1,x2,y1,y2), nrow=2)); # print(c("branch_euclid",branch_euclid))
	#		if(dont_plot_short_branches == 1 & branch_euclid < 0.000001){print_line <- 0}
	#		}

		if(dont_plot_any_branches == 1){print_line <- 0}

		if(dont_plot_tiny_points == 1 & current_frame_x_limit >= 0.50)
			{
			if(current_frame_x_limit-current_x1 >= tiny_points_cutoff){print_point<-0}
			}

		####################################################################################################################
		if(print_line == 1)
			{
			#   BRANCH

			startalpha <- 0.01;startlwd <- 10; startcolor <- 0.1;lwd_shift_factor <- 4 # default=2
			fancy_branch <- 2
			if(fancy_branch == 1) # not used
				{
				for(i in 1:20)
					{
				#	colors2<- rgb(colorfunc(startcolor) , maxColorValue=255)
					command_string <- paste( "colors2 <- rgb(colorfunc(" , startcolor , ") , maxColorValue=255);segments(", 
							 x1 , "," , y1 , "," , x2 , "," , y2 , ",alpha( colors2 , ", startalpha , "), lwd = " ,
							 startlwd , ") # BRANCH," , i , sep = "")
					write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")
					startalpha <- startalpha * 1.2;startlwd <- startlwd - lwd_shift_factor; startcolor <- startcolor + 0.04
					}

				}else{ # default

				#	note very thin branches will appear to flicker in the video
					# 0.5 = red; 0.8=yellow; 0.95=white
					command_string <- paste( "colors2 <- rgb(colorfunc(" , 0.5 , ") , maxColorValue=255);segments(", 
							 x1 , "," , y1 , "," , x2 , "," , y2 , ",alpha( colors2 , ", 0.999 , "), lwd = " ,
							 1.5 , ") # BRANCH," , 1 , sep = "")
				#	write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")
					command_string <- paste( "colors2 <- rgb(colorfunc(" , 0.8 , ") , maxColorValue=255);segments(", 
							 x1 , "," , y1 , "," , x2 , "," , y2 , ",alpha( colors2 , ", 0.999 , "), lwd = " ,
							 branchlwd.outer , ") # BRANCH," , 2 , sep = "")
				#	write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")
					command_string <- paste( "colors2 <- rgb(colorfunc(" , 0.95 , ") , maxColorValue=255);segments(", 
							 x1 , "," , y1 , "," , x2 , "," , y2 , ",alpha( colors2 , ", 0.999 , "), lwd = " ,
							 branchlwd.inner , ") # BRANCH," , 3 , sep = "")
					write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")

				} 


			# this branch is within current frame (total, not just column).
			# whether to play tone depends on if the branch is within column,
			# and the volume of the tone is dependent on the length from start of branch to current x limit.

			# record bifurcations at current x column, their y positions will be turned to tones.
			current_column_LHS <- current_frame_x_limit - bin_size
			# is branch RHS within column:
			if(x2 >= current_column_LHS ) 
				{
				      frame_speciations[frame_speciations_count] <- point_y; 
				length_after_speciation[frame_speciations_count] <- x2-x1
				frame_speciations_count <- frame_speciations_count + 1;
				}

			}
		####################################################################################################################
		if(print_point == 1)
			{
			#    POINT DENOTING EDGE OF CURRENT FRAME

			# diminishing flickering effect for branch tips:
			if(frame_number < total_number_frames)
				{
				branchlength <- x2-x1;
				}else{
				branchlength <- (bin_size * frame_number)-x1;
				}
			# initiate point size for branch terminus from sampling a distribution
			x<-rlnorm(101);x<-x/max(x);x<-sort(x, decreasing=T); flicker_val <- x[round(branchlength*100)+1]
			# also, gets a bit ugly at Right hand side (reduces to colored bars), so diminish cex along frame sequence
			# flicker_val <- flicker_val * (1-(current_frame_x_limit/pointcex_x_diminish_factor)); 
			# column 5 of vid.table1 is number of terminals decended from branch

			if(terminus_point_scale_type == 1) # not in use
				{
				# scale to number of terminals derived from node:
			 	flicker_val <- (flicker_val * vid.table1[ row , 5]) * point_diminish; 
				};

			if(terminus_point_scale_type == 2) # default
				{
				# OR, scale by how many branch terminals present at current frame
				# point_diminish gets lower every frame during end phase
				flicker_val <- (flicker_val / log(current_frame_terminal_count) ) * point_diminish;  
				};

			for( i in 1:10 ) # note point sizes go large to small: 20, 16, 12, 8, 6, 4, 3, 2, 1, 0.8)
				{
				startalpha <- alphas[i];  colorindex <- colorindices[i]; current_offset <- xoffsets[i]
				startcex <- point_cexes[i] * flicker_val;

				colors2<- rgb(colorfunc(colorindex) , maxColorValue=255)
				command_string <- paste( "colors2 <- rgb(colorfunc(" , colorindex , ") , maxColorValue=255); " , 
						"points(" , current_frame_x_limit , "," , point_y , ", cex= " , startcex , 
						", col=alpha( colors2 ," ,  startalpha , ") , pch=16)  # FRAME_EDGE," , i , sep = "")
				write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")
				}
			}
		####################################################################################################################


		} # for(row in 1:number_rows) # for current frame, look through entire tree and see what fits limit.

#	onlythese <- sample(splits_in_current_bin, 5); sample_rows[onlythese] <- 1


	# error encountered. massive number of y's recorded, but all the same value. check for this here.
	# print(c("x1 & x2 both <= current_frame_x_limit:",case1,"both x's outside (RHS) frame:",case2,"case3",case3,"case4",case4,"case5",case5))
	# print(c("current_frame_points_printed", current_frame_points_printed ))



	if(length(unique(frame_speciations)) <= 1 & length(frame_speciations) >= 1000){print ("ERROR 2974");break()}
	# bug occuring in last frame only ... due to coopting variable not assigned in some cases.



	##########################################################################################
	# for FRAME PLOT. 
	# Perl script, sort plotting commands into order which makes image nice looking,
	# current frame x, and current age, are input seperatly.
	frame_age <- round(  (1-current_frame_x_limit) * root_age , digits = 0); 
	if(frame_age <= 1){frame_age = 0};
	if(frame_number >= total_number_frames){frame_age = 0};

	frame_age <- format(frame_age, scientific=F);

	command_string <- paste( "perl sort_commands_for_plotting_frame.pl commands_for_current_frame commands_sorted " , current_filename , " " , current_frame_x_limit , " " , frame_age , " ", pan_speed, sep = "")
	print(c("running:",command_string))
	system(command_string)
	print("plotting in R")
	system( paste( "rm " , current_filename , sep = "") )
	system("R < commands_sorted --vanilla --slave")
	##########################################################################################


	##########################################################################################
	# for AUDIO, 
	# y positions of all bifurcations in current x column have been recorded, if there are not many, record all. if many, subset.	
	frame_speciations_count <- frame_speciations_count - 1;	print(c("frame_speciations_count" , frame_speciations_count))
#	if(frame_speciations_count <= max_notes)
#		{
#		audify1[ 1:frame_speciations_count , frame_number ] <- frame_speciations # insert list of y's
#		audify2[ 1:frame_speciations_count , frame_number ] <- length_after_speciation # insert length after birurcation
#		}else{
#		sample_indices <- sample(1:length(frame_speciations), max_notes, prob = 1-length_after_speciation, replace = FALSE)
#		audify1[ , frame_number ] <- frame_speciations[sample_indices]
#		audify2[ , frame_number ] <- length_after_speciation[sample_indices]
#		}
	##########################################################################################


	#stop();
	} # for (frame_number in total_number_frames)

##################################################################################################################################



