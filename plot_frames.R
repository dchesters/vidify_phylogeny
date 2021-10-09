
###################################################################################################################################
#
# 
#	plot_frames.R, part of the 'vidify_phylogeny' tool for depicting a phylogeny as audio and video.
#		  
#    	Copyright (C) 2021  Douglas Chesters
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
#	
#	CHANGE LOG
#	 
#	2021-SEP-11 (version 1.01)
#	2021-OCT-03 (version 1.02): checks for dependent Perl script in working directory.
#	 
#
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
table2 <- read.table("vidify_table", header=F, row.names = 1, sep = "\t")
head(table2)


colorfunc = colorRamp(c(
	"darkblue" ,"darkblue" ,"darkblue" ,  "red" ,"red" ,"red" ,  "yellow","yellow",  "white","white"
	))


# each have column 9 which is count of terminals derived from each brnach
 # VIDEO:
vid.table1 <- table2[ , c(1,2,3,4,9)]
 # AUDIO:
aud.table1 <- table2[ , c(5,6,7,8,9)]

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
temp1 <- vid.table1[  , 1] +  - min(vid.table1[  , c(1,3)]) ; temp1b <- temp1 / max(vid.table1[  , c(1,3)])
temp2 <- vid.table1[  , 3] +  - min(vid.table1[  , c(1,3)]) ; temp2b <- temp2 / max(vid.table1[  , c(1,3)])
temp3 <- vid.table1[  , 2] +  - min(vid.table1[  , c(2,4)]) ; temp3b  <- temp3 / max(vid.table1[  , c(2,4)])
temp4 <- vid.table1[  , 4] +  - min(vid.table1[  , c(2,4)]) ; temp4b <- temp4 / max(vid.table1[  , c(2,4)])
vid.table1[  , 1] <- temp1b; vid.table1[  , 3] <- temp2b; vid.table1[  , 2] <- temp3b; vid.table1[  , 4] <- temp4b


 # Parameters 

# NOTE, 2000 is very slow, test on much less.
# one minute of tree traversal plus 10 seconds out
# 25 fps * 60 seconds is 1500
total_number_frames <- 1500 # default 1500
# 25 * 5 seconds
finishing_frames <- 125 # default; 125 after tree terminals reached, video should not immediatly cut off, add a few frames, stationary.

# equivelent for audio:
# this will give audio length of 1 second, so:
lengthener <- 60 # default 60
postsegment <- 5 # default 5
bin_size <- 1 / total_number_frames # 0.01
number_rows <- length(vid.table1[ , 1 ] );print(number_rows)





# for video only


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
branchlwd.outer <- 4;branchlwd.inner <- 1;point_diminish <- 1



##################################################################

print(total_number_frames)

system("rm video_frame_commands*");system("rm video_frame.*")

for (frame_number in 1:(total_number_frames+finishing_frames)) # total_number_frames)
	{
#	frame_number <- 50

	if(frame_number < total_number_frames)
		{
		current_frame_x_limit <- bin_size * frame_number
	#	branchlength <- x2-x1;
		}else{
		############################################################################################
		# end phase, 
		current_frame_x_limit <- bin_size * (total_number_frames-1)
	#	branchlength <- (bin_size * frame_number)-x1;

		# glow of points denoting branch edge will fade in end phase		
	#	for( i in 1:10 ){point_cexes[i] <- point_cexes[i]*phaseout}
		point_diminish <- point_diminish * 0.95;
		branchlwd.outer <- branchlwd.outer * 0.92; branchlwd.inner <- branchlwd.inner *0.92
		############################################################################################
		}
	current_filename <- paste ( "video_frame." , frame_number , ".jpeg" , sep = "")
	system("rm commands_for_current_frame");
	print(c("frame_number" , frame_number, "current_frame_x_limit" , current_frame_x_limit))
	frame_speciations <- 0; length_after_speciation <- 0; frame_speciations_count <- 1;
	case1<-0;case2<-0;case3<-0;case4<-0;case5<-0;
	splits_in_current_bin <-0; splits_in_current_count <- 1; 

	for(row in 1:number_rows) # for current frame, look through entire tree and see what fits limit.
		{
		current_x1 <- vid.table1[ row, 1 ];current_x2 <- vid.table1[ row , 3 ];	print_line <- 0; print_point <- 0; point_y <- NA

		

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


		####################################################################################################################
		if(print_line == 1)
			{
			#   BRANCH

			startalpha <- 0.01;startlwd <- 30; startcolor <- 0.1;lwd_shift_factor <- 4 # default=2
			fancy_branch <- 2
			if(fancy_branch == 1)
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
				}else{
				#	note very thin branches will appear to flicker in the video
					# 0.5 = red; 0.8=yellow; 0.95=white
					command_string <- paste( "colors2 <- rgb(colorfunc(" , 0.5 , ") , maxColorValue=255);segments(", 
							 x1 , "," , y1 , "," , x2 , "," , y2 , ",alpha( colors2 , ", 0.999 , "), lwd = " ,
							 1.5 , ") # BRANCH," , 1 , sep = "")
				#	write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")
					command_string <- paste( "colors2 <- rgb(colorfunc(" , 0.8 , ") , maxColorValue=255);segments(", 
							 x1 , "," , y1 , "," , x2 , "," , y2 , ",alpha( colors2 , ", 0.999 , "), lwd = " ,
							 branchlwd.outer , ") # BRANCH," , 2 , sep = "")
					write(command_string , file = "commands_for_current_frame" , ncolumns = 1,  append = T, sep = "\t")
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
			x<-rlnorm(101);x<-x/max(x);x<-sort(x, decreasing=T); flicker_val <- x[round(branchlength*100)+1]
			# also, gets a bit ugly at Right hand side (reduces to colored bars), so diminish cex along frame sequence
			# flicker_val <- flicker_val * (1-(current_frame_x_limit/pointcex_x_diminish_factor)); 
			# maybe more informative way to do this is scale to number of terminals derived fro node:
			flicker_val <- (flicker_val * vid.table1[ row , 5]) * point_diminish; 

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
	print(c("x1 & x2 both <= current_frame_x_limit:",case1,"both x's outside (RHS) frame:",case2,"case3",case3,"case4",case4,"case5",case5))
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

	command_string <- paste( "perl sort_commands_for_plotting_frame.pl commands_for_current_frame commands_sorted " , current_filename , " " , current_frame_x_limit , " " , frame_age , sep = "")
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



	} # for (frame_number in total_number_frames)

##################################################################################################################################



