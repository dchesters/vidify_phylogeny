



###########################################################################################################################################
#
#
#
#	audify_phylogeny.R, part of the 'vidify_phylogeny' tool for depicting a phylogeny as audio and video.
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
#
###########################################################################################################################################
#
#
#
# 	There is 1 input file, named: 'vidify_table'
# 	script will delete tracks from previous run, if present (audio_track*)
# 	to run:
# 	R < audify_phylogeny.R --vanilla --slave
#
#
#
# 	CHANGE LOG
# 
# 	2021-10-04 (v 1.01): option to tonify all nodes, not just main ones (appropriate only for very large trees)
# 
# 
# 
# 
###########################################################################################################################################


library(scales)
library(tuneR)


 # Parameters 

tonify_all_nodes <- 1 # large tree use 0; smaller tree use 1

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

###########################################################################################################################################


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
quitebig <- quantile(vid.table1[  , 5],  probs = 0.995)
bignodes <- (1:length(vid.table1[  , 5]))[ vid.table1[  , 5] >= quitebig ]
vid.table1[ bignodes , 5] <- quitebig;vid.table1[  , 5] <- vid.table1[  , 5] / quitebig
aud.table1[ bignodes , 5] <- quitebig;aud.table1[  , 5] <- aud.table1[  , 5] / quitebig

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




bin_size <- 1 / total_number_frames # 0.01
number_rows <- length(vid.table1[ , 1 ] );print(number_rows)




########################################################################################################################################

########################################################################################################################################




# distortion when trying to play so many tones at once, sample ...
# current sampling at random, though would be better to remove tonally similar.
# 2021-08-14: distortion arose from loudness with lots of notes, can be adressed, so following might not be neccessary.


temp1 <- aud.table1[  , 1] +  - min(aud.table1[  , c(1,3)]) ; temp1b <- temp1 / max(aud.table1[  , c(1,3)])
temp2 <- aud.table1[  , 3] +  - min(aud.table1[  , c(1,3)]) ; temp2b <- temp2 / max(aud.table1[  , c(1,3)])
temp3 <- aud.table1[  , 2] +  - min(aud.table1[  , c(2,4)]) ; temp3b <- temp3 / max(aud.table1[  , c(2,4)])
temp4 <- aud.table1[  , 4] +  - min(aud.table1[  , c(2,4)]) ; temp4b <- temp4 / max(aud.table1[  , c(2,4)])
aud.table1[  , 1] <- temp1b; aud.table1[  , 3] <- temp2b; aud.table1[  , 2] <- temp3b; aud.table1[  , 4] <- temp4b

sample_amount <- 10
bin_multiplier <- 1
sample_rows <- rep(0, length.out=length(aud.table1[ , 1]))

# for (frame_number in 1:(total_number_frames+finishing_frames)) # total_number_frames)
#		just up to total_number_frames? remaining should have no bracnhes
for (frame_number in 1:total_number_frames)
	{

	if(frame_number < total_number_frames)
		{
		current_frame_x_limit <- bin_size * frame_number
		}else{
		current_frame_x_limit <- bin_size * (total_number_frames-1)
		}

	print(c(frame_number  , " out of " , total_number_frames  ))
	print (c((current_frame_x_limit-bin_size) , "to" , current_frame_x_limit,"total rm:" , sum(sample_rows)))

	splits_in_current_bin <-0; splits_in_current_count <- 1; 
	for(row in 1:number_rows) # for current frame, look through entire tree and see what fits limit.
		{
		current_x1 <- aud.table1[ row, 1 ];current_x2 <- aud.table1[ row , 3 ];	print_line <- 0; print_point <- 0; point_y <- NA

		# testing just initiation of branch (not where it ends)
		if(current_x1 > (current_frame_x_limit-(bin_size*bin_multiplier)) & current_x1 <= current_frame_x_limit)
			{
			splits_in_current_bin[splits_in_current_count] <- row; splits_in_current_count <- splits_in_current_count+1
			}
		} # for(row in 1:number_rows) # for current frame, look through entire tree and see what fits limit.

	if(splits_in_current_count >= sample_amount)
		{
		remove_count <- splits_in_current_count - sample_amount
		print("sampling current x slice")
		print(c("removing" ,remove_count))
		rmthese <- sample(splits_in_current_bin, remove_count, replace=F); 
		sample_rows[rmthese] <- 1
		}

	}
sum(sample_rows)
# 75370
length(sample_rows)
# 84416




########################################################################################################################################

########################################################################################################################################




	# first working version of several attempts at methods for audifying,
	# pick out branches, for each make wave, silence, decaying tone, silence,
	# then combine waves 

# tuneR::sine
#	freq: The frequency (in Hertz) to be generated.
#	duration: Duration of the ‘Wave’ in ‘xunit’.
#	from: Starting value of the ‘Wave’ in ‘xunit’.

# complete:
 branches.table <- aud.table1
# or sampled:
# branches.table <- aud.table1[ (1:length(sample_rows))[sample_rows == 0] , ]



bin_size <- 1 / (total_number_frames+finishing_frames) # 
number_rows <- length(branches.table[ , 1 ] );print(number_rows)

newmin <- 100; newmax <- 900
branches.table[  , 2] <- (branches.table[  , 2] - min(branches.table[  , c(2,4)]))/(max(branches.table[  , c(2,4)])-min(branches.table[  , c(2,4)])) * (newmax - newmin) + newmin
branches.table[  , 4] <- (branches.table[  , 4] - min(branches.table[  , c(2,4)]))/(max(branches.table[  , c(2,4)])-min(branches.table[  , c(2,4)])) * (newmax - newmin) + newmin

# also standardize x values (corresponds to audio length)
# [is this repeated?]
temp1 <- branches.table[  , 1] +  - min(branches.table[  , c(1,3)]) ; temp1b <- temp1 / max(branches.table[  , c(1,3)])
temp2 <- branches.table[  , 3] +  - min(branches.table[  , c(1,3)]) ; temp2b <- temp2 / max(branches.table[  , c(1,3)])
branches.table[  , 1] <- temp1b; branches.table[  , 3] <- temp2b; 


# test how many

long_tones <- 0;short_tones <- 0;max_tone_length <-0

for(row in 1:number_rows) 
	{
	current_x1 <- branches.table[ row, 1 ];current_x2 <- branches.table[ row , 3 ];
	current_y1 <- branches.table[ row, 2 ];current_y2 <- branches.table[ row , 4 ];
	tone_length <- current_x2 - current_x1;
	if(max_tone_length <= tone_length){max_tone_length <- tone_length}
	if(current_x1 < current_x2 & current_x1 > 0)
		{
		if(tone_length >= 0.2)
			{
			long_tones <- long_tones+1
			}
		if(tone_length < 0.1 & tone_length >= 0.05)
			{
			short_tones <- short_tones+1
			};
		}
	}

print(c("long_tones" , long_tones))
print(c("short_tones" , short_tones))
# "long_tones" "104"       
 print(c("short_tones" , short_tones))
# "short_tones" "235"    
print(c("max_tone_length" , max_tone_length))
# "max_tone_length"   "0.787508646425822"

sampling.rate <- 3000
bits <- 32
start_adjustment <- 0 # audio and video dont quite match

system("rm audio_track*")

# length 1:22
long_tones <- 0;short_tones <- 0;max_tone_length <-0

# early loops will have some very long bracnhes which are slower. will speed up later.
file_count_A <- 1;file_count_B <- 1;
	###################################################################################################
for(row in 1:number_rows) # no rows of branches.table.
	{
	# row <- 4
	current_x1 <- branches.table[ row, 1 ];current_x2 <- branches.table[ row , 3 ]; nodesize <- branches.table[ row , 5 ];
	current_y1 <- branches.table[ row, 2 ];current_y2 <- branches.table[ row , 4 ];
	print(c("row" , row  , "of" , number_rows))

	tone_length <- current_x2 - current_x1

	# the more tones combined, the loss of quality of individual tones,
	# thus do 2 major files, small number of longer tones that will be noticed,
	# and large number of short tones that will mostly form the mass, and for which quality decrease is unimportant.

	if(current_x1 < current_x2 & current_x1 > 0)
		{

		test1<-0;	if(tone_length >= 0.2) 		{test1 <- 1};
				if(current_x1 <= 0.33) 		{test1 <- 1};
				if(nodesize >= 0.1) 		{test1 <- 1};
				if(tonify_all_nodes == 1)	{test1 <- 1};

		##########################################################################################################
		if( test1 == 1 )
			{
			print(c("key branch:" ,"partA",current_x1, current_x2, "time_to_start_tone:" , (current_x1*lengthener)+start_adjustment))

			# pre-tone silence:
			w <- tuneR::silence(duration = (current_x1*lengthener)+start_adjustment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)

			# tone: default is mono ....  stereo=TRUE
			wobj1 <- tuneR::sine( current_y1 , duration=(tone_length*lengthener)+postsegment, from = 0, xunit = c("samples", "time")[2], 
					bit=bits, samp.rate=sampling.rate, stereo=TRUE);

			# decay:
	 		len <- length(wobj1@left); current_divis <- 1; decay_factor <- 0.005 # default long = 0.005; short = 0.1
			for (j in 1:len)
				{
				wobj1@left[j] <- wobj1@left[j] / current_divis;wobj1@right[j] <- wobj1@right[j] / current_divis; 
				current_divis <- current_divis + decay_factor
				}
			wobj <- wobj1;  wobj <- prepComb(wobj, where = "start");  w <- tuneR::bind( w, wobj); # tuneR::bind does concatenation

			# post-tone silence:
		#	if(1-current_x2 <= 0.000001)
		#		{
		#		print ("branch finishes at end, no need for post silence")
		#		}else{
				post_silence <- tuneR::silence(duration = ((1-current_x2)*lengthener)+postsegment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
				# Warning messages:1: In min(x) : no non-missing arguments to min; returning Inf
				# error message if small value is input to duration (branch end x effecivly at end)
				w <- tuneR::bind( w, post_silence);
		#		}

				
			#	track_name <- paste( "audio_track_A." , row  , sep = "")
				track_name <- paste( "audio_track_A." , file_count_A  , sep = ""); file_count_A <- file_count_A+1
				w <- tuneR::normalize(w,unit = "32") # this can change volume if all quite

				# scale volume after normalization:
				nodesize <- nodesize*3; if(nodesize >= 1){nodesize <- 1};w <- w * nodesize

				writeWave(w,track_name);long_tones <- long_tones+1

			}else{
		print(c(current_x1, current_x2))

		skip_little_nodes <- 1
		if(skip_little_nodes == 0)
		{
		# pre-tone silence:
		w <- tuneR::silence(duration = (current_x1*lengthener)+start_adjustment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
		# tone:
		wobj1 <- tuneR::sine( current_y1 , duration=(tone_length*lengthener)+postsegment, from = 0, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE);
		# decay:
 		len <- length(wobj1@left); current_divis <- 1; 
		for (j in 1:len)
			{wobj1@left[j] <- wobj1@left[j] / current_divis;wobj1@right[j] <- wobj1@right[j] / current_divis;current_divis <- current_divis + 0.005};
		wobj <- wobj1;  wobj <- prepComb(wobj, where = "start");  
		w <- tuneR::bind( w, wobj); # tuneR::bind does concatenation
		# post-tone silence:
		post_silence <- tuneR::silence(duration = ((1-current_x2)*lengthener)+postsegment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
		w <- tuneR::bind( w, post_silence);
		# if previous command didnt run: Error in FUN(X[[i]], ...) :   One 'Wave' object is mono, the other one stereo
		track_name <- paste( "audio_track_B." , file_count_B  , sep = ""); file_count_B <- file_count_B+1
		w <- tuneR::normalize(w,unit = "32") # this can change volume if all quite
		writeWave(w,track_name)
		}

			}; # short branches



		}




	} # for(row in 1:number_rows) #  look through entire tree
	###################################################################################################

print(c("long_tones" , long_tones))

system("ls -l audio_track_A.* | wc -l")

##########################################################











