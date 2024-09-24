



###########################################################################################################################################
#
#
#
#	audify_phylogeny.R, tool for depicting a phylogeny as audio.
#		  
#    	Copyright (C) 2021-2024  Douglas Chesters
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
# 	2021-10-04 (v 1.01): 	Option to tonify all nodes, not just main ones (appropriate only for very large trees)
# 	2024-08-12 (v 2.01): 	Implementing pitchshifting tones for more information content.
#				Bugfix changing Y to audio frequency range
# 	2024-08-24 (v 2.02):	Bugfix stop crashing on full length branches
# 	2024-08-31 (v 2.03): 	Bugfix tone subsampling according to clade size.
# 	2024-09-04 (v 2.04):	To prevent abrupt finish of tones at tree terminus, 
#				tones persist beyond terminus though decay fast
# 	2024-09-12 (v 2.05):	Automatically adjust parameters according to number of branches input.
# 
# 
# 
# 
# 
# 
# 
# 	NOTEs
# 	how might stereo be used, 
# 
# 
# 
# 
###########################################################################################################################################


library(scales)
library(tuneR)


 # Parameters 

# 25 * 5 seconds
finishing_frames 	<- 125 	# default; 125 after tree terminals reached, video should not immediatly cut off, add a few frames, stationary.
postsegment 		<- 5 	# default 5

pitchshift 		<- 1 	# default==0, 1== pitch ascends or decends according to child branch delta Y
				# uses different input table




###########################################################################################################################################


# video and audio columns

colorfunc = colorRamp(c(
	"darkblue" ,"darkblue" ,"darkblue" ,  "red" ,"red" ,"red" ,  "yellow","yellow",  "white","white"
	))


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
	};

head(table2)

number_branches <- length(table2[,1]); # about 2x number of terminals.
print(c("number_branches",number_branches))





#################################################################################################


 # optimize parameters according to size of the phylogeny


 # if < 800 tones, decay_factor <- 0.005 ; decay_factor.terminals <- 0.01
 # 0.2 makes too short tones that dont convey frequency
 # decay_factor of 0.05 was too sharp for 2 minute 
tonelength_cutoff <- 0.01;x1_cutoff<-0.50;nodesize_cutoff<-15 # these not used for smaller trees, but still need defining

if(number_branches < 160)
{
tonify_all_nodes<-1;track_length<-30
decay_factor<-0.002;decay_factor.terminals<-0.005;print("tree < 160 branches")
}else{

if(number_branches < 600)
{
tonify_all_nodes 	<- 1 	# large tree use 0; smaller tree use 1
track_length 		<- 60	# default is 60 seconds, for terminal-rich trees, make longer
decay_factor 		<- 0.002 	# default = 0.005; 	shorten due to many splits = 0.05
decay_factor.terminals 	<- 0.005 	# default 0.01; 	shorten due to many splits = 0.05
print("tree < 600 branches")
}else{

if(number_branches < 1200)
{
tonify_all_nodes <- 1;track_length <- 60;decay_factor <- 0.002;decay_factor.terminals <- 0.01;print("tree < 1200 branches")
}else{

if(number_branches < 1600)
{
tonify_all_nodes<-1;track_length<-60;decay_factor<-0.002;decay_factor.terminals<-0.01;print("tree < 1600 branches")
}else{

if(number_branches < 2000)
{
tonify_all_nodes<-1;track_length<-60;decay_factor<-0.005;decay_factor.terminals<-0.01;print("tree < 2000 branches")
}else{

if(number_branches < 2400) #
{
tonify_all_nodes<-0;track_length<-60;decay_factor<-0.01;decay_factor.terminals<-0.02;print("tree < 2400 branches")
tonelength_cutoff <- 0.01;x1_cutoff<-0.50;nodesize_cutoff<-15 # these where tonify_all_nodes==0
}else{

if(number_branches < 3200)
{
tonify_all_nodes<-0;track_length<-60;decay_factor<-0.01;decay_factor.terminals<-0.02;print("tree < 3200 branches")
tonelength_cutoff <- 0.01;x1_cutoff<-0.50;nodesize_cutoff<-15
}else{

if(number_branches < 8000)
{
tonify_all_nodes<-0;track_length<-120;decay_factor<-0.05;decay_factor.terminals<-0.05;print("tree < 8000 branches")
tonelength_cutoff <- 0.03;x1_cutoff<-0.35;nodesize_cutoff<-50
}else{

if(number_branches < 10000)
{
tonify_all_nodes<-0;track_length<-120;decay_factor<-0.05;decay_factor.terminals<-0.05;print("tree < 10000 branches")
tonelength_cutoff <- 0.04;x1_cutoff<-0.35;nodesize_cutoff<-50
}else{

if(number_branches >= 10000)
{
tonify_all_nodes<-0;track_length<-120;decay_factor<-0.06;decay_factor.terminals<-0.07;print("tree > 10000 branches")
tonelength_cutoff <- 0.08;x1_cutoff<-0.32;nodesize_cutoff<-70
}

}}}}}}}}}





# NOTE, 2000 is very slow, test on much less.
# one minute of tree traversal plus 10 seconds out
# 25 fps * 60 seconds is 1500
total_number_frames <- track_length*25; # default 1500
# equivelent for audio:
# this will give audio length of 1 second, so:
lengthener <- track_length # default 60; 1500 terminals then 120; 3000 terminals then 180



# process both tables

root_age <- max(vid.table1[  , c(1,3)]); # in millions years
print(c("root_age", root_age))

max_number_terminals <- max(vid.table1[  , 5])
quitebig <- quantile(vid.table1[  , 5],  probs = 0.995)
bignodes <- (1:length(vid.table1[  , 5]))[ vid.table1[  , 5] >= quitebig ]
vid.table1[ bignodes , 5] <- quitebig;vid.table1[  , 5] <- vid.table1[  , 5] / quitebig
# aud.table1[ bignodes , 5] <- quitebig;aud.table1[  , 5] <- aud.table1[  , 5] / quitebig

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






temp1 <- aud.table1[  , 1] +  - min(aud.table1[  , c(1,3)]) ; temp1b <- temp1 / max(aud.table1[  , c(1,3)])
temp2 <- aud.table1[  , 3] +  - min(aud.table1[  , c(1,3)]) ; temp2b <- temp2 / max(aud.table1[  , c(1,3)])
temp3 <- aud.table1[  , 2] +  - min(aud.table1[  , c(2,4)]) ; temp3b <- temp3 / max(aud.table1[  , c(2,4)])
temp4 <- aud.table1[  , 4] +  - min(aud.table1[  , c(2,4)]) ; temp4b <- temp4 / max(aud.table1[  , c(2,4)])
aud.table1[  , 1] <- temp1b; aud.table1[  , 3] <- temp2b; aud.table1[  , 2] <- temp3b; aud.table1[  , 4] <- temp4b


# distortion when trying to play so many tones at once, sample ...
# current sampling at random, though would be better to remove tonally similar.
# 2021-08-14: distortion arose from loudness with lots of notes, can be adressed, so following might not be neccessary.

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

	print(c(frame_number  , " out of " , total_number_frames , "frames" ))
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

# set audio frequency range here, average about 440
oldmin <- min(branches.table[  , c(2,4)]); oldmax <- max(branches.table[  , c(2,4)])
newmin <- 200; newmax <- 900
branches.table[  , 2] <- (branches.table[  , 2] - oldmin)  /  (oldmax-oldmin) * (newmax - newmin) + newmin
branches.table[  , 4] <- (branches.table[  , 4] - oldmin)  /  (oldmax-oldmin) * (newmax - newmin) + newmin

# also standardize x values (corresponds to audio length)
# [is this repeated?]
temp1 <- branches.table[  , 1] +  - min(branches.table[  , c(1,3)]) ; temp1b <- temp1 / max(branches.table[  , c(1,3)])
temp2 <- branches.table[  , 3] +  - min(branches.table[  , c(1,3)]) ; temp2b <- temp2 / max(branches.table[  , c(1,3)])
branches.table[  , 1] <- temp1b; branches.table[  , 3] <- temp2b; 

# to assist interpretation, plot lines corresponding to all tones 
pdf(file = "tones_plotted.pdf", 6, 6);
plot(c(0,1),c(100,900), type = "n");



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

head(branches.table)
#           V2       V3         V4       V5 V6
# 2  0.00000000 500.0000 0.03751352 100.2141 14
# 4  0.03751352 873.1959 0.23624106 100.2210  2

# stop()


# early loops will have some very long bracnhes which are slower. will speed up later.
file_count_A <- 1;file_count_B <- 1;file_count_C <- 1;
	###################################################################################################
for(row in 1:number_rows) #number_rows) # no rows of branches.table.
	{
	# row <- 4
	current_x1 <- branches.table[ row, 1 ];current_x2 <- branches.table[ row , 3 ]; nodesize <- branches.table[ row , 5 ];
	current_y1 <- branches.table[ row, 2 ];current_y2 <- branches.table[ row , 4 ];
	print(c("row" , row  , "of" , number_rows, "Xs:",current_x1,current_x2))

	tone_length <- current_x2 - current_x1

	# the more tones combined, the loss of quality of individual tones,
	# thus do 2 major files, small number of longer tones that will be noticed,
	# and large number of short tones that will mostly form the mass, and for which quality decrease is unimportant.

	

	if(current_x1 < current_x2 & current_x1 >= 0)
		{
		
		test1<-0;	if(tone_length 		>= tonelength_cutoff) 		{test1 <- 1}; # 0.2
				if(current_x1 		<= x1_cutoff) 		{test1 <- 1}; # 0.33

				# do need to consider clade size, these can be short branches 
				# but due to significace can be noticable ommisions
				# previous default (proportion) didnt make sense, these look like descendent terminal counts
				if(nodesize 		>= nodesize_cutoff) 	{test1 <- 1}; #
				if(tonify_all_nodes 	== 1)			{test1 <- 1};

		print(c(tone_length,current_x1,nodesize,tonify_all_nodes,test1))
		time_to_start_tone <- (current_x1*lengthener)+start_adjustment

		##########################################################################################################
		if( test1 == 1 )
			{
			print(c("key branch:" ,"partA",current_x1, current_x2, "time_to_start_tone:" , (current_x1*lengthener)+start_adjustment))

			if(time_to_start_tone <= 0.00000001) # if value 0 here, will crash when very tiny start branch encountered (Struck2023Annelida1)
				{}else{
				# pre-tone silence:
				silence_lhs <- tuneR::silence(duration = (current_x1*lengthener)+start_adjustment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
				silence_lhs <- prepComb(silence_lhs, where = "end");
				}




			if(current_x1 >= 0.92 & current_x2 > 0.9999)
				{

				###########################################################################
					# near end of X traversal, dont bother with pitchshifting,
					# also final tones will go beyond terminus to avoid abrupt ends.
					# thus, just need 2 items, pre-tone silence and strait tone
				wave6 <- tuneR::sine( current_y1 , duration=((1-current_x1)*lengthener)+postsegment, from = 0, xunit = c("samples", "time")[2], 
					bit=bits, samp.rate=sampling.rate, stereo=TRUE);
		 		len <- length(wave6@left); current_divis <- 1; 
				for (j in 1:len)
					{
					wave6@left[j] <- wave6@left[j] / current_divis;wave6@right[j] <- wave6@right[j] / current_divis; 
					current_divis <- current_divis + decay_factor.terminals
					}
				track_name <- paste( "pitchshifting_terminal." , file_count_C  , sep = ""); file_count_C <- file_count_C+1
				wave6 <- prepComb(wave6, where = "start");  
				wave6 <- tuneR::bind( silence_lhs, wave6);
				wave6 <- tuneR::normalize(wave6,unit = "32") # this can change volume if all quite
				# scale volume after normalization:
			#	nodesize <- nodesize*3; if(nodesize >= 1){nodesize <- 1};wave_static <- wave_static * nodesize
				writeWave(wave6,track_name);
				###########################################################################


				}else{



			# print("line 350.")
			# tone: default is mono ....  stereo=TRUE
			wobj1 <- tuneR::sine( current_y1 , duration=(tone_length*lengthener)+postsegment, from = 0, xunit = c("samples", "time")[2], 
					bit=bits, samp.rate=sampling.rate, stereo=TRUE);



			#######################################################################################################################
			
			# 20240812: try different tone style, pitch shift between start and end of branch according to child node delta Y 			
			# an argument for stepping pitchshift is audio cue for matching child pairs
			# notes , try working with only one channel
			# s1 <- sndObj@left

			if(pitchshift == 1)
				{
				# key columns are  1)sum_branchlength  2)y1_proportion  3)assign_new_x  4)y2_proportion
				print(c("tone_length",tone_length))
				pitchshift_segment_length <- 0.004
				if(tone_length > (pitchshift_segment_length*2))
					{
					print(c("current_x1",current_x1,"current_y1",current_y1,"current_x2",current_x2,"current_y2",current_y2))	
					shift_segments <- floor(tone_length/pitchshift_segment_length)
					# print("line 372.")

					wave1 <- tuneR::sine( current_y1 , duration=(pitchshift_segment_length*lengthener), from = 0, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE);
					track_duration <- pitchshift_segment_length;
					for (segment in 1:(shift_segments-1))
						{
						segment_Xl <- current_x1+(segment*pitchshift_segment_length);
						m <- (current_y2 - current_y1) / (current_x2 - current_x1) # slope
						b <- current_y1 - m * current_x1; 			# y-intercept using one of the points
						segment_Y <- m * segment_Xl + b; 			# y-coordinate for given X, using the line equation
						# print(c("segment_Y",segment_Y))
						wave2 <- tuneR::sine( segment_Y , duration=(pitchshift_segment_length*lengthener), from = 0, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE);
						wave1 <- prepComb(wave1, where = "end");wave2 <- prepComb(wave2, where = "start");
						wave1 <- tuneR::bind( wave1,wave2);track_duration <- track_duration+pitchshift_segment_length;
						segments(segment_Xl,segment_Y,segment_Xl+pitchshift_segment_length,segment_Y, col="black",lwd=0.8);
						}
					print(c("track_duration",track_duration))	

					# attach final part, which is shorter than segment threshold
					segment_Xl <- current_x1+(shift_segments*pitchshift_segment_length);
					m <- (current_y2 - current_y1) / (current_x2 - current_x1) # slope
					b <- current_y1 - m * current_x1; 			# y-intercept using one of the points
					segment_Y <- m * segment_Xl + b; 			# y-coordinate for given X, using the line equation
					# print("line 395.");

					if((tone_length-track_duration)<=0.0000001)
						{
						# print(c(tone_length, track_duration))
						# print("skip final segment");
						}else{
						# print(c(tone_length, track_duration))
						# print("line 406...");
						wave2 <- tuneR::sine( segment_Y , duration=((tone_length-track_duration)*lengthener), from = 0, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE);
						wave1 <- prepComb(wave1, where = "end");wave2 <- prepComb(wave2, where = "start");
						wave1 <- tuneR::bind( wave1,wave2);
						}

			 		len <- length(wave1@left); current_divis <- 1; 
					for (j in 1:len)
						{
						wave1@left[j] <- wave1@left[j] / current_divis;wave1@right[j] <- wave1@right[j] / current_divis; 
						current_divis <- current_divis + decay_factor
						}

					# print("line 408.");
					# silence before branch starts
					if(time_to_start_tone <= 0.00000001 )
						{}else{
					wave1 <- prepComb(wave1, where = "start");  
					wave1 <- tuneR::bind( silence_lhs, wave1); # tuneR::bind does concatenation
						}

					# print("line 416.")
					# post-tone silence:
					post_silence <- tuneR::silence(duration = ((1-current_x2)*lengthener)+postsegment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
					# Warning messages:1: In min(x) : no non-missing arguments to min; returning Inf
					# error message if small value is input to duration (branch end x effecivly at end)
					wave1 <- tuneR::bind( wave1, post_silence);

					# write pitchshifting tone for current branch
					wave1 <- tuneR::normalize(wave1,unit = "32") # this can change volume if all quite
					PS_trackname <- paste( "pitchshifting_A." , file_count_A  , sep = "");file_count_A <- file_count_A+1
					writeWave(wave1,PS_trackname)


					}else{
					print("not pitchshifting such a short branch")	
					}

				}else{
			#######################################################################################################################
			# non-shifting tones

			# decay:
	 		len <- length(wobj1@left); current_divis <- 1; 
			for (j in 1:len)
				{
				wobj1@left[j] <- wobj1@left[j] / current_divis;wobj1@right[j] <- wobj1@right[j] / current_divis; 
				current_divis <- current_divis + decay_factor
				}
			wobj <- wobj1;  
			
			if(time_to_start_tone ==0)
				{
				wave_static <- wobj
				}else{
				wobj <- prepComb(wobj, where = "start");  wave_static <- tuneR::bind( silence_lhs, wobj); # tuneR::bind does concatenation
				};

			# post-tone silence:
				post_silence <- tuneR::silence(duration = ((1-current_x2)*lengthener)+postsegment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
				# Warning messages:1: In min(x) : no non-missing arguments to min; returning Inf
				# error message if small value is input to duration (branch end x effecivly at end)
				wave_static <- tuneR::bind( wave_static, post_silence);
				
			#	track_name <- paste( "audio_track_A." , row  , sep = "")
				track_name <- paste( "audio_track_A." , file_count_A  , sep = ""); file_count_A <- file_count_A+1
				wave_static <- tuneR::normalize(wave_static,unit = "32") # this can change volume if all quite

				# scale volume after normalization:
				nodesize <- nodesize*3; if(nodesize >= 1){nodesize <- 1};wave_static <- wave_static * nodesize

				writeWave(wave_static,track_name);long_tones <- long_tones+1


				};




				};










			# stop("stopping. 390");

			}else{ # 		if( test1 == 1 )

			# print(c(current_x1, current_x2))
			skip_little_nodes <- 0
			if(skip_little_nodes == 0)
			{


			# tone [ removed +postsegment from duration ]
			wobj1 <- tuneR::sine( current_y1 , duration=(tone_length*lengthener), from = 0, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE);
			# decay:
	 		len <- length(wobj1@left); current_divis <- 1; 
			for (j in 1:len)
				{wobj1@left[j] <- wobj1@left[j] / current_divis;wobj1@right[j] <- wobj1@right[j] / current_divis;current_divis <- current_divis + 0.005};

			wobj <- wobj1;  

			if(time_to_start_tone ==0)
				{
				w <- wobj
				}else{
				# pre-tone silence:
				pre_tone_silence <- tuneR::silence(duration = (current_x1*lengthener)+start_adjustment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
				wobj <- prepComb(wobj, where = "start");  
				w <- tuneR::bind( pre_tone_silence, wobj); # tuneR::bind does concatenation
				}

			# post-tone silence:
			post_silence <- tuneR::silence(duration = ((1-current_x2)*lengthener)+postsegment, xunit = c("samples", "time")[2], bit=bits, samp.rate=sampling.rate, stereo=TRUE)
			w <- tuneR::bind( w, post_silence);
			# if previous command didnt run: Error in FUN(X[[i]], ...) :   One 'Wave' object is mono, the other one stereo
			track_name <- paste( "audio_track_B." , file_count_B  , sep = ""); file_count_B <- file_count_B+1
			w <- tuneR::normalize(w,unit = "32") # this can change volume if all quite
			writeWave(w,track_name)
			segments(current_x1,current_y1,current_x2,current_y1, col="red",lwd=0.8);
			}

			}; # short branches



		}else{
		print(c("ignoring row",row,current_x1,current_x2))
		}




	} # for(row in 1:number_rows) #  look through entire tree
	###################################################################################################

print(c("long_tones" , long_tones));
dev.off();


system("ls -l audio_track_A.* | wc -l")

##########################################################











