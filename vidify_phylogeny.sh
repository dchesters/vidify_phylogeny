
#	
# 'vidify_phylogeny', Linux tools for depicting a phylogeny as audio and video.
#		  
# Copyright (C) 2021  Douglas Chesters
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# contact address dc0357548934@live.co.uk
#
#
#############################################################################################################
#
#
# Dependencies:
#	ffmpeg
#	audacity
# 	process_newick.pl
# 	plot_frames.R
#	audify_phylogeny.R
# 	sort_commands_for_plotting_frame.pl
# 	R packages plotrix, scales, tuneR
#
#
#############################################################################################################


  # Read phylogeny, interpret each node as X and Y position, then tabulate.

# Tested on process_newick.pl v2.26
# Following settings should be used in script:
#   $highlight_tips_from_file = 0;
#   $reduce_very_long_branches_to_this = 1000.3;
#   $print_internal_node_labels_on_rectangular_tree = 1;
intree=filename.nwk
rm vidify_table
perl process_newick.pl -intree $intree -outprefix $intree -ncbi_tax_node 6960 -plot_tree -plot_trees_with_branchlengths


#############################################################################################################


  # Read table of nodes and plot sequential frames in root to tip direction.

# Following command requires Perl script (sort_commands_for_plotting_frame.pl) in working directory.
# Needs to have installed R libraries: plotrix, scales, tuneR
# Input: vidify_table
# It will remove if present, all previous video frames (video_frame.*)
R < plot_frames.R --vanilla --slave

# Output are loads of video frames, each jpeg image. check these are present with following command. 
# Under default settings should be about 1621 frames (for 70 second video)
ls video_frame.*.jpeg | wc -l


#############################################################################################################


  # make video from frames

# Uses default framerate of 25 (use -framerate switch to change)
ffmpeg -i video_frame.%d.jpeg Ephemeroptera_Video.mpeg


#############################################################################################################


  # Convert each bifurcation of phylogeny, to a tone (silence up to x position followed by diminishing tone), 
  #  and save to indivdual files.

# Two R libraries should be installed: scales, tuneR
# There is 1 input file, named: 'vidify_table'
# script will delete tracks from previous run, if present (audio_track*)
R < audify_phylogeny.R --vanilla --slave

# How many tones (each in seperate file), this number depends on number of nodes in phylogeny:
ls audio_track_A.* | wc -l


#############################################################################################################


  # Combine tones (currently each in seperate file) to single track.

# Manual steps required here. 
#  (A couple of command line options i tried gave terrible results. Thus opted for Audacity)
#  Open AUDACITY and do the following:
#	Import audio tracks (no more than abuot 700 at a time), 
#		for each set imported: Select all, Tracks/mix/mix and render.
#	When all imported: 
#	Select track, effect, normalize, normalize peak amplitude to -10 db, click OK.
#	Effect, base and treble, manually increase bass, click apply, close.
# 	File, Export, export as WAV (filename audio.wav).


#############################################################################################################


  # Combine audio and video using linux command line tool.

ffmpeg -i Ephemeroptera_Video.mpeg -i audio.wav -c:v copy -c:a aac AudioVideo.mp4


#############################################################################################################


