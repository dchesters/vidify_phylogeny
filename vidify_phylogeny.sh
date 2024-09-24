
# 
# Get latest version from github.com/dchesters/vidify_phylogeny
# 
# Assuming Unix-like environment.
# Dependencies:
#   R
#   R libraries plotrix, scales, tuneR
#   Perl scripts process_newick.pl, sort_commands_for_plotting_frame.pl
#   R scripts audify_phylogeny.R, plot_frames.R
#   ffmpeg
#   Audacity
#   a text controlfile named controlfile.txt
#   a text subtitles file named subtitles.srt
# Tree viewers recommended for processing (e.g. format conversion, rooting, outgroup removal),
# I suggest you have installed both Figtree (tree.bio.ed.ac.uk/software/figtree/) and 
# Archaeopteryx (phylosoft.org/archaeopteryx/).
#
# Input is Newich format machine reable phylogeny (Nexus format will need converting to Newick first).
# Tree must contain branch lengths. Node labels best removed (certainly complex format ones).
# Can be either ultrametric or phylogram.
# Tree should be rooted, clade-ordered, best to remove outgroup if practical (if clearly labelled as such).
# There should be no duplicate terminal IDs (which will crash many applications).
# Thus, your Newick string should look something like this:
# (Apis_florea:0.5,((Apis_cerana:0.1,Apis_mellifera:0.1):0.2,(Apis_dorsata:0.2,Apis_laboriosa:0.2):0.1):0.2);

# Read phylogeny, interpret each branch as vector (X and Y positions for start and end), then tabulate.
intree=your_tree.nwk
control=controlfile.txt
rm vidify_table2
perl process_newick.pl -controlfile $control -intree $intree -outprefix $intree -audify
# Output file of branch vectors is vidify_table2 (other files produced can be ignored).

# Branch vectors converted to waveforms.
# Process will take anywhere between one minute and several hours depending on tree size.
# The script will automatically change some variables according to number of branches input,
# such as track length (currently one of either 30, 60, 120 seconds) and amplitude decay.
R < audify_phylogeny.R --vanilla --slave
# Output is a set of Wav files, corresponding one for each branch.

# Manual steps required here for combining Wav files. 
# Open Audacity and do the following:
#	Import audio tracks (no more than about 1000 at a time), 
#	Can take some minutes and there is no progress bar, but scroll bar appears on RHS when done.
#	For each set imported: Select all, Tracks, mix, mix and render.
#	When all imported: 
#	Effect, base and treble, manually increase bass, click apply, close.
#	Select track, effect, normalize, normalize peak amplitude to -4 db, click OK.
# 	File, Export, export as Wav, filename Audio.wav

# Optional, make a video corresponding to the audio track. Input file again vidify_table2.
# Needs to have installed R libraries: plotrix, scales, tuneR
# Image dimensions have been precisely set, close to widescreen ratio 16/9,
# though slightly off to ensure divisible dimensions (a requirement for adding subtitles).
# It will remove if present, all previous video frames (video_frame.*)
# Following command requires Perl script (sort_commands_for_plotting_frame.pl) in working directory.
R < plot_frames.R --vanilla --slave
# Output is a load of JPEGs, numbering 25(fps) * track length

# Combine JPEG frames to video.
ffmpeg -i video_frame.%d.jpeg Video.mpeg

# Combine audio and video.
ffmpeg -i Video.mpeg -i Audio.wav -c:v copy -c:a aac AudioVideo.mp4

# Add subtitle, I recommend having this additional identifier on the video itself to avoid mixups.
# Standard format for subtitles used, e.g. this will add a single line to start of video:
# 1
# 00:00:02,000 --> 00:00:07,000
# Chesters et al 2023, Diptera
# Default font size doesnt look great, so is reduced.
ffmpeg -i AudioVideo.mp4 -vf subtitles=subtitles.srt:force_style='Fontsize=8' AudioVideo.srt.mp4

