#!/usr/bin/perl
my $arg_string  = join ' ', @ARGV;

# use Math::Trig;


#######################################################################################################3
# 
# 
# 
# 
# 
# 
# OPTIONS:
# 
# -label_tips [tip_label_string] [tip_label_file]
# 	takes members found in file, then finds them in the tree and adds the specified string to them
# 
# plot_circular_tree_with_branchlengths
# plot_cartoonized_clades 
# 
# 
# 
# 
# 
# 
# Jan 2017: decided to make a new script for processing newly inferred newick trees ,
# 	always after making tree i will process to add taxonomic information, perhaps re-root.
#	Previously this was a bit limited (for example cant read branchlengths, only key tax ranks used),
# 	Would be nice also to have option to label internal nodes.
# 	besides, i have a new Newick reader which makes re-rooting easy, so time to start anew.
# 
# 
# 	contains some code from bagpipe_phylo (getting shared taxonomic names)
# 
# 
# 	prints node labels with unique node ids (that can be used for rooting and other things),
# 		counts of terminals, taxonomic names
# 
# 	note, reroots at node not internode, so this will make a trifurcation at new base
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 	change log
# 	2017 JAN 22, v1.01: Got it working on a large complex-labelled tree.
# 	2017 JAN 25, v1.02: Streamline retrieval of terminals (ignore those without tax info),
# 		and their sampling for inference of shared taxa.
# 	2017 FEB 02, v2.01: included tree plotter
# 	2017 FEB 18, v2.02: Option to label user specified subset of tips with a string 
#		(output is .process_newick4)
# 	2017 MAR 07, v2.03: Hack working for applying higher taxa where only genus names given in labels
# 	2017 MAR 09, v2.04: Bugfix tree drawing, some tip branches were overlapping
# 		plot rectangular in addiition to circluar tree
# 	2017 MAR 20, v2.05: Ignore lineage inforamtino only given till ordinal rank (BOLDs), 
#		when inferring taxonomies
# 		option to fix newick strings, fully bifurcate, remove non-splitting internal nodes
# 	2017 AUG 12, v2.06: Option to reduce length of very long branches.
# 	2017 AUG 13, v2.07: Polytomies are now correctly plotted.
# 	2017 AUG 15, v2.08: Minor, terminating semicolon was missing from rerooted newick outfile
# 	2017 SEP 02, v2.09: Prints log with features of the tree
# 		new option available, -replace_scientific_notation_in_branchlengths
# 	2017 OCT 03, v2.10: Prints simple list of terminals in input tree (with taxonomic lineage info),
# 		sounds like silly little thing, but for huge trees can be useful,
# 	2018 FEB 11, v2.11: Multiply_all_branchlengths given in command, i got bored hashing out and not so
# 	2018 MAR 07, v2.12: New function to combine node labels from 2 different trees.
# 		works even if rerooting ...
# 		and should work in partially overlapping trees, it matches nodes via lists of terminals
#		new node labels have support value from tree given in -intree ,
#		followed by underscore, followed by support given in tree in -combine_node_labels
#	2018 MAR 19, v2.13: Plots support from second user tree, including bootstrap_IC
# 	2018 APR 07, v2.14: Couple of things transfered to command line, since, irritating if forget to change them
# 	2018 APR 09, v2.15: Prints tip labels, majority taxon for 3 ranks. also prints colorful tips.
# 			command option to draw circles at tips, if tip has string match: -highlight_tip_string_match
# 	2018 DEC 20, v2.16: Option to do some minor processing of tips
# 	2019 SEP 21, v2.17: Corrects minor root format error (made by bioperl reroot function i think)
# 	2019 SEP 22, v2.18: Function for coloring lineages from user specified terminals
# 	2019 OCT 23, v2.19: NA
# 	2019 DEC 04, v2.20: Prints trees with family or subfam as node labels
# 	2020 MAY 09, v2.21: Bugfix tip coloring with circles.
# 	2020 JUL 23, v2.22: Produces output file for depicting phylogeny as audio (sonification)
# 	2021 JAN 18, v2.23: On tree figure, plot pie charts for trait probabilities
# 	2021 MAY 29, v2.24: Table output for making video frame, corresponding to audification
# 	2021 JUN 06, v2.25: Some annotations made while transfering code to a seperate script (plotting function to relational_constraints).
# 	2021 AUG 11, v2.26: prints counts of dervied terminals for each node to vidify table, useful parameter in prettyfying
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#######################################################################################################3


# OPTIONS:



$randomly_resolve_polytomies = 0; # not yet working


 
# for labelling tips:
# @ranksarray 	= (" family");#, " family"
# why no family, subfamily tribe for flies?
# @ranksarray 	= ( " superfamily" , " family", " subfamily", " tribe");#, " family"
# @ranksarray 	= (" order" , " family");#, " family"
 @ranksarray 	= (" family"," subfamily", " tribe");#, " family"


my $remove_branchlengths_from_input_tree = 0; 

my $remove_node_support_from_input_tree = 1;

my $process_terminals = 0; # remove quotes 'Acropteris_sp.'

# for nodes with lots of decendents, reduce the number of them which are used to infer the taxonomic name:
my $reduced_comparisons = 10;
$reduce_by_proportion = 0.2;


# for tip labels, might not work if = 0
$uppercase_tax_strings = 1;


$verbose_log = 0;  # 1 == a few more details, but runs slower on big tree


# $plot_tree	= 0;

# for plotting trees:
$branch_plot_type 		= 2; # default = 2; 1= shortcut single line branching
$multipleier 			= 6; # Circumferal range


# $image_settings 		= "png\(\"outfile_name.png\"\, 2000\, 2000\)\n";
 $image_settings 		= "png\(\"outfile_name.png\"\, 9000\, 9000\)\n"; # upper limit is 10,000*10,000


	
						# aug 2017: 0.02
$reduce_very_long_branches_to_this = 1000; # to inactivate , put large number (1000).

# **** $multiply_all_branchlengths = 200; # hash out to inactivate *** now given in command

# settings for circular plotted tree:

$default_node_color = "black"; # grey black gray33
$clade_size_limit_for_interal_label = 10;

					# if zooming into huge tree, need low number here
$notch_size 			= 0.001; # size of segments for circumferal branches. 
						# 0.001 was too big in zoomed in 90000 sp tree (some noticable as not plotted) 
						# for huge ITOL2 tree, chhose 0.00002

$big_clade_color 		= "pink";
$truncate_new_taxonomic_strings		= 0;
# $plot_circular_tree_with_branchlengths = 1;# transfered to command arguments
# $print_internal_labels_circular_tree = 1; # now in command
$scale_branch_thickness = 0; # since greater density of branches to the tips, decrease thickness so can see them
$branch_root_to_tip_thickness_scale_factor = 1.01; # the degree to which thickness is decreased.
							# val=1 was indistinguishable in 350,000 speices tree
						# large number == very thin branches






######################################################################

# main parameters to change for plotting sub-trees for ITOL2

$print_tips_circular_tree = 0;
$all_tip_labels_eqidistant = 0;
$plot_tip_label_offset	= 0.10; # dont start tip label but against node, little space

$ref_tree_internal_label_cexB 	= 8;
$width 				= 25.0;	# BRANCH WIDTH; default = 1
	# if no branch lengths, probably both -1,1 
$xlim_circ	= "c(-2,2)"; $ylim_circ	= "c(-2,2)"; # 

# xlim = c(0.2,1), ylim = c(1,1.8)

$reference_tip_label_x = 1.9; # majority tax label given to tip, how far out to plot the label
	$append_terminal_label_with_fam_name = 1;
$tip_tax_label_cex = 9.0;

$plot_rainbow = 0;



	$rainbow_line_width = 12; # 2
 	$num_labs_to_print_rainbow = 230;
$plot_root_node_branchlength = 0.00025; # the very first branch, would be far left of a rectangualr 
					# < AUG2017 = 0.2; 0.75

######################################################################


# $internal_names_circular_tree

				# 202108: $pie_radius 		= 0.007;$pie_radius_border 	= 0.009;
$pie_radius 		= 0.004;
$pie_radius_border 	= 0.005;# border, different color for prior known states, and predicted states


		# TT project:
		#	 $trait_state_colors{"A"} = "red";
		#	$trait_state_colors{"B"} = "blue";
		#	$trait_state_colors{"C"} = "green";
		#	$trait_state_colors{"D"} = "yellow";
		#	$trait_state_colors{"E"} = "red";
		#	$trait_state_colors{"U"} = "yellow";
		#	$trait_state_colors{"S"} = "blue";

		# hive MBC project:
			$trait_state_colors{"H"} = "green";
			$trait_state_colors{"T"} = "brown";
			$trait_state_colors{"S"} = "blue";
			$trait_state_colors{"V"} = "red";
			$trait_state_colors{"G"} = "yellow";



# $highlight_tip_string_match = 0; # find tips with string below, then draw circles at these
# $tip_string_to_highlight = "CK_OTU";

$highlight_tip_string_cex = 5; # circular point size, 
$highlight_tips_from_file = 0;

# not sure what following is for, the bar x is set by: $rainbow_X1 = 1.5;$rainbow_X2
$highlight_tips_X = 1.9;# 0.992

$highlight_tips_eqidistant = 0;#  0 = highlight tip at tip, 1= hoghligt tip outwards a bit and all at same radius
# $highlight_file = "/home/douglas/scripted_analyses/BEF_China_and_TreeDi/MQ_paper_2/data/OTU_colors"; # OTU_abund_mt30_2
# $highlight_file = "/home/douglas/scripted_analyses/nCov19/data/haplotype_colors.txt";
# $highlight_file = "/home/douglas/scripted_analyses/tingting/analysis4/current_trait";

      


# $highlight_file = "specieslist_for_treeplot.U_BodyLength"; $rainbow_X1 = 1.60;$rainbow_X2 = 1.63;
# $highlight_file = "specieslist_for_treeplot.U_ITD"; $rainbow_X1 = 1.64;$rainbow_X2 = 1.67;
# $highlight_file = "specieslist_for_treeplot.U_HairLength"; $rainbow_X1 = 1.68;$rainbow_X2 = 1.71;
# $highlight_file = "specieslist_for_treeplot.Sociality"; $rainbow_X1 = 1.72;$rainbow_X2 = 1.75;
# $highlight_file = "specieslist_for_treeplot.Nest_location"; $rainbow_X1 = 1.76;$rainbow_X2 = 1.79;
# $highlight_file = "specieslist_for_treeplot.Parasitism"; $rainbow_X1 = 1.80;$rainbow_X2 = 1.83;
 $highlight_file = "specieslist_for_treeplot.Lecty"; $rainbow_X1 = 1.84;$rainbow_X2 = 1.87;

# $highlight_file = "subject_species_list";


# file looks like:
# OTU63	red
# OTU64	red
# OTU65	gray

$highlight_tips_type = 2; # 1=circle, 2=bar
$colored_tip_bar_length = 0.005;



# $tip_bar_color = "brown"; 	# color now defined in the file, can be differnt for different otu

				# green, XSBN, 1.02-1.04; 
				# red,     BJ, 1.04-1.06; 
				# blue,   QTP, 1.06-1.08; 
				# brown,   XJ, 1.08-1.10; 

# $highlight_file = "CN_obs.XJ"
#	$rainbow_X1 = 1.08;$rainbow_X2 = 1.10;
# $tip_bar_color = "brown";






###############################################################################################


	# RECTANGUALR:


$xlim 				= "c(0.1,6)";# regular:c(-16, 16)
# $image_settings_RCT 		= "png\(\"outfile_name_RECTANGULAR.png\"\, 2000\, 2500\)\n";
# $image_settings_RCT 		= "postscript\(\"outfile_name_RECTANGULAR.eps\"\, width=3000\, height=4800\)\n";
# $image_settings_RCT 		= "jpeg\(\"Global_sp_accum.jpg\"    , width = 1200 , height = 2400, res=2000\)\n";
$image_settings_RCT 		= "tiff\(filename = \"ITOL2_species_level.tiff\", width = 8000 , height = 9000, units=\"px\", res=600, compression = \"jpeg\"\)\n";




$rectangular_plot_branch_width = 2.0;
$rectangular_plot_tip_label_cex = 0.2;	# WP=0.6
$rectangular_plot_support_circle_size = 1.2;
$plot_tip_label_offset_rectangular	= 0.05; # dont start tip label but against node, little space
$rectangular_tree_internal_label_cex = 10;# WP=0.6
$print_internal_node_labels_on_rectangular_tree = 1;
$print_tips_rectangular_tree = 0;
$cartooned_clade_external_label_cex = 1.4;
$rect_branch_color_default = "gray33";


$ylim 	= "c(2,10)";

$branchnumber = 1; # for vidify

#############################################################################################


# $collapse_nodes{"INTERNAL_NODE_997"} = "Beetles";
# $collapse_nodes{"INTERNAL_NODE_1107"} = "Flies";
# $collapse_nodes{"INTERNAL_NODE_1362"} = "Moths";
# $collapse_nodes{"INTERNAL_NODE_1467"} = "Wasps";


$plot_ITOL2_9_trees = 0;

if($plot_ITOL2_9_trees == 1)
	{plot_ITOL2_9_trees()};




$plot_climate_tree = 0;

if($plot_climate_tree == 1)
	{

$draw_clades_climate = "
text(-0.0114318983366303, 0.00546884464086701, col = \"black\", labels=\"Collembola\",cex=5, font=2)
text(-0.00883818726704874, 0.00681024632191072, col = \"black\", labels=\"Paraneoptera\",cex=5, font=2)
text(-0.00662308169591456, 0.011915531503297, col = \"black\", labels=\"Coleoptera\",cex=5, font=2)
text(0.00433795792349205, 0.014285396862067, col = \"black\", labels=\"Hymenoptera\",cex=5, font=2)
text(0.0161155129122418, -0.000173161468182245, col = \"black\", labels=\"Diptera\",cex=5, font=2)
text(0.0111938503802918, -0.0129320928596986, col = \"black\", labels=\"Trichoptera\",cex=5, font=2)
text(-0.00904688199796051, -0.0145153449212619, col = \"black\", labels=\"Lepidoptera\",cex=5, font=2)
";

	};



###############################################################################################




#####################################
read_command_arguments($arg_string);#
#####################################




%paint_clades;

$paint_clades_circular_tree = 0;	# these are a series of user defined nodes, subtree of which to color
$paint_clades_from_file = 0;

$big_clade_size 		= 2;# deactivated

# my $paint_nodes_file = "/home/douglas/scripted_analyses/insect_temp_shifts/bayestraits/koppen_A_nodes";
# my $paint_nodes_file = "/home/douglas/scripted_analyses/insect_temp_shifts/bayestraits/koppen_nodes";
my $paint_nodes_file = "/home/douglas/scripted_analyses/insect_temp_shifts/bayestraits/koppen_nodes_unfiltered_gt0.95";



# alternative way to color branches, define terminal subset, then coleor each all the way up to root

$color_lineages_from_file = 0;
my $color_lineages_file = "color_lineage_tax_list.blue";  # color_lineage_tax_list.blue 
$color_lineages_this = "red";





$label_painted_nodes = 0;
$painted_branch_label_cex = 1.5;

$draw_tree_R_commands6 = "";
$draw_tree_R_commands7 = "";


if($paint_clades_circular_tree == 1){paint_clades_circular_tree()};



if($highlight_tips_from_file == 1)
	{
	open(HF, $highlight_file) || die "\nerror 342. you specified highlight tips from file, but cant open file:$highlight_file\n";
	while (my $line = <HF>)
		{
		$line =~ s/\n//;$line =~ s/\r//;
		# hack so can read this file, dont need another parser
# CONSTR number:1473, termnial:Xizicus, taxon_assigned:Xizicus, count_new:3 new_members_incl:Xizicus_fascipes	Xizicus_howardi	Xizicus_sp_SM_2011_HQ609311
		if($line =~ s/CONSTR number:\d+\, termnial:\w+\, taxon_assigned:(\w+)\, count_new:\d+ new_members_incl:(\w+).+/$2	$1/)
			{print ""};
		# normally require two columns
		if($line =~ /^(\S+)\t([\w\d]+)/)
			{
			my $terminal = $1; my $colour_to_plot = $2;$highlight_tips{$terminal} = $colour_to_plot;
			
			}elsif($line =~ /^(\S+)/)
			{
			my $terminal = $1; $highlight_tips{$terminal} = "red";

			$errorprintting++;if($errorprintting <= 3){print "tip to highlight given but no color, using red\n"};

		#	die "\nneed 2 colomns in highlight tip file, give termainl and color to plot. quitting.\n";
			};

		# check later is that hash $highlight_tips{... has something in it

		};
	close HF;

my @highlight_keys = keys %highlight_tips;print "will try find $#highlight_keys tips to higlight\n";
# 12_BD00927758 Ibidoecus_bisignatu Pielomastax_zheng Ancistronycha_abdominalis Brachythemis_contaminat Hydrochus_BBCCM42410_JN290349 Scaphidium_quadriguttatum Chorotypus_fenestratu Agabus_affinis Campanulotes_bidentatu Hesperobaenus_c
	};





# $plot_cartoonized_clades = 1; # transfered to command arguments

if($plot_cartoonized_clades == 1)
{
$cartoonize_clades{'INTERNAL_NODE_347'} = 'Collembola';

# $cartoonize_clades{'INTERNAL_NODE_375'} = 'Phasmatidae';

$cartoonize_clades{'INTERNAL_NODE_27'} = 'Vespidae';
# $cartoonize_clades{'INTERNAL_NODE_51'} = 'Sphecinae';
$cartoonize_clades{'INTERNAL_NODE_43'} = 'Crabroninae';
#$cartoonize_clades{'INTERNAL_NODE_69'} = 'Apidae';
#$cartoonize_clades{'INTERNAL_NODE_85'} = 'Halictidae';
$cartoonize_clades{'INTERNAL_NODE_170'} = 'Braconidae';
$cartoonize_clades{'INTERNAL_NODE_177'} = 'Ichneumonidae';
$cartoonize_clades{'INTERNAL_NODE_137'} = 'Chrysidoidea';
#$cartoonize_clades{'INTERNAL_NODE_252'} = 'Noctuoidea';
#$cartoonize_clades{'INTERNAL_NODE_244'} = 'Bombycoidea';
$cartoonize_clades{'INTERNAL_NODE_313'} = 'Schizophora';
$cartoonize_clades{'INTERNAL_NODE_325'} = 'Culicomorpha';

$cartoonize_clades{'INTERNAL_NODE_370'} = 'Dictyoptera';
$cartoonize_clades{'INTERNAL_NODE_390'} = 'Orthoptera';
$cartoonize_clades{'INTERNAL_NODE_14'} = 'Tenthredinoidea';
 $cartoonize_clades{'INTERNAL_NODE_206'} = 'Cucujiformia';
# $cartoonize_clades{'INTERNAL_NODE_3'} = 'Aphidomorpha';
$cartoonize_clades{'INTERNAL_NODE_424'} = 'Panheteroptera';

$cartoonize_clades{'INTERNAL_NODE_231'} = 'Neuropterida';
$cartoonize_clades{'INTERNAL_NODE_272'} = 'Ditrysia';
$cartoonize_clades{'INTERNAL_NODE_292'} = 'Trichoptera';
 $cartoonize_clades{'INTERNAL_NODE_379'} = 'Phasmatodea';
$cartoonize_clades{'INTERNAL_NODE_436'} = 'Sternorrhyncha';
$cartoonize_clades{'INTERNAL_NODE_433'} = 'Euhemiptera';
$cartoonize_clades{'INTERNAL_NODE_148'} = 'Chalcidoidea';
$cartoonize_clades{'INTERNAL_NODE_209'} = 'Scarabaeidae';
$cartoonize_clades{'INTERNAL_NODE_93'} = 'Apoidea';

$cartoonize_clades{'INTERNAL_NODE_53'} = 'Sphecoidea';

};







###############################################################################################









#######################################################################################################3


#######################################################################################################3










#################################################################################################

if($remove_polytomies	== 1){remove_polytomies($intree)};

#################################################################################################

if($replace_scientific_notation_in_branchlengths == 1)
	{
	replace_scientific_notation_in_branchlengths(); exit;
	};

#################################################################################################








# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;
%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object
%assign_these_new_terminals_to_taxon;

###############
store_nodes();#
parse_namesfile();#
###################

$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;
print "name for the starting node ($starting_node) which has been specified:$starting_name\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


#######################################################################################################3



$interal_node	= 0;
$root_node;
$bls_read;
$count_terminals;
$minimum_branchlength = 9999;
$maximum_branchlength = 0;
$count_scientific_notation_branchlengths_removed=0;
$count_regular_format_branchlengths_removed=0;
$count_zero_length_branchlengths_removed=0;
$bifurcations = 0;
$polytomies = 0;


#########################
record_tree3($intree);#
#########################


my @polytomys = keys %polytomy_IDs;
my $polytomy_string = join( ',',  @polytomys); $polytomy_string =~ s/^(.{0,100}).+/$1 . "......"/e;

print "polytomies_randomly_resolved:$polytomies_randomly_resolved\n";

			# test:"INTERNAL_NODE_2346"; # $root_node;
			# Pol rooted at INTERNAL_NODE_2609




# $combine_node_label_tree = $1; 
if($combine_node_labels == 1)
	{
	combine_node_labels();
	};



if( $color_lineages_from_file == 1){color_lineages_from_file()};



if($arg_string =~ /\-reroot_at\s+(\S+)/)
	{
	$traversal_start_node = $1;
#	$count_terminals_node_lable{$traversal_start_node} = $count_terminals;
	}else{
	$traversal_start_node = $root_node; 
	};



$internal_nodes_with_tax_assigned = 0;
$count_sub_calls;

open(LOG, ">$outprefix.process_newick_LOG") || die "";
print LOG "current_node_ID\tcount_termnials\tshared_tobycode\tshared_lineage\n";
print "\ncalling assign_names_to_internal_nodes with $traversal_start_node , New_Root\n";

# In addition to names, this sub makes terminal counts and write in object count_terminals_node_lable. These are used extensivly in plotting.

#######################################################################
assign_names_to_internal_nodes($traversal_start_node , "New_Root");#
#######################################################################

close LOG;




if($identify_nodes_via_decendents == 1)
	{
	if($decendent_defined_nodes_found_in_tree >= 1)
		{
		print "of the decendent defined nodes specified in file, $decendent_defined_nodes_found_in_tree were found in the tree\n";
		}else{
		die "\ndid not find any of the decendent defined nodes.\n";
		};
	};






if($arg_string =~ /\-reroot_at\s+(\S+)/)
	{
	$count_terminals_node_lable{$1} = $count_terminals;
	};



print "
internal_nodes_with_tax_assigned:$internal_nodes_with_tax_assigned
";

print "traverse tree from start node $traversal_start_node\n";

# simple newick:
$new_newick_string = "($traversal_start_node)";
# with internal node labesls:
$new_newick_string2 = "($traversal_start_node)";
# will get higher tax names put on tips:
$new_newick_string3 = "($traversal_start_node)";
$new_newick_string4 = "($traversal_start_node)";
# for combined node labels:
$new_newick_string5 = "($traversal_start_node)";

# label internal nodes only if family
$new_newick_string_family = "($traversal_start_node)";

# label internal nodes only if SUBfamily
$new_newick_string_subfamily = "($traversal_start_node)";


if($sonify == 1)
	{
	open(SONIFYFILE, ">$sonify_file" ) || die "\nerror 621.\n";

	};


if($plot_tree == 1)
	{
	print "\nuser opted to PLOT TREE. intiating newick strings and plotting commands.... ";
	# plotting:
	$draw_tree_R_commands = "";$draw_tree_R_commands3 = "";$draw_tree_R_commands4 = "";$draw_tree_R_commands4b = "";
$draw_tree_R_commands_tip_labels = "";
	$draw_tree_R_commands_RECTANGL = "";$draw_tree_R_commands_RECTANGL2 = "";
	$draw_tree_R_commands_RECTANGL_internal_node_labels = "";


	open(OUT8 , ">r_commands_filename") || die "\nerror 31. cant open output file for writing r script ($r_commands_filename)\n";
	print OUT8 "
library(\"plotrix\")
library(mapplots)

colorfunc = colorRamp(c(\"blue\",\"white\",\"red\"))
";
	print OUT8 $image_settings;
	print OUT8 "plot(0, type = \"n\", xlim = $xlim_circ, ylim = $ylim_circ,";
	print OUT8 " xlab = \"\",ylab = \"\",xaxt = \"n\", yaxt = \"n\", bty = \"n\")\n";

	open(OUT_RECTANGL , ">r_commands_filename_RECTANGL") || die "\nerror 241. cant open output file for writing r script ($r_commands_filename)\n";
	print OUT_RECTANGL $image_settings_RCT;
	print OUT_RECTANGL "colorfunc = colorRamp(c(\"red\" , \"green\"))
		";
 # colorfunc = colorRamp(c("red" , "gray", "gray", "gray", "gray" , "green"))
 # colorfunc = colorRamp(c("red" , "white", "white", "white", "white", "white", "white", "white", "white" , "green"))
 # colorfunc = colorRamp(c("red" , "gray", "gray", "gray", "gray", "gray", "gray" , "green"))


	print OUT_RECTANGL "plot(0, type = \"n\", xlim = $xlim, ylim = $ylim,";
	print OUT_RECTANGL " xlab = \"\",ylab = \"\",xaxt = \"n\", yaxt = \"n\", bty = \"n\")\n";

	open (RECT_PLOT_LOG , ">RECT_PLOT_LOG"); # prints positions of tips so you can decide where boxes go
	print ".. done\n"; 
	};


print "calling sub traverse_tree, start node:$traversal_start_node\n";# builds up new newick strings

open(NODE_LOG, ">process_newick.node_log") || die "\nerror 373.\n";

open(VIDIFY, ">vidify_table") || die "";
# open(VIDIFY2, ">vidify_table2") || die "";

# if rerooting, count of number of terminals from new root, needs redefining. 
# $count_terminals_node_lable


######################################################
	#     $current_node;          $from_parent;     $sum_branchlength;                  $node_plot_y 
traverse_tree($traversal_start_node , "New_Root" ,      $plot_root_node_branchlength ,      1 ,           $default_node_color  ,  0 );
######################################################

print "
\tnodes_traversed:$nodes_traversed\n";

close NODE_LOG;
close VIDIFY;
# close VIDIFY2;



unless($new_newick_string4 =~ s/^\((.+\))[^\)]+\)$/$1;/){print "\nunexpected. 137.\n"};
open(OUT1, ">$outprefix.process_newick1") || die "\nerror 48\n";
print OUT1 "$new_newick_string4\n";
close OUT1;


unless($new_newick_string =~ s/^\((.+\))[^\)]+\)$/$1;/){print "\nunexpected, 260 .... continuing\n"};
open(OUT0, ">$outprefix.process_newick0") || die "\nerror 48\n";
print OUT0 "$new_newick_string;\n";
close OUT0;

unless($new_newick_string_family =~ s/^\((.+\))[^\)]+\)$/$1;/){print "\nunexpected, 260 .... continuing\n"};
open(OUT01, ">$outprefix.process_newick_fam") || die "\nerror 48\n";
print OUT01 "$new_newick_string_family;\n";
close OUT01;

unless($new_newick_string_subfamily =~ s/^\((.+\))[^\)]+\)$/$1;/){print "\nunexpected, 260 .... continuing\n"};
open(OUT02, ">$outprefix.process_newick_subfam") || die "\nerror 48\n";
print OUT02 "$new_newick_string_subfamily;\n";
close OUT02;

if($combine_node_labels == 1)
	{
unless($new_newick_string5 =~ s/^\((.+\))[^\)]+\)$/$1;/){print "\nunexpected. 137.\n"};
open(OUT5, ">$outprefix.process_newick5") || die "\nerror 48\n";
print OUT5 "$new_newick_string5\n";
close OUT5;	
	};


if($sonify == 1)
	{

	# probably best thing to record would be:
	#	for each node, x position of bifurcation, y positions of each child branch.
	# 	also child branch length, maybe just the shortest one.
	#	store according to node id.
	# 	sort according to y, omit nodes at exactly the same y.
	
#	$sonify_node_x{$current_node} = $sum_branchlength;
#	unless($sonify_child_nodes{$current_node} =~ /\t$test\t/){$sonify_child_nodes{$current_node} .= "	$test	"};

	my @sonify_nodes = keys %sonify_node_x; # KEY is nodeID, ENTRY = y position
	@sonify_nodes = sort @sonify_nodes;
	
	foreach my $sonify_node(@sonify_nodes) # nodeIDs
		{
		my $x_start = $sonify_node_x{$sonify_node};my $y_start = $sonify_node_y{$sonify_node};
	#	print "sonify_node:$sonify_node x_start:$x_start\n";

		if($sonify_child_nodes{$sonify_node} =~ /([^\t]+)\t+([^\t]+)/)
			{
			my $child1 = $1; my $child2 = $2; # print "\tchild1:$child1 child2:$child2\n";
	
			if($sonify_child_xy{$sonify_node}{$child1} =~ /(.+)\t(.+)/)
				{#	$sonify_child_xy{$current_node}{$test} = "$assign_new_x\t$y2_proportion";
				my $child1_x_end = $1; my $child1_y = $2;

				if($sonify_child_xy{$sonify_node}{$child2} =~ /(.+)\t(.+)/)
					{
					my $child2_x_end = $1; my $child2_y = $2;
					my $shortest_x = $child1_x_end;if($child1_x_end > $child2_x_end){$shortest_x = $child2_x_end};

					#	starts with NODE_ID       
					# then 4 values for each node, x start and end which is same for both child nodes, and y of each child.
					print SONIFYFILE "$sonify_node\t$x_start\t$shortest_x\t$child1_y\t$child2_y\t$y_start\n";
					# print "sonify_node:$sonify_node\tx_start:$x_start\tshortest_x:$shortest_x\t$child1_y\t$child2_y\t$y_start\n";
					};

				};



			}else{
			print "warning cant sonify node\n";
			};

		};


	


#	my @sonify_coomands = keys %sonify_commands;@sonify_coomands = sort @sonify_coomands;
#	foreach my $sonify_command(@sonify_coomands)
#		{
#		print "sonify_command:$sonify_command\n";
#		if($sonify_command =~ /(.+)__(.+)/)
#			{
#			my $x = $1; my $y = $2;print SONIFYFILE "\"$x\",\"$y\"\n";
#			};
#		};

	print "sonify positions written to file $sonify_file\n";

	close SONIFYFILE;

	};




if($plot_tree == 1)
	{

	if( $print_terminal_tax_labels == 1)
		{
		print_terminal_tax_labels();
		if($plot_rainbow == 1)
			{print_rainbow()}


		};



	print OUT8 $draw_tree_R_commands;
	print OUT8 $R_commandtip;
	print OUT8 $draw_tree_R_commands3;
	print OUT8 $draw_tree_R_commands4;
	print OUT8 $draw_tree_R_commands4b;
	print OUT8 $draw_tree_R_commands6;
	print OUT8 $draw_tree_R_commands7;
	print OUT8 $draw_tree_R_commands_tip_labels;
	print OUT8 $draw_clades_climate;

	print OUT8 "dev\.off\(\)\n";close OUT8;

	print OUT_RECTANGL $draw_tree_R_commands_RECTANGL3; # discontinued
	print OUT_RECTANGL $draw_tree_R_commands_RECTANGL;# segments, points, text
	print OUT_RECTANGL $draw_tree_R_commands_RECTANGL2; # node labels
	print OUT_RECTANGL $draw_tree_R_commands_RECTANGL_internal_node_labels; # large node labels

	print OUT_RECTANGL "dev\.off\(\)\n";close OUT_RECTANGL;

	close OUT_RECTANGL; close OUT8; close RECT_PLOT_LOG;

	print "\n\nplot_tree == 1, written r commands:
	R < r_commands_filename --vanilla --slave\n\tOR:
	R < r_commands_filename_RECTANGL --vanilla --slave\n\tquitting.\n";
	exit;
	};





# put strings on specific user defined tips:
if($tip_label_string =~ /./)
	{

		###############################################
	$tr = label_user_defined_tips($new_newick_string4);#
		###############################################

	open(OUT4, ">$outprefix.process_newick4") || die "\nerror 48\n";
	print OUT4 "$tr\n";
	close OUT4;

	print "\nlabelled subset of tree tips as specified by user. FIN.\n\n";
	exit;
	};





open(LOG2, ">outprefix.process_newick_LOG2") || die "";
print LOG2 "
intree:$intree
	outprefix:$outprefix
	starting_node:$starting_node
	tip_label_file:$tip_label_file

remove_branchlengths_from_input_tree:$remove_branchlengths_from_input_tree
	remove_node_support_from_input_tree:$remove_node_support_from_input_tree
	reduce_very_long_branches_to_this:$reduce_very_long_branches_to_this
	number_of_branches_over_user_specified_limit:$number_of_branches_over_user_specified_limit
	multiply_all_branchlengths:$multiply_all_branchlengths
	remove_polytomies:$remove_polytomies

traversal_start_node:$traversal_start_node
	internal_nodes_with_tax_assigned:$internal_nodes_with_tax_assigned
	nodes_traversed:$nodes_traversed

details on tree read:
	minimum_branchlength:$minimum_branchlength
	maximum_branchlength:$maximum_branchlength
	count_scientific_notation_branchlengths:$count_scientific_notation_branchlengths_removed
	count_regular_format_branchlengths:$count_regular_format_branchlengths_removed
	count_zero_length_branchlengths:$count_zero_length_branchlengths_removed

newick string stored, 
	bifurcations:$bifurcations
	monotomies:$monotomies
	count of polytomies:$polytomies
		polytomy node IDs:$polytomy_string
	count of internal nodes:$interal_node
	count terminals:$count_terminals
	branchlengths read:$bls_read
	newick_string:$newick_string
	duplicated_terminal_labels:$duplicated_terminal_labels
	nodes_collapsed:$nodes_collapsed
\n";
close LOG2;



if($write_subtrees == 1)
	{
	print "\nwrite_subtrees subtrees_written:$subtrees_written\n";

#	fisher_yates_shuffle( \@list_of_bayestraits_commands );
#	open(BAYESTRAITS_BASH_COMMANDS, ">bayestraits_bash_commands") || die "";
#	print BAYESTRAITS_BASH_COMMANDS @list_of_bayestraits_commands;
#	close BAYESTRAITS_BASH_COMMANDS;


	die "\nyou specified write subtree, presuming you dont need tip modified tree writing. exiting.\n";
	};



# family names on termials:
print "putting taxonomic names ontp terminals\n";
$lineages_assigned =0;


################################################
my $tr = proccess_tree($new_newick_string4);#
################################################

################################################
my $tr_df = proccess_tree_differently($new_newick_string4);#
################################################



open(LOG2, ">>$outprefix.process_newick_LOG2") || die "";
print LOG2 "binomials_found_in_newick:$binomials_found_in_newick
	of these, lineages_assigned:$lineages_assigned
	newick_binomials_for_which_taxonomies_not_found:$newick_binomials_for_which_taxonomies_not_found
";
close LOG2;


# unless($tr =~ s/^\((.+\))[^\)]+\)$/$1;/){print "\nunexpected. 137.\n"};


open(OUT3, ">$outprefix.process_newick3") || die "\nerror 48\n";
print OUT3 "$tr\n";
close OUT3;

open(OUT6, ">$outprefix.process_newick6") || die "\nerror 48\n";
print OUT6 "$tr_df\n";
close OUT6;


print "

details on tree read:
	minimum_branchlength:$minimum_branchlength
	maximum_branchlength:$maximum_branchlength

 process_newick0 has no extra labeling, but is re-rooted if specified by user.
 process_newick1 has internal node labels in format: nodeID.count_terminal.internal_taxname
 process_newick3 has the same, and higher tax names added to terminal labesl
 process_newick6 same as 3, but might be more reliable

";




print "\n\nFIN\n\n";
exit;








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub record_tree3
{
my $treefile = shift;
print "sub record_tree3\n";
my $tree1= "";
open(IN, $treefile) || die "\n\nerror 561 cant open file named:$treefile\n\n";
while (my $line = <IN>)
	{if($line =~ /(.+\(.+)/){$tree1 = $1}};
close(IN);
$tree1 =~ s/ //g;$tree_length = length($tree1);

print "\nnewick string of length $tree_length has been read from file:$treefile\n\n";

# this format prodcued by soe rerooting function, correct for it:
# 438):0.199726792634503)Acerentomon_microrhinus);

if($tree1 =~ s/(\d\)\:\d\.\d+\))([A-Z][a-z]+_[a-z]+\)\;)$/$1:0.000001,$2/)
	{
	print "warning, corrected minor rooting format error in your tree\n";
	};



if($plot_trees_with_branchlengths == 1)
	{
	if($tree1 =~ /\:/)
		{
		print "user opted to plot_trees_with_branchlengths, and branch lengths have been found\n";
		}else{
		die "\nerror, you opted to print circular tree with branch lengths,\n\thowever there are no brnach lengths in your tree. quitting.\n"
		};
	};




if($remove_branchlengths_from_input_tree == 1)
	{
	print "\n820. user chosen to remove branch-lengths in input tree ...\n";
	# print "test 821\n";
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	print " scientific notation .. \n";
	while($tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/)
		{$count_scientific_notation_branchlengths_removed++}; 
	# remove regular branchlengths: 0.02048
	print "$count_scientific_notation_branchlengths_removed. regular branchlengths .. \n";
	while($tree1 =~ s/\:\-*\d+\.\d+//)
		{$count_regular_format_branchlengths_removed++}; 
	# remove 0 length branchlengths
	print "$count_regular_format_branchlengths_removed. 0 length .. \n";
	while($tree1 =~ s/\:\d+//)
		{$count_zero_length_branchlengths_removed++};

	print "\tdone\n$count_scientific_notation_branchlengths_removed count_scientific_notation_branchlengths_removed\n";
	print "\t$count_regular_format_branchlengths_removed count_regular_format_branchlengths_removed\n";
	print "\t$count_zero_length_branchlengths_removed count_zero_length_branchlengths_removed\n\n";
	if($tree1 =~ /..\:../){die "\nhuh\n"};

	}else{

	print "user chosen to retain branch-lengths in input tree ...\n";

	if($verbose_log == 1)
		{
		print "user opted for verbose log, getting branch length details ..... ";
		my $tree2 = $tree1;
		# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
		while($tree2 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/)
			{$count_scientific_notation_branchlengths_removed++}; 
		# remove regular branchlengths: 0.02048
		while($tree2 =~ s/\:\-*\d+\.\d+//)
			{$count_regular_format_branchlengths_removed++}; 
		# remove 0 length branchlengths
		while($tree2 =~ s/\:\d+//)
			{$count_zero_length_branchlengths_removed++};
		print " done\n";
		print "\tcount_scientific_notation_branchlengths:$count_scientific_notation_branchlengths_removed\n" ,
		"\tcount_regular_format_branchlengths:$count_regular_format_branchlengths_removed\n" ,
		"\tcount_zero_length_branchlengths:$count_zero_length_branchlengths_removed\n";
		};

	};


if($remove_node_support_from_input_tree == 1)
	{
	print "\n868. user chosen to remove branch-support in input tree ...\n";
	while($tree1 =~ s/(\))\d\.\d+/$1/)
		{$count_proportional_node_support_removed++}; 
 # ictus:0.05948)0.903:0.02798,(Cryptocercus_GMGSG37012_BD00817042:0.03212,Cryp

	print "
	count_proportional_node_support_removed:$count_proportional_node_support_removed
	";

	}else{
	if($tree1 =~ /(\))\d\.\d+/)
		{
		print "\nWARNING. Node support found, and script set to retain.  Uncertain behaviour.\n";
		}; 
	
	};


my $newick_string 	= $tree1;



#############################################################################################################

# third generation Newick reader, mainly to overcome limitation of reading direction. 

print "\nreading Newick tree ....\n";

while ($newick_string =~ s/\(([^\(\)]+)\)([0-9\.]*)/INTERNAL_NODE_$interal_node/)
	{# no point reading the adjacent branchlength, 
		# since the identity of the node to which it connects, is not known at this stage,
		# so read lengths to child nodes.

	my $node = $1; my $support = $2;my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";
	my @child_nodes = split /\,/ , $node;
	$child_counts{$nodeID} = $#child_nodes;#	print "\nnodeID:$nodeID \@child_nodes:@child_nodes\n";
	$node_support{$nodeID} = $support;#	print "\nnodeID:$nodeID \@child_nodes:@child_nodes\n";
	if($#child_nodes >= 2)
		{
		$polytomies++;$polytomy_IDs{$nodeID} = 1
		}elsif($#child_nodes == 0){$monotomies++
		}elsif($#child_nodes == 1){$bifurcations++}else{die "\nHuh? child_nodes:$#child_nodes\n"};

	unless($#child_nodes == 1)
		{
		print "non bifurcating, node $interal_node, count child nodes $#child_nodes\n";
		if(length($node)< 100){print "\toffending node:$node\n"};
		};

	if($interal_node =~ /0000$/){print "\tnodeID:$nodeID count_child_nodes:" , scalar @child_nodes , " bls_read:$bls_read minimum_branchlength:$minimum_branchlength maximum_branchlength:$maximum_branchlength\n"};

	for $i(0 .. $#child_nodes)
		{
		my $current_child = $child_nodes[$i];

		if($process_terminals == 1)
			{
			$current_child =~ s/^\'([A-Za-z\.\_]{2,30})\.\'/$1/; # 'Acropteris_sp.'
			$current_child =~ s/^\'([A-Za-z\.\_]{2,30})\'/$1/; # 'Jana_nr._eurymas'
			};
			

		my $current_child_branchlength = "NA";# length from node to child
		if($current_child =~ s/\:([0-9\.\-E]+)$//) # allows for scientific notation lengths sometimes given:7.59E-4
			{
			$current_child_branchlength = $1;$bls_read++;
			if($current_child_branchlength >= $maximum_branchlength){$maximum_branchlength = $current_child_branchlength};
			if($current_child_branchlength <= $minimum_branchlength){$minimum_branchlength = $current_child_branchlength};

			if($multiply_all_branchlengths >= 1.1)
				{
				$current_child_branchlength *= $multiply_all_branchlengths;$bls_multiplied++;
				};


			if($current_child_branchlength >= $reduce_very_long_branches_to_this)
				{

				if($number_of_branches_over_user_specified_limit <= 15)
					{
					print "warning, you specified to reduce branches over length $reduce_very_long_branches_to_this, current branch is $current_child_branchlength, modifying your tree\n";
					};
				if($number_of_branches_over_user_specified_limit == 16)
					{print " ....... storping printing warning ...\n"};

				$current_child_branchlength = $reduce_very_long_branches_to_this;
				$number_of_branches_over_user_specified_limit++;
				};


			}elsif($current_child =~ /\:/)
			{
			die "\nerror 1184. cant parse branchlength from:$current_child, quitting\n\n"
			};
		
#die "nodeID:$nodeID node:$node branchlength:$branchlength current_child:$current_child\n";

		my $current_child_number_connections;		
		if($child_counts{$current_child} =~ /\d/)
			{$current_child_number_connections = $child_counts{$current_child} + 1
			}else{$current_child_number_connections = 0
			};

		# record the child nodes for current node:
		$nodes{$nodeID}{$i} 			= $current_child;
		$child_counts{$current_child} = $current_child_number_connections;
		$parent_node_IDs{$current_child} = $nodeID;

		# record length of branch between the current two nodes,
		# store in both directions
		$branchlengths{$nodeID}{$current_child} = $current_child_branchlength;
		$branchlengths{$current_child}{$nodeID} = $current_child_branchlength;

		$nodes{$current_child}{$current_child_number_connections} = $nodeID;
	#	print "\t\tcurrent_child_number_connections:$current_child_number_connections\n";
	#	print "\ti:$i $current_child\n";

		unless($current_child =~ /^INTERNAL_NODE_/)
			{
			$terminal_labels{$current_child}++;
			if($terminal_labels{$current_child}>= 2)
				{
				$duplicated_terminal_labels++;
				print "\nwarning, duplicated terminal label:$current_child\n";
				};
			};

		};

	if(exists($collapse_nodes{$nodeID}))
		{
		$newick_string =~ s/$nodeID(\D)/$collapse_nodes{$nodeID}$1/;
		$nodes_collapsed++;
		};


	if(length($newick_string) <= 200)
		{
	#	print "newick_string:$newick_string\n";
		};	

	if($randomly_resolve_polytomies == 1)
		{
		# first get string within parentheses, next that will be encountered in main loop
		if($newick_string =~ /\(([^\(\)]+)\)/)
			{
			my $segment_next_loop  = $1;
			# for only this segment, check if it is trifurcation:
			if($segment_next_loop =~ /^[^\(\)\,]+\,[^\(\)\,]+\,[^\(\)\,]+$/)	
				{
				print "polytomy randomly resolved, INTERNAL_NODE_$interal_node\n";
				$polytomies_randomly_resolved++;
				# if so, swap it in whole newick string:
				$newick_string =~ s/\(([^\(\)\,]+)\,([^\(\)\,]+)\,([^\(\)\,]+)\)/"((" , $1 , "," , $2 , ")," , $3 , ")"/e;
				};

			};
		};

	$root_node = $nodeID;
	$interal_node++;

	};#while ($newick_string =~

#############################################################################################################


my @terminal_array = keys %terminal_labels; @terminal_array = sort @terminal_array;
$count_terminals = scalar @terminal_array;

open(TIPS, ">list_of_terminals_in_your_tree");

foreach my $tip(@terminal_array)
	{
	my $copy_id = $tip;$copy_id =~ s/_temp\d+\.\d+$//;
	my $completelineage = "NA";
	if(exists($complete_lineage_for_this_species{$copy_id}))
		{
		$completelineage = $complete_lineage_for_this_species{$copy_id};
		}else{
		$genus = $copy_id;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
	#	print "no sucess looking for binomial:$binomial, trying genus name:$genus\n";
		$completelineage = $complete_lineage_for_this_species{$genus};
		}


	print TIPS "$tip\t$completelineage\n";
	};
close TIPS;


print "
printed file list_of_terminals_in_your_tree (there are " , scalar @terminal_array , ")


newick string stored, 
	count of polytomies:$polytomies
	count of internal nodes:$interal_node
	count terminals:$count_terminals
	branchlengths stored:$bls_read
	newick_string:$newick_string
";

if($multiply_all_branchlengths >= 1.1)
	{
	print "\nNOTE! you have opted to multiply all branchlengths (by * $multiply_all_branchlengths)\n";
	print "\tthus $bls_multiplied branchlengths have been multiplied\n";
	};



unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree2




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################







#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub traverse_tree
{

# what is input:
# 0        1               2               3              4             5
# $test , $current_node , $assign_new_x , $assign_new_y, $node_colour, $cartoonize_clade

my $current_node = $_[0];
my $from_parent = $_[1];
my $sum_branchlength = $_[2]; 
my $node_plot_y = $_[3]; 	#  The previous node's assign_new_y; or simply 1 if the root node
my $node_colour = $_[4];
my $cartoonize_clade = $_[5];
my $count_connections = $child_counts{$current_node};
$subtree_newick_string = "($current_node)";
$nodes_traversed++;




if(exists($paint_clades{$current_node})){$node_colour = $paint_clades{$current_node}};





if($color_lineages_from_file == 1)
	{
	if($color_these_branches{$current_node} =~ /\w/)
		{
		$node_colour = $color_these_branches{$current_node};
		}else{
		$node_colour = $default_node_color;
		};
	
	
	};





if($cartoonize_clade == 1){$cartoonize_clade = 2};

if(exists($cartoonize_clades{$current_node}))
	{
#	$node_colour = $paint_clades{$current_node}
	$min_distance_to_tip_for_subtree = 999999;$max_distance_to_tip_for_subtree = 0;
	##################################################
	traverse_subtree($current_node, $from_parent, 0);#
	##################################################
	print "min_distance_to_tip_for_subtree:$min_distance_to_tip_for_subtree
max_distance_to_tip_for_subtree:$max_distance_to_tip_for_subtree
";
	$cartoonize_clade = 1;
	};







# print "\nNEW NODE ($current_node), from_parent:$from_parent count connecting nodes:$count_connections\n";



# defined child nodes, these are all connecting nodes except that from direction of origin in traversal.
my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	unless($connecting_node =~ /[\w\d]/)
		{die "\n\nerror 330. NODE ($current_node), from_parent:$from_parent " , 
			" count connecting nodes:$count_connections no:$all_connections connecting_node:$connecting_node\n\n"};

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};

my @replace_array = @next1;
my %branchlengths_for_current_children = ();

if($bls_read >= 2)
	{
	my @node_with_bls;
	foreach my $child(@next1)
		{
		$default_branchlength = 0.00001;
		my $bl = $default_branchlength;
		if(exists($branchlengths{$current_node}{$child}))
			{$bl = $branchlengths{$current_node}{$child}}else{die "\ncant find branchlength for $current_node to $child \n"};
		push @node_with_bls, "$child:$bl";
		if(exists($nodes_to_tip_from_current{$child}))
			{
		#	print "have step count for $child\n";
			}else{
		#	print "\nwhy no step count for $child\n"
			};
	# $nodes_to_tip_from_current{$current_node} = $max_steps_to_tip;


		};
	@replace_array = @node_with_bls;
	};

my $internal_taxname = "_" . $taxonomic_node_IDs{$current_node};

my $internal_taxname_raw = $taxonomic_node_IDs{$current_node};
$internal_taxname_raw =~ s/\d+$//;

my $rank_of_internal_taxname = $taxonomic_node_IDs_rank{$current_node};
# print "rank_of_internal_taxname:$rank_of_internal_taxname ";


my $abreviate_nodeID = $current_node;$abreviate_nodeID =~ s/INTERNAL_NODE_/IN/;
my $count_termnials = "_" . $count_terminals_node_lable{$current_node};


my $combined_node_label = $combined_node_labels{$current_node};


my $join_the_child_nodes = join ',', @replace_array;
$swap_string = "($join_the_child_nodes)";
$new_newick_string =~ s/$current_node(\W)/$swap_string$1/;
# $new_newick_string2 =~ s/$current_node(\W)/$swap_string$current_node$1/;
$new_newick_string3 =~ s/$current_node(\W)/$swap_string$internal_taxname$1/;
$new_newick_string4 =~ s/$current_node(\W)/$swap_string$abreviate_nodeID$count_termnials$internal_taxname$1/;
$new_newick_string5 =~ s/$current_node(\W)/$swap_string$combined_node_label$1/;

if($rank_of_internal_taxname eq "family")
	{
	$new_newick_string_family =~ s/$current_node(\W)/$swap_string$internal_taxname_raw$1/
	}else{
	$new_newick_string_family =~ s/$current_node(\W)/$swap_string$1/
	};
if($rank_of_internal_taxname eq "subfamily")
	{
	$new_newick_string_subfamily =~ s/$current_node(\W)/$swap_string$internal_taxname_raw$1/
	}else{
	$new_newick_string_subfamily =~ s/$current_node(\W)/$swap_string$1/
	};


unless($combined_node_label =~ /../){$combined_node_label = "NA"};


# $sec_tree_suport
# $node_support{$nodeID}


if($join_the_child_nodes =~ /_2288/)	
	{
#	print "\ncurrent_node:$current_node from_parent:$from_parent sum_branchlength:$sum_branchlength\n" , 
#		"\tnode_plot_y:$node_plot_y color:$node_colour count_connections:$count_connections\n" ,
#		"\t\@next1:@next1\nctnl0:$count_terminals_node_lable{$next1[0]} ctnl1:$count_terminals_node_lable{$next1[1]} ctnl2:$count_terminals_node_lable{$next1[2]}\n";
#	die "\nquit 12345\n"
	};


# order child nodes by number of decendents ... should result in laderized tree.
my @next2 = @next1;

# bifurcation working:
if ($#next2 == 1)
{
if($count_terminals_node_lable{$next1[1]} >= $count_terminals_node_lable{$next1[0]})
	{
	$next2[0] = $next1[0];$next2[1] = $next1[1];
	}else{
	$next2[0] = $next1[1];$next2[1] = $next1[0];
	};
}else{

# quick hack for partial laderisation of trichotomies, dont have time to figure out how to do perfectly for any polytomies
if($count_terminals_node_lable{$next1[2]} >= $count_terminals_node_lable{$next1[0]} && 
	$count_terminals_node_lable{$next1[2]} >= $count_terminals_node_lable{$next1[1]}	)
	{
	# last one already biggest, do nothing.	
	}elsif($count_terminals_node_lable{$next1[1]} >= $count_terminals_node_lable{$next1[0]} && 
	$count_terminals_node_lable{$next1[1]} >= $count_terminals_node_lable{$next1[2]}	){
	# 1 is biggest, needs to go to last one
	$next2[2] = $next1[1];$next2[1] = $next1[2];	
	}else{
	# 0 is bigggest, needs to go last
	$next2[2] = $next1[0];$next2[0] = $next1[2];	

	};

};
# print "\@next2:@next2
# count_terminals_node_lable{next2[0]}:$count_terminals_node_lable{$next2[0]}
# count_terminals_node_lable{next2[1]}:$count_terminals_node_lable{$next2[1]}
# count_terminals_node_lable{next2[2]}:$count_terminals_node_lable{$next2[2]}
# \n";



###########################################3

# here  print subtrees code ....

$list_terminals_in_subtree = "";

if($write_subtrees == 1 && 
	$count_terminals_node_lable{$current_node} <= $write_subtrees_tip_limit && 
	$count_terminals_node_lable{$current_node} >= 10 )
	{
	open(NODE_DESCRIPTIONS, ">$outprefix.subtree.$current_node.tree.bayestraits_commands_file") || die "\nerror 940.\n";
	print NODE_DESCRIPTIONS "1\n1\n";

	};

if($write_subtrees == 1)
	{
	##################################################
	traverse_subtree($current_node, $from_parent, 0);#
	##################################################
	$subtree_newick_string =~ s/^\((.+\))\)$/$1;/;
	};

if($write_subtrees == 1 && 
	$count_terminals_node_lable{$current_node} <= $write_subtrees_tip_limit && 
	$count_terminals_node_lable{$current_node} >= 10 )
	{
	$subtrees_written++;
	if($subtrees_written =~ /00$/)
		{
		print "printing subtree ($subtrees_written written), node:$current_node has tips $count_terminals_node_lable{$current_node}, which is < user limit $write_subtrees_tip_limit\n";
		};
	open(SUBTREE, ">$outprefix.subtree.$current_node.tree") || die "\nerror 937\n";
	print SUBTREE "$subtree_newick_string\n";
	close SUBTREE;
	print NODE_DESCRIPTIONS "run\n";
	close NODE_DESCRIPTIONS;

# 	push @list_of_bayestraits_commands , 
my $print_command = "rm bayestraits_tree_input scriptOUT; perl ~/usr_scripts/bayestraitswrap.pl $outprefix.subtree.$current_node.tree; rm community_ecol.bayestraits.log.txt; /home/douglas/software/bayestraits/BayesTraits32bit/BayesTraits bayestraits_tree_input community_ecol.bayestraits < $outprefix.subtree.$current_node.tree.bayestraits_commands_file > unwanted_screenout; mv community_ecol.bayestraits.log.txt $outprefix.subtree.$current_node.OUT;\n";
	open(BAYESTRAITS_BASH_COMMANDS, ">>bayestraits_bash_commands") || die "";
	print BAYESTRAITS_BASH_COMMANDS $print_command;
	close BAYESTRAITS_BASH_COMMANDS;


	}; # if($write_subtrees == 1 &&


my $count_terminals_from_node = "NA";
if($count_terminals_node_lable{$current_node} =~ /\d/){$count_terminals_from_node = $count_terminals_node_lable{$current_node}};

my $taxon_assigned_to_node = "NA";
if($taxonomic_node_IDs{$current_node} =~ /\w{3,}/){$taxon_assigned_to_node = $taxonomic_node_IDs{$current_node}};

my $lineage_assigned_to_node = "NA";
if($lineages_assigned_to_nodes{$current_node} =~ /\w{3,}/){$lineage_assigned_to_node = $lineages_assigned_to_nodes{$current_node}};

if(
	$count_terminals_node_lable{$current_node} <= $write_subtrees_tip_limit && 
	$count_terminals_node_lable{$current_node} >= 1 )
	{
	print NODE_LOG "$current_node\t$sum_branchlength\t$count_terminals_from_node\t$taxon_assigned_to_node\t$lineage_assigned_to_node\t$list_terminals_in_subtree\n";
	};

###########################################3






for my $index(0 .. $#next2)
	{
	my $test = $next2[$index];
	my $terminals_from_child = $count_terminals_node_lable{$test};
#	print "$test $terminals_from_child\n";


	#################################################################

	# a few thing for tree plotting

	unless($terminals_from_child =~ /\d/){$terminals_from_child = 1};
#	print "count_terminals_node_lable{current_node}:$count_terminals_node_lable{$current_node} terminals_from_child:$terminals_from_child\n";
	my $terminals_from_parent = $count_terminals_node_lable{$current_node};



	# change in position on the Y axis, calculated from the number of terminals from each child node
		#                                           PARENT                     CHILD
	my $delta_y = (($count_terminals_node_lable{$current_node} - 0) - $terminals_from_child ) / 2; # scalar @next1;


	my $assign_new_y;

#	if($#next2 == 1)
#	{
#	if($index == 0)# figure out later how to plot polytomies!
#		{
#		$assign_new_y = $node_plot_y+$delta_y; 
#		}else{$assign_new_y = $node_plot_y - $delta_y;
#		};
#	};

#	if($from_parent eq "New_Root")
#		{
#		print "index $index of $#next2, terminals_from_parent:$terminals_from_parent terminals_from_child:$terminals_from_child\n";

		# calculate positions for polytomies ...
		my $current_sum = 0; my $distance;
		for my $index_again(0 .. $index)
			{
			my $test_again = $next2[$index_again];
			my $terminals_from_child_again = $count_terminals_node_lable{$test_again};
			unless($terminals_from_child_again =~ /\d/){$terminals_from_child_again = 1};

			my $halfway = $terminals_from_child_again / 2;
			$distance = $current_sum + $halfway;
		#	print "\tindex_again:$index_again of $#next2, count terminals:$terminals_from_child_again\n";
			$current_sum += $terminals_from_child_again;
			};
		$delta_y = ($terminals_from_parent / 2) - $distance;
		$assign_new_y = $node_plot_y + $delta_y;

#		};

	my $branchlength_to_child;
	my $count_nodes_to_tip_from_current = $nodes_to_tip_from_current{$test} + 1;
	my $remaining = 1 - $sum_branchlength;


	if($plot_trees_with_branchlengths == 1)
		{

		my $default_bl = 0.00001;
		my $bl = $default_bl;
		if(exists($branchlengths{$current_node}{$test}))
			{$bl = $branchlengths{$current_node}{$test}};
		$branchlength_to_child = $bl;
		
		}else{
		if($count_nodes_to_tip_from_current >= 1)
			{
			$branchlength_to_child = $remaining / $count_nodes_to_tip_from_current;
			}else{
			$branchlength_to_child = 0.01;
			};

		};


	my $assign_new_x = $sum_branchlength+$branchlength_to_child;

	#################################################################


	if($plot_tree == 1)
		{
		###################################################################################################################################################
			# Variable names when inside following sub (very similar):
			#                   0      1             2               3       4             5                   6               7               8                   9                   10             11                  12                                13                              14                      15
			#                   test   current_node  scalar_next1    index   node_plot_y   sum_branchlength    assign_new_x    assign_new_y    abreviate_nodeID    internal_taxname    node_colour    cartoonize_clade    min_distance_to_tip_for_subtree   max_distance_to_tip_for_subtree   combined_node_label, terminals_from_child
		get_tree_plotting_commands($test, $current_node, scalar @next1, $index, $node_plot_y, $sum_branchlength , $assign_new_x , $assign_new_y , $abreviate_nodeID , $internal_taxname , $node_colour , $cartoonize_clade , $min_distance_to_tip_for_subtree, $max_distance_to_tip_for_subtree, $combined_node_label, $terminals_from_child);#
		###################################################################################################################################################
		};




	if($test =~ /^INTERNAL_NODE_/)
		{

	#	unless(exists($collapse_nodes{$test}))
	#		{
		#####################################
		traverse_tree($test , $current_node , $assign_new_x , $assign_new_y, $node_colour, $cartoonize_clade);#	# recurse
		#####################################
	#		};
		}else{

		# looks like code for printing all terminals:
#		my $R_command = "par(srt=0)\ntext($assign_new_x,$assign_new_y," . 
#			"adj=c(0, 0.5),labels=\"$test\",cex = 4)\n";
#		$draw_tree_R_commands_RECTANGL .= $R_command; 

		# probably too many terminals to print all,
		# inserted code from ITOL plotting script, prints most frequent tax for user specified number of tip labels
		$store_terminals_at_each_y_position{$assign_new_y} = $test;


		};





	}; # for my $index(0 .. $#next2)



# die "";


return();

};






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub store_nodes
{

my %rank_hash;

# nodes.dmp contains each taxon, described in a tree-like structure. 
# so each node has a name (eg polyphaga) and the higher group node to which it belongs (beetles). 

open(NODES , "nodes.dmp") || die "cant open nodes.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";

print "\nreading files from NCBI taxonomy database .... nodes.dmp .... ";


my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

			$ncbi_nodes{$tax_id}{rank} = $rank;
		#	$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode;
			$ncbi_nodes{$tax_id}{parent} = $parent_tax_id;
			$ncbi_nodes{$parent_tax_id}{child_nodes} .= "\t" . $tax_id;

		}else{
		print "line_counter:$line_counter line:$line";
		die "UNEXPECTED LINE:$line\nquitting\n";
		}
	$line_counter++;
	}

close(NODES);

#my @ranks = keys %rank_hash;@ranks = sort @ranks;
#print "ranks found in nodes.dmp:\n";
#print LOG "ranks found in nodes.dmp:\n";

#foreach my $rank(@ranks){print "$rank\t" , $rank_hash{$rank} , "\n";print LOG "$rank\t" , $rank_hash{$rank} , "\n"};

my @all_nodes = keys %ncbi_nodes;@all_nodes = sort @all_nodes;

print scalar @all_nodes , " nodes have been read.\n";


}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , "names.dmp") || die "cant open names.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";


print "\nnames.dmp, parsing 'scientific name', ignoring others ... ";

my $names_line_counter=0;
while (my $line = <NAMES>)
	{
# 24	|	Shewanella putrefaciens	|		|	scientific name	|

	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]*)\t[\|]\tscientific name/)
		{
		my $tax_id = $1;my $name = $2;#my $rank = $3;
		# print "tax_id:$tax_id name:$name\n";

		# if you want to remove non-alphanumerical characters from assigned species names:
		$name =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$name =~ s/\s\s+/ /g;$name =~ s/\s+$//;$name =~ s/^\s+//;
		$ncbi_nodes{$tax_id}{name} = $name;

		$names_line_counter++;#print "$names_line_counter\n";


		}else{
		if($line =~ /^(\d+).+scientific name/){die "UNEXPECTED LINE:\n$line\nquitting\n"}
		}

	}

close(NAMES);

print "$names_line_counter names parsed.\n";




}



#####################################################################################################
#
#
#
#####################################################################################################




sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated

my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);



if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	}



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;

	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

	if($name_string =~ /Zorochros/)
		{
		#print "name_string:($name_string) child:$child complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";
		};


	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;

	my $name_assignment_to_taxnumber = "";

	if($ignore_subspecies == 1)
		{
		if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
			{
			my $parentoriginalname = $ncbi_nodes{$current_node}{name};
			if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
				{
				$parentoriginalname =~ s/[\s\t]/_/g;
				$name_assignment_to_taxnumber = $parentoriginalname;
			#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
			#	print "\tassigning parent name instead:$parentoriginalname\n\n";
				}else{$name_assignment_to_taxnumber = $originalname}
			}else{
			$name_assignment_to_taxnumber = $originalname
			}

		}else{
		$name_assignment_to_taxnumber = $originalname
		}

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
	}}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################





sub proccess_tree
{
my $tree = shift;

unless($tree =~ /\(\w/){die "\nerror, sub process_tree called with no tree.\n"}

# print "tree:$tree\n";






##############################################################################################################


# if selecting one of these 3, also need to select corresponding string replace regex which is further below


# tree with binomials:
 while($tree =~ /([\(\,])([A-Z][a-z][a-z]+_[a-z][a-z][a-z]+)([\)\,\:])/){
# deleted close braket from first block in regex, that would be a node label not a terminal
# also, i cant see there would be name followed by open-bracket

# hack to work on tree with just genus names:
# while($tree =~ /([\(\,])([A-Z][a-z]+)([\(\)\,\:])/){print "warning hack (for genus only labels) apllied ... \n";

# this format:
# 3890,Simulium_verecundum_temp12.3:0.00028566965298361722):0.00
# while($tree =~ /([\(\,])([A-Z][a-z]+_[a-z][a-z][a-z]+_[a-z0-9\.]+)([\(\)\,\:])/){print "";

	
	my $before = $1; my $binomial  = $2; my $after = $3;
	$binomials_found_in_newick++; #print "\n";

	if($binomials_found_in_newick > $count_terminals*2)
		{
		print "\n\nerror, unexpected number of tips encountered. count_terminals:$count_terminals binomials_found_in_newick:$binomials_found_in_newick\n";
		print "printing outputs anyway for debug\n";
		return($tree);

		};

	my $add_tax_string = "NA_";
	my $look_for_ID = $binomial;

	my $completelineage;my $genus;
	if(exists($complete_lineage_for_this_species{$look_for_ID}))
		{
		$completelineage = $complete_lineage_for_this_species{$look_for_ID};
		}else{
		$genus = $look_for_ID;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
	#	print "no sucess looking for binomial:$binomial, trying genus name:$genus\n";
		$completelineage = $complete_lineage_for_this_species{$genus};
		}

# $binomial $genus

	if($completelineage =~ /\w/)
		{
		$add_tax_string = "";my $tax_added = 0;

		foreach my $current_rank (@ranksarray)
			{
			if($completelineage =~ /($current_rank):(\S+)/i)
				{
				my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
				$rank_assignment_counts{$rank}++;

				if($truncate_new_taxonomic_strings == 1)
					{
					$tax =~ s/^(\w\w\w\w\w\w)\w+/$1/;
					};
				if($uppercase_tax_strings == 1)
					{$tax  =uc($tax)};

				$add_tax_string .= "$tax" . "_";$tax_added++;
				};
			};

		unless($tax_added >= 1)
			{$add_tax_string = "NA_"; 
			if($warning_print_limit_1 < 50)
				{
			print "warning 1759, maybe inapproriate taxonomic ranks selected, because non found in current lineage\n";
			print "\tcompletelineage:$completelineage\n";
				}elsif($warning_print_limit_1 == 50)
				{print "warning 1759 print limit reached\n"};
			$warning_print_limit_1++;
			};

		#$add_tax_string =~ s/_$//;
		$lineages_assigned++;

		unless($add_tax_string =~ /\w/)
			{
		#	foreach my $current_rank (@ranksarray)
		#		{print "\tcurrent_rank:$current_rank\n"};
		#	die "\n\nwhat happened? completelineage:$completelineage\n";
			};

		}else{
		$newick_binomials_for_which_taxonomies_not_found++;	
		print "no lineage found for binomial:$binomial or genus:$genus\n";
		};


		# strict binomial:
		 if($tree =~ s/([\(\,])([A-Z][a-z][a-z]+_[a-z][a-z][a-z]+)([\)\,\:])/$1$add_tax_string$2$3/)
		# genus-only label hack:
		# if($tree =~ s/([\(\,])([A-Z][a-z][a-z]+)([\)\,\:])/$1$add_tax_string$2$3/)
		# general:
 		# if($tree =~ s/([\(\,])([A-Z][a-z]+_[a-z][a-z][a-z]+_[a-z0-9\.]+)([\(\)\,\:])/$1$add_tax_string$2$3/)
			{
			my $a = $1; my $b = $2; my $c = $3;

			if($binomials_found_in_newick =~ /000$/)
				{
				print "\nproccesed $binomials_found_in_newick binomial TIPs in newick\n" , 
				"\tbefore:$before binom:$binomial after:$after\n" , 
				"\ta:$a add_tax_string:$add_tax_string b:$b c:$c\n";
				};

			
			}else{
			die "\n\nerror 767, cant excise tip ... before:$before binom:$binomial after:$after\n\n"
			};

## hack part 2:
#		$tree =~ s/([\(\,])([A-Z][a-z]+)([\(\)\,\:])/$1$add_tax_string$2$3/;



	}; #  while($tree =~



print "
binomials_found_in_newick:$binomials_found_in_newick
of these, lineages_assigned:$lineages_assigned
newick_binomials_for_which_taxonomies_not_found:$newick_binomials_for_which_taxonomies_not_found
";

unless($lineages_assigned >= 1)
	{
	print "seems error, no lineages found for binomials of your tree\n";

	if($tree =~ /\([A-Z][a-z]{5,14}\:1.0\,/)
		{
		print "\tprobably cause you have only species names,\n\t" , 
			"if required, check the script, there is a hack-around\n";
		};

	};

return($tree);




}


##############################################################################################################




##############################################################################################################






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

# perl ~/usr_scripts/process_newick.pl -intree inNUCL.RMmulti.FtBl -outtree inNUCL.RMmulti.FtBl.tax -ncbi_tax_node 6960


print "

process_newick.pl
	processes newick tree labels, 
	adds taxonomic information to terminals and internal nodes,
	re-roots
	can plot tree

syntax:
perl process_newick.pl -intree tree.nwk -outtree tree2.nwk -ncbi_tax_node 33208

";





if($arguments =~ /-intree\s+(\S+)/)
	{
	$intree = $1;print "intree:$intree\n";
	}else{
	die "\ncommand error\n";
	};



if($arguments =~ /-combine_node_labels\s+(\S+)\s+(INTERNAL[A-Z\d\_]+)/)
	{
	$combine_node_label_tree = $1; $second_tree_reroot_node = $2; $combine_node_labels = 1;
	print "combine_node_labels\n";
	}elsif($arguments =~ /-combine_node_labels\s+(\S+)/)
	{
	$combine_node_label_tree = $1; $combine_node_labels = 1;
	};



if($arguments =~ /-outprefix\s+(\S+)/)
	{
	$outprefix = $1;print "outprefix:$outprefix\n";
	}else{
	die "\ncommand error\n";
	};

$sonify = 0;$sonify_file;
if($arguments =~ /-sonify\s+(\S+)/) # intree already specified, here give the output file. 
	{
	$sonify_file = $1; $sonify = 1;
	};


if( $arguments =~ /-ncbi_tax_node\s+(\d+)/ )
	{
	$starting_node = $1;print "starting_node:$starting_node\n";
	}else{
	print "user did not given NCBI taxonomic number. using default 33208 (Metazoa)
if you have references outside this default, you need to get the appropriate NCBI taxonomy number from http://www.ncbi.nlm.nih.gov/taxonomy/
and input this here with the option -node taxon_number
";
	$starting_node = 33208;
	}


if($arguments =~ /-plot_tree/)
	{
	$plot_tree	= 1;print "plot_tree\n";
	}

if($arguments =~ /-multiply_all_branchlengths\s+(\d+)/)
	{
	$multiply_all_branchlengths = $1;print "multiply_all_branchlengths:$multiply_all_branchlengths\n";
	}

if($arguments =~ /-highlight_tip_string_match\s+(\S+)/)
	{
	$tip_string_to_highlight = $1; # "CK_OTU";
	$highlight_tip_string_match = 1; # find tips with string below, then draw circles at these
	print "highlight_tip_string_match:$tip_string_to_highlight\n";
	};



# -piechart_trait_prob BOLD_speciestraits2.U_ITD BayesTraits.tip_predictions.U_ITD.0.5.out2
if($arguments =~ /-piechart_trait_prob\s+(\S+)\s+(\S+)/)
	{
	$traitfileA = $1;$traitfileB = $2;
	$piechart_trait_prob = 1; 
	print "\nuser specified to make piecharts from trait probabilities\n";

	#####################################################################
	# this file gives known states
	open(IN_TRAIT_A, $traitfileA) || die "\nerror 2266. cant open file ($traitfileA).\n";
	while(my $line = <IN_TRAIT_A>)
		{
		# OTU_HCL0293	?
		# OTU_HCL0295	?
		# OUTGROUP	-
		# Osmia_aglaia	C
		$line =~ s/\n//;$line =~ s/\r//;
		if($line =~ /^(\S+)\t([A-Z])/) # used to have [A-Z]+, but this picks up where multiple states, and complicates things
			{my $sp = $1; my $st = $2; $known_state{$sp} = $st;$all_states_for_pie{$st} = 1 };
		};
	close IN_TRAIT_A;
# must sort following, because pie plotting function default sorts, thus mixing colors if this object unsorted. 
#  will need to check works ok when irregular characters used
@all_states_for_pie_array = keys %all_states_for_pie;@all_states_for_pie_array = sort @all_states_for_pie_array;
	#####################################################################
	# this file gives predicted states.
	open(IN_TRAIT_B, $traitfileB) || die "\nerror 2283. cant open file ($traitfileB).\n";
	while(my $line = <IN_TRAIT_B>)
		{
		$line =~ s/\n//;$line =~ s/\r//; # print "$line\n";
		if($line =~ /^(\S+)\t(.+)/)
			{
			my $species = $1; my $state_probs = $2;my @split_state_probs = split /\t+/ , $state_probs;$predicted_state{$species} = "A";
			foreach my $state(@split_state_probs)
				{
				if($state =~ /^([A-Z])\:([0-9\.]+)/)
					{
					my $st = $1; my $pr = $2; $state_probabilities_for_pie{$species}{$st} = $pr; # print "\t$s,$p\n"
					};
				};			
			};
		};
	close IN_TRAIT_B;
	#####################################################################
	};# if($arguments =~ /-piechart_trait_prob\s+(\S+)\s+(\S+)/)

# $known_state{$sp}; $predicted_state{$sp} = "A"; $state_probabilities_for_pie{$sp}{$pr} = $p;


if($arguments =~ /-remove_polytomies/)
	{
	$remove_polytomies	= 1;print "remove_polytomies\n";
	}

$write_subtrees = 0;
if($arguments =~ /-subtrees\s+(\d+)/)
	{
	$write_subtrees_tip_limit = $1;
	$write_subtrees = 1;

	open(BAYESTRAITS_BASH_COMMANDS, ">bayestraits_bash_commands") || die "";
	close BAYESTRAITS_BASH_COMMANDS;
	};

if($arguments =~ /-print_terminal_tax_labels\s+(\d+)/)
	{
	$num_labs_to_print = $1;
	$print_terminal_tax_labels = 1;
	print "print_terminal_tax_labels:$num_labs_to_print\n";
	};


if($arguments =~ /-replace_scientific_notation_in_branchlengths/)
	{
	$replace_scientific_notation_in_branchlengths = 1;
	print "replace_scientific_notation_in_branchlengths\n";
	}

if($arguments =~ /-plot_circular_tree_with_branchlengths/)
	{
	die "\nerror, command depracated: -plot_circular_tree_with_branchlengths, use instead -plot_trees_with_branchlengths\n";
	};

if($arguments =~ /-plot_trees_with_branchlengths/) # previously plot_circular_tree_with_branchlengths, though should apply to rectangular also
	{
	$plot_trees_with_branchlengths = 1;
	print "plot_trees_with_branchlengths\n";
	};

if($arguments =~ /-plot_cartoonized_clades/)
	{
	$plot_cartoonized_clades = 1;
	print "plot_cartoonized_clades\n";
	};

if($arguments =~ /-print_internal_labels_circular_tree/)
	{
	$print_internal_labels_circular_tree = 1;
	print "print_internal_labels_circular_tree\n";
	};

if($arguments =~ /-internal_names_STA\s+(.+)\s+-internal_names_FIN/)
	{
	my $internal_names_circulartree = $1;print "internal_names_circulartree:$internal_names_circulartree\n";
	my @splitem = split /\s+/ , $internal_names_circulartree;
	foreach my $splititems(@splitem)
		{
		if($splititems =~ /(.+)\:(.+)/)
			{
			my $nodecode = $1; my $printthis = $2;$internal_names_circular_tree{$nodecode} = $printthis;
			};
		};
	};



if($arguments =~ /-identify_nodes_via_decendents\s+(\S+)/)
	{
	$identify_nodes_via_decendents_file = $1;
	$identify_nodes_via_decendents = 1;
	print "\nyou have opted to identify certain nodes in the plotted tree, via lists of their decendents.\n";
	print "looking for decendent lists in the specified file:$identify_nodes_via_decendents_file\n";



	open(INFILE, $identify_nodes_via_decendents_file) || die "";
	my $line_count92 = 0;
	while (my $line  = <INFILE>)
		{
		$line =~ s/\n//;$line =~ s/\r//;$line =~ s/\s+$//;
		my @splitline = split /\s+/ , $line;
		# print "line count:$line_count92 node has decendenets:$#splitline\n";

		if($#splitline >= 2)
			{
			@splitline = sort @splitline;#print "\@splitline:@splitline\n";
			my $join_sorted_line = join ' ' , @splitline;
			$store_decendent_defined_nodes{$join_sorted_line} = 1;
			$decendent_defined_nodes_stored++;
			}else{
			print "warning, node specified only has $#splitline decendetns. ignoring\n";
			};

		$line_count92++;
		};
	close INFILE;
	unless($line_count92 >= 1){die "\nerror, nothing defined in your file ....\n"};

	unless($decendent_defined_nodes_stored >= 1){die "no significant nodes found in file.\n"};
	print "file read, has $line_count92 lines. $decendent_defined_nodes_stored nodes will be looked for in tree.
\n";

	};





if($arguments =~ /-label_tips\s+(\S+)\s+(\S+)/)
	{
	$tip_label_string = $1;$tip_label_file = $2;
	print "you have opted to add the string:$tip_label_string " , 
		"to terminal of the tree of members from the file:$tip_label_file\n";
	open(IN, $tip_label_file) || die "\nerror 1095\n";
	while (my $line = <IN>)
		{
	#	print $line;
		$line =~ s/\n//;$line =~ s/\r//;
		if($line =~ /^>(.+)/){$add_string_to_this_member{$1} = 1;$members_to_add_string_to++};
		};
	close IN;

	print "will add string to $members_to_add_string_to, assuming these are on the tree\n";

	}else{
	};







#$output_filename	= "$treefile.query_clades";

print "\n";

}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub assign_names_to_internal_nodes
{
my $current_node = $_[0];my $from_parent = $_[1];
my $count_connections = $child_counts{$current_node};
$count_sub_calls++;

# print "\nsub assign_names\n\tcurrent_node:$current_node from_parent:$from_parent\n";



# defined child nodes, these are all connecting nodes except that from direction of origin in traversal.
my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};

if($current_node =~ /INTERNAL_NODE/)
	{
	$terminals_belonging_to_current_node = "";
	$terminals_belonging_to_current_node_which_have_tax_info = "";
	$max_steps_to_tip = 0;

	###############################################################
	get_terminals_from_this_node($current_node , $from_parent , 0);#
	###############################################################

	my $shared_tobycode = "";my $shared_rank = ""; my $shared_lineage = ""; my $sub_results;
	# particularly with OTUs in the data, there will be cases where a node has lots of terminals,
	# but, one or even none have taxonomic information. no point calling get_shared for these.
	my $count_those_with_tax_info = count_how_many_terminals_have_lineage_information($terminals_belonging_to_current_node);
#	print "terminals_belonging_to_current_node:$terminals_belonging_to_current_node\n";
#	print "\ncount_those_with_tax_info:$count_those_with_tax_info";


	if($count_those_with_tax_info >= 2)
		{
		
		##################################################
		$sub_results = get_shared_taxonomic_name($terminals_belonging_to_current_node_which_have_tax_info);#
		##################################################

# returnstring = $most_inclusive_name . "	" .  $rank_of_most_inclusive_name . "	" . $most_inclusive_lineage;


		}else{
		};

	$terminals_belonging_to_current_node =~ s/(\t)\t+/$1/;
	$terminals_belonging_to_current_node =~ s/^\t+//;$terminals_belonging_to_current_node =~ s/\t+$//;
	my @count_terms_array = split /\t/ , $terminals_belonging_to_current_node;
	my $count_termnials = scalar @count_terms_array;


	if($identify_nodes_via_decendents == 1)
		{
		my @sorted_terminals_array = sort @count_terms_array;
		my $join_sorted_line = join ' ' , @sorted_terminals_array;
		if($store_decendent_defined_nodes{$join_sorted_line} == 1)
			{
			$decendent_defined_nodes_found_in_tree++;
			$label_these_decendent_defined_nodes_found{$current_node} = 1;
			};
		};

	if($combine_node_labels == 1)
		{
		my @sorted_terminals_array = sort @count_terms_array;
		my $join_sorted_line = join ' ' , @sorted_terminals_array;
			
		if($store_second_tree_support{$join_sorted_line} =~ /\d/)
			{
			my $support_on_equivelent_node_of_other_tree = $store_second_tree_support{$join_sorted_line};
			my $support_ths_tree = $node_support{$current_node};


			$combined_node_labels{$current_node} = $support_ths_tree . "_" . $support_on_equivelent_node_of_other_tree;

		#	print "matched node. Support other tree:$support_on_equivelent_node_of_other_tree, this tree:$support_ths_tree\n";
			}elsif(length($join_sorted_line) <= 600)
				{
				print "cant match node $current_node, termnials:$join_sorted_line\n";
				};
		};


	if($write_subtrees == 1)
		{
		$store_terminals_derived_from_this_node{$current_node} = join ' ' , @count_terms_array
		};


	if($sub_results =~ /^([^\t]+)\t(.+)\t(.+)/){$shared_tobycode = $1; $shared_rank = $2; $shared_lineage = $3};
	print LOG "$current_node\t$count_termnials\t$shared_tobycode\t$shared_lineage\n";


	if($count_sub_calls =~ /000$/)
		{
		print "\nsub assign_names called $count_sub_calls times\n\tassigning name to NODE($current_node),\n" , 
		"\tfrom_parent:$from_parent	count_those_with_tax_info:$count_those_with_tax_info\n" , 
		"\tcount connecting nodes:$count_connections count_termnials:$count_termnials\n",
		 "\tterminals_belonging_to_current_node_which_have_tax_info:$terminals_belonging_to_current_node_which_have_tax_info\n",
		 "\tassigning name ($shared_tobycode) to NEW NODE ($current_node), from_parent:$from_parent\n" , 
		"\tcount connecting nodes:$count_connections count_termnials:$count_termnials\n";
		};


	#print "";

#	if($count_those_with_tax_info >= 2 && $count_those_with_tax_info <=5 && $shared_tobycode =~ /Hymenoptera/){die ""};

#	if($count_those_with_tax_info >= 2){
#	unless($shared_tobycode =~ /\w/){die ""};
#		};

	if($shared_tobycode =~ /\w/)
		{$internal_nodes_with_tax_assigned++;
		}else{
	############	die "\nno tax name returned\n"
		};

	my $assigned_name = $shared_tobycode;
	my $increment =1;
	while (exists($all_names_assigned{$assigned_name}))
		{$increment += 1;
		my $try_name = $shared_tobycode . $increment;
		$assigned_name = $try_name;
		};
	$all_names_assigned{$assigned_name} = 1;
	$taxonomic_node_IDs{$current_node} = $assigned_name;
	$taxonomic_node_IDs_rank{$current_node} = $shared_rank;



	$lineages_assigned_to_nodes{$current_node} = $shared_lineage;

	$count_terminals_node_lable{$current_node} = $count_termnials;
	$nodes_to_tip_from_current{$current_node} = $max_steps_to_tip;


	}else{
	$nodes_to_tip_from_current{$current_node} = 0;
	}; # if($current_node =~ /INTERNAL_NODE/)




for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		assign_names_to_internal_nodes($test , $current_node );#	# recurse
		#####################################
		}
	}

return();

}

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub get_terminals_from_this_node
{
my $current_node = $_[0];my $from_parent = $_[1]; my $steps_to_tip = $_[2];$steps_to_tip++;
my $count_connections = $child_counts{$current_node};

my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};


for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
	#	unless(exists($collapse_nodes{$test}))	
	#		{
		#####################################
		get_terminals_from_this_node($test , $current_node , $steps_to_tip);#	# recurse
		#####################################
	#		};
		}elsif($test =~ /[\w\d]/)# 20180307: some suspected errors here, thus require a name
		{

		if($steps_to_tip > $max_steps_to_tip){$max_steps_to_tip = $steps_to_tip};

	#	print "CN:$current_node FP:$from_parent steps_to_tip:$steps_to_tip max:$max_steps_to_tip index:$index test:$test\n";

		# store all terminals:
		$terminals_belonging_to_current_node .= "$test\t";

		# store only terminal with tax info:
		my $test_species = $test;$test_species =~ s/^([A-Z][a-z]+)_.+/$1/;
		my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
		if($test_taxonomy =~ /\w+/)
			{$terminals_belonging_to_current_node_which_have_tax_info .= "$test\t"};

		};
	}




return();

};#sub get_terminals_from_this_node


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_shared_taxonomic_name
{
my $tax_list = $_[0];# my $members_with_tax_info = $_[1];


# jan2017: some changes made to this long-used subroutine

if($tax_list =~ s/(\t)\t+/$1/g){print "\nwhy missing entries in this object?\n"};
my @tax_array = split /\t/ , $tax_list;
fisher_yates_shuffle( \@tax_array );

if($#tax_array == 0){die "\nwhy sub get_shared_taxonomic_name called when theres only 1?\n"};

my $shared_tax_substring = "";

if($verbose >= 1){
print "\nfinding largest common substring  ....\n\n";
}


my $index=0;

# my $test_species = $tax_array[0];
# my $test_taxonomy = $complete_lineage_for_this_species{$test_species};

my $test_taxonomy;
my $x = 0; my $j = 0;
while($x == 0)
	{
	my $test_species = $tax_array[$j];		
	$test_taxonomy = $complete_lineage_for_this_species{$test_species};
	if($test_taxonomy =~ /\w+/)
		{
		unless($test_taxonomy =~ / order\:$test_species $/)# 20170320: ignore lineage information up to order only
			{$x = 1};
	
		}else{
		$test_species =~ s/^([A-Z][a-z]+)_.+/$1/;# try genus only ...
		$test_taxonomy = $complete_lineage_for_this_species{$test_species};
		if($test_taxonomy =~ /\w+/)
			{
			unless($test_taxonomy =~ / order\:$test_species $/)# 20170320: ignore lineage information up to order only
				{$x = 1};

			}else{
			# should keep loop ...
			# $test_taxonomy = " no_rank:cellular_organisms genus:$test_species";#print "\nerror 770, no taxonomy found for this:($test_species)\n";
			};
		};	
#	print "finding initial lineage for comparisons x:$x j:$j test_species:$test_species test_taxonomy:$test_taxonomy\n";

	$j++;
	};

if($x == 0 ){die "\nerror 1118. could not find any lineage.\n"};
unless($test_taxonomy =~ /\w+/){die "\nerror 1119. could not find any lineage.\n"};


# print "test_taxonomy:$test_taxonomy\n";
my $most_inclusive_name = "";
my $rank_of_most_inclusive_name = "";
my $most_inclusive_lineage = "";


while($test_taxonomy =~ s/^\s*([^:]+):(\w+)//)
	{
	my $current_rank = $1;my $current_taxname = $2;

	my $all_members_have_this_name =1;
	my $end;

	if($#tax_array >= 200)
		{
		$end = $#tax_array*$reduce_by_proportion; $end =~ s/\.\d+//;
		}elsif($#tax_array >= 30){
		$end = $j+$reduced_comparisons;
		}else{
		$end = $#tax_array; 
		};	

# taking sample or propirtion based on how many terminals does not always work,
# try something else, sample based on count of how many have tax info
# 


	my $found_comparisons =0;
	foreach my $index($j .. $end)
		{
		$tax = $tax_array[$index];

		$tax =~ s/^([A-Z][a-z]+)\_.+$/$1/;
		my $test_taxonomy2 = $complete_lineage_for_this_species{$tax};

		unless($test_taxonomy2 =~ / order\:$tax $/)# 20170320: ignore lineage information up to order only
			{

			if($test_taxonomy2 =~ /\w+/)
				{
				$found_comparisons = 1;
				if($test_taxonomy2 =~ /\:$current_taxname\s/)
					{}else{
					$all_members_have_this_name =0;#print "0";
					};
				};

			};

	#	print "\t\tindex:$index tax:$tax all_members_have_this_name:$all_members_have_this_name\n";
		};

	unless($found_comparisons == 1){
		print "\nWARNING 1153, did not find any (of $j to $end;" , 
			" \$\#tax_array:$#tax_array) with lineages.\n"};

	if($all_members_have_this_name == 1)
		{
		$most_inclusive_name = $current_taxname;$rank_of_most_inclusive_name = $current_rank;
		$most_inclusive_lineage .= "$current_rank:$current_taxname ";
		};


	}



my $returnstring = $most_inclusive_name . "	" .  $rank_of_most_inclusive_name . "	" . $most_inclusive_lineage;
return($returnstring);


}





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


	 sub fisher_yates_shuffle 	# http://perldoc.perl.org/perlfaq4.html
		{
	 	my $deck = shift;my $i = @$deck;
 		while (--$i) 
			{my $j = int rand ($i+1);
 			@$deck[$i,$j] = @$deck[$j,$i]}
 		}





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub count_how_many_terminals_have_lineage_information
{
my $tax_list = shift;

if($tax_list =~ s/(\t)\t+/$1/g){print "\nwhy missing entries in this object?\n"};
my @tax_array = split /\t/ , $tax_list;

my $count_number_with_tax_info = 0;

for my $j(0 .. $#tax_array)
	{
	my $found=0;
	my $test_species = $tax_array[$j];		
#	$test_taxonomy = $complete_lineage_for_this_species{$test_species};
#	if($test_taxonomy =~ /\w+/)
#		{
#		$found = 1;
#		}else{
		$test_species =~ s/^([A-Z][a-z]+)_.+/$1/;# try genus only ...
		$test_taxonomy = $complete_lineage_for_this_species{$test_species};
		if($test_taxonomy =~ /\w+/)
			{
		#	if($test_species =~ /Hymenoptera/)
		#		{die "\ntest_species:$test_species test_taxonomy:$test_taxonomy\n"};

			unless($test_taxonomy =~ / order\:$test_species $/)# 20170320: ignore lineage information up to order only
				{$found = 1};
			}
#		};	

	if($found == 1)
		{
		$count_number_with_tax_info++;
		};

	};

return($count_number_with_tax_info);


}




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


# this is invoked for both child nodes of current node, tax name is for current node, not child

sub get_tree_plotting_commands
{
my $test = $_[0]; 	 # child node
my $current_node = $_[1]; # node
my $scalar_next1 = $_[2]; my $index = $_[3];
my $node_plot_y = $_[4]; my $sum_branchlength = $_[5];
my $assign_new_x = $_[6]; my $assign_new_y = $_[7];
my $abreviate_nodeID = $_[8];my $internal_taxname = $_[9];my $node_colour = $_[10];
my $cartoonize_clade  = $_[11]; my $min_distance_to_tip_for_subtree  = $_[12];
my $max_distance_to_tip_for_subtree  = $_[13]; my $combined_node_label =  $_[14];my $count_terminalsfromnode =  $_[15];

unless($count_terminalsfromnode =~ /\d/){die "\nerror 2986\n"};

#######################################################################################
#
#######################################################################################

# note $count_terminals only takes one value, count of all terminals in the backbone made earlier

my $y1_proportion = ($node_plot_y / $count_terminals)+1;	$y1_proportion *= $multipleier;
my $y2_proportion = ($assign_new_y / $count_terminals)+1;	$y2_proportion *= $multipleier;

# note $y2_proportion is ordered along y axis differently than output by .process_newick

#	my $y1_proportion = ($node_plot_y / $count_terminals);	$y1_proportion *= $multipleier;
#	my $y2_proportion = ($assign_new_y / $count_terminals);	$y2_proportion *= $multipleier;



# start position.
my $Xstart 	= $sum_branchlength * cos($y1_proportion);# X ax
my $Ystart 	= $sum_branchlength * sin($y1_proportion);# Y ax
# end position.
my $Xend 	= $assign_new_x * cos($y2_proportion);# X
my $Yend 	= $assign_new_x * sin($y2_proportion);# Y
# scale branch thickness:
my $width3 = $width;

# unless($node_colour =~ /grey/ || $node_colour =~ /gray/ || $node_colour eq "black"){$width3 = $big_clade_size};



if($scale_branch_thickness == 1)
	{
	$width3 = $width / ( $assign_new_x * $branch_root_to_tip_thickness_scale_factor);
	};

if($nodes_traversed =~ /0000$/)
	{
	print "nodes_traversed:$nodes_traversed child:$test parent:$current_node\n" , 
	"\tcount_nodes_to_tip_from_current:$count_nodes_to_tip_from_current X:$assign_new_x Y:$assign_new_y\n";
	print "\tnode plot Y:$node_plot_y coutn terms:$count_terminals\n";
	print "\ty1_proportion:$y1_proportion y2_proportion:$y2_proportion\n"
	};



############################################################
# rectangular: 		$sum_branchlength $assign_new_x $y1_proportion $y2_proportion

if($cartoonize_clade =~ /[12]/)
	{
	
	}else{

if($sonify == 1)
	{

	# probably best thing to record would be:
	#	for each node, x position of bifurcation, y positions of each child branch.
	# 	also child branch length, maybe just the shortest one.
	#	store according to node id.
	# 	sort according to y, omit nodes at exactly the same y.
	

#	$sonify_commands{$sum_branchlength  . "__"  .  $y1_proportion} = 1;
#	$sonify_commands{$sum_branchlength  . "__"  .  $y2_proportion} = 1;
#	$sonify_commands{$assign_new_x  . "__"  .  $y2_proportion} = 1;

	$sonify_node_x{$current_node} = $sum_branchlength;# print "sonify node x:$sum_branchlength\n";
	$sonify_node_y{$current_node} = $y1_proportion;

	unless($sonify_child_nodes{$current_node} =~ /\t$test\t/){$sonify_child_nodes{$current_node} .= "	$test	"};
	$sonify_child_xy{$current_node}{$test} = "$assign_new_x\t$y2_proportion";

	};# if($sonify == 1)


# for one child invoked for sub:::

# parent to one child, vertical:
my $R_command =  "segments(" . "$sum_branchlength, $y1_proportion , $sum_branchlength, $y2_proportion , col = \"$rect_branch_color_default\", lwd = $rectangular_plot_branch_width)\n";
$draw_tree_R_commands_RECTANGL .= $R_command;


########################################################
#                    first 5 columns  VIDEO, 6-9 audio     
 print VIDIFY "$branchnumber\t$sum_branchlength\t$y1_proportion\t$sum_branchlength\t$y2_proportion\t" , 
			"$sum_branchlength\t$y1_proportion\t$sum_branchlength\t$y1_proportion\t$count_terminalsfromnode\n";
$branchnumber++;
# print VIDIFY2 "$branchnumber\t$sum_branchlength\t$y1_proportion\t$sum_branchlength\t$y1_proportion\n";$branchnumber++;
########################################################


# parent to one child, horizontal
my $R_command =  "segments(" . "$sum_branchlength, $y2_proportion , $assign_new_x, $y2_proportion , col = \"$rect_branch_color_default\", lwd = $rectangular_plot_branch_width)\n";
$draw_tree_R_commands_RECTANGL .= $R_command;


########################################################
#                    first 5 columns  VIDEO, 6-9 audio     
 print VIDIFY "$branchnumber\t$sum_branchlength\t$y2_proportion\t$assign_new_x\t$y2_proportion\t" , 
		"$sum_branchlength\t$y1_proportion\t$assign_new_x\t$y1_proportion\t$count_terminalsfromnode\n";
$branchnumber++;
# print VIDIFY2 "$branchnumber\t$sum_branchlength\t$y1_proportion\t$assign_new_x\t$y1_proportion\n";$branchnumber++;
########################################################


	
	};



if($cartoonize_clade == 1)
	{
	# my $y1_proportion = ($node_plot_y / $count_terminals)+1;	$y1_proportion *= $multipleier;
	my $y_cartoon_upper = ( (($node_plot_y-0.5) + ($count_terminals_node_lable{$current_node}/2)) / $count_terminals ) +1;$y_cartoon_upper *= $multipleier;
	my $y_cartoon_lower = ( (($node_plot_y+0.5) - ($count_terminals_node_lable{$current_node}/2)) / $count_terminals ) +1;$y_cartoon_lower *= $multipleier;
	my $calc22 = $sum_branchlength+$min_distance_to_tip_for_subtree;
	my $calc33 = $sum_branchlength+$max_distance_to_tip_for_subtree;

	my $R_command =  "segments(" . "$sum_branchlength, $y1_proportion , $calc22, $y_cartoon_upper, col = \"black\", lwd = $rectangular_plot_branch_width)\n";	
	$draw_tree_R_commands_RECTANGL .= $R_command;
	my $R_command =  "segments(" . "$sum_branchlength, $y1_proportion , $calc22, $y_cartoon_lower, col = \"black\", lwd = $rectangular_plot_branch_width)\n";	
	$draw_tree_R_commands_RECTANGL .= $R_command;

	my $R_command =  "segments(" . "$calc22, $y_cartoon_upper, $calc33 , $y_cartoon_upper, col = \"black\", lwd = $rectangular_plot_branch_width)\n";	
	$draw_tree_R_commands_RECTANGL .= $R_command;
	my $R_command =  "segments(" . "$calc22, $y_cartoon_lower, $calc33 , $y_cartoon_lower, col = \"black\", lwd = $rectangular_plot_branch_width)\n";	
	$draw_tree_R_commands_RECTANGL .= $R_command;


	my $R_command = "par(srt=0)\ntext($calc33,$y1_proportion," . 
		"adj=c(0, 0.5),labels=\"$cartoonize_clades{$current_node}\",cex = $cartooned_clade_external_label_cex)\n";
	$draw_tree_R_commands_RECTANGL .= $R_command;


	};



if($test =~ /^INTERNAL_NODE_/)
	{
	if($node_support{$test} >= 0.75)
		{
		my $R_command = "points($assign_new_x,$y2_proportion,col=\"green\", pch=16, cex = $rectangular_plot_support_circle_size)\n";
		$draw_tree_R_commands_RECTANGL2 .= $R_command; 
		};

	if($combine_node_labels == 1)
		{
		print "combined_node_label:$combined_node_label\n";
		if($combined_node_label eq "NA")
			{
			
			}elsif($combined_node_label =~ /^_([\d\.]+)_([\d\.]+)$/){
				{
				my $support1 = $1;my $support2 = $2;

				my $R_command = "colors1<- rgb(colorfunc($support1) , maxColorValue=255)\n" . 
					 "points($assign_new_x,$y2_proportion,col=colors1, pch=16, cex = 2)\n";
				$draw_tree_R_commands_RECTANGL2 .= $R_command; 

				my $R_command = "colors2<- rgb(colorfunc($support2) , maxColorValue=255)\n" . 
					 "points($assign_new_x,$y2_proportion,col=colors2, pch=16, cex = 1)\n";
				$draw_tree_R_commands_RECTANGL2 .= $R_command; 


				};

			
			};
		
		};



	}else{# TIP, write tip name:

	######################################################
	my $prestring = return_tax_string_for_this($test);#
	######################################################

	my $plot_current_label_X = $assign_new_x + $plot_tip_label_offset_rectangular;
	my $R_command = "par(srt=0)\ntext($plot_current_label_X,$y2_proportion," . 
		"adj=c(0, 0.5),labels=\"$prestring$test\",cex = $rectangular_plot_tip_label_cex)\n";

	if($print_tips_rectangular_tree == 1)
		{$draw_tree_R_commands_RECTANGL .= $R_command};



	print RECT_PLOT_LOG "$prestring$test\t$plot_current_label_X\t$y2_proportion\n";




	};
############################################################


if($branch_plot_type == 1)
	{
	# single line branches:
	# tree BRANCHES,      x1, y1, x2, y2
	my $R_command =  "segments(" . "$Xstart, $Ystart , $Xend, $Yend , col = \"$node_colour\", lwd = $width3)\n";
	$draw_tree_R_commands .= $R_command;


#	if($internal_taxname =~ /Elateroidea/)
#		{

#	my $y1_proportion = ($node_plot_y / $count_terminals)+1;	$y1_proportion *= $multipleier;


my $upper = $node_plot_y + ($count_terminals_node_lable{$current_node} / $scalar_next1);
$upper = ($upper / $count_terminals)+1; $upper *= $multipleier;
			my $lower = $node_plot_y - ($count_terminals_node_lable{$current_node} / $scalar_next1);
$lower = ($lower / $count_terminals)+1; $lower *= $multipleier;


	my $Xstart 	= $sum_branchlength * cos($y1_proportion);# X ax
	my $Ystart 	= $sum_branchlength * sin($y1_proportion);# Y ax
	# end position.
		my $Xend 	= 1 * cos($upper);# X
		my $Yend 	= 1 * sin($upper);# Y
		my $R_command =  "segments(" . "$Xstart, $Ystart , $Xend, $Yend , col = \"$big_clade_color\", lwd = $big_clade_size)\n";
		$draw_tree_R_commands .= $R_command;
		my $Xend 	= 1 * cos($lower);# X
		my $Yend 	= 1 * sin($lower);# Y
		my $R_command =  "segments(" . "$Xstart, $Ystart , $Xend, $Yend , col = \"$big_clade_color\", lwd = $big_clade_size)\n";
		$draw_tree_R_commands .= $R_command;
			

		my $lower_val;my $higher_val;
		if($upper > $lower){$higher_val = $upper; $lower_val = $lower}else{$higher_val = $lower; $lower_val = $upper};

			my $old_part = $lower_val; my $new_part = $lower_val+$notch_size;
			while($new_part <= $higher_val)
				{
				# Y lower for curent segment = $old_part; Y upper is $new_part	
				my $SEG_XS = 1 * cos($old_part);my $SEG_YS = 1 * sin($old_part);
				my $SEG_XE = 1 * cos($new_part);my $SEG_YE = 1 * sin($new_part);
				#
				my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"$big_clade_color\", lwd = $big_clade_size)\n";
				$draw_tree_R_commands .= $R_command;
				$old_part = $new_part; $new_part += $notch_size;


				if($new_part > $higher_val)
					{
					# remaining segment before break loop
					my $SEG_XS = 1 * cos($old_part);my $SEG_YS = 1 * sin($old_part);
					my $SEG_XE = 1 * cos($higher_val);my $SEG_YE = 1 * sin($higher_val);
					my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"$big_clade_color\", lwd = $big_clade_size)\n";
					$draw_tree_R_commands .= $R_command;
					};
				}

			
			
		#	};





	};#if($branch_plot_type == 1)# single line branches:



if($branch_plot_type == 2) # nice lines ... 
	{
	my $current_loop_draw_tree_R_commands = "";


	#	if(($y1_proportion+$notch_size) >= $y2_proportion )
	#		{
		#	# just draw one line
		#	my $SEG_XS = $sum_branchlength * cos($y1_proportion);my $SEG_YS = $sum_branchlength * sin($y1_proportion);
		#	my $SEG_XE = $sum_branchlength * cos($y2_proportion);my $SEG_YE = $sum_branchlength * sin($y2_proportion);
		#	my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"black\", lwd = $width3)\n";
		#		$draw_tree_R_commands .= $R_command;
	
		#	}else{

			if($y1_proportion < $y2_proportion)
			{
			my $old_part = $y1_proportion; my $new_part = $y1_proportion+$notch_size;
			while($new_part <= $y2_proportion)
				{
				# Y lower for curent segment = $old_part; Y upper is $new_part	
				my $SEG_XS = $sum_branchlength * cos($old_part);my $SEG_YS = $sum_branchlength * sin($old_part);
				my $SEG_XE = $sum_branchlength * cos($new_part);my $SEG_YE = $sum_branchlength * sin($new_part);
#
				my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"$node_colour\", lwd = $width3)\n";

				# $draw_tree_R_commands .= $R_command;
				$current_loop_draw_tree_R_commands .= $R_command;

				$old_part = $new_part; $new_part += $notch_size;
				if($new_part > $y2_proportion)
					{
					# remaining segment before break loop
					my $SEG_XS = $sum_branchlength * cos($old_part);my $SEG_YS = $sum_branchlength * sin($old_part);
					my $SEG_XE = $sum_branchlength * cos($y2_proportion);my $SEG_YE = $sum_branchlength * sin($y2_proportion);
					my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"$node_colour\", lwd = $width3)\n";
					# $draw_tree_R_commands .= $R_command;
					$current_loop_draw_tree_R_commands .= $R_command;

					};
				}
			}else{
			my $old_part = $y2_proportion; my $new_part = $y2_proportion+$notch_size;
			while($new_part <= $y1_proportion)
				{
				# Y lower for curent segment = $old_part; Y upper is $new_part	
				my $SEG_XS = $sum_branchlength * cos($old_part);my $SEG_YS = $sum_branchlength * sin($old_part);
				my $SEG_XE = $sum_branchlength * cos($new_part);my $SEG_YE = $sum_branchlength * sin($new_part);
#
				my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"$node_colour\", lwd = $width3)\n";
				# $draw_tree_R_commands .= $R_command;
				$current_loop_draw_tree_R_commands .= $R_command;

				$old_part = $new_part; $new_part += $notch_size;
				if($new_part > $y1_proportion)
					{
					# remaining segment before break loop
					my $SEG_XS = $sum_branchlength * cos($old_part);my $SEG_YS = $sum_branchlength * sin($old_part);
					my $SEG_XE = $sum_branchlength * cos($y1_proportion);my $SEG_YE = $sum_branchlength * sin($y1_proportion);
					my $R_command =  "segments(" . "$SEG_XS, $SEG_YS , $SEG_XE, $SEG_YE , col = \"$node_colour\", lwd = $width3)\n";
					# $draw_tree_R_commands .= $R_command;
					$current_loop_draw_tree_R_commands .= $R_command;

					};
				}
			
			}

		#	};


	# horizontal line of tree, adjusted for circular:
	# start position.
	my $Xstart 	= $sum_branchlength * cos($y2_proportion);# X ax
	my $Ystart 	= $sum_branchlength * sin($y2_proportion);# Y ax
	# end position.
	my $Xend 	= $assign_new_x * cos($y2_proportion);# X
	my $Yend 	= $assign_new_x * sin($y2_proportion);# Y
	my $R_command =  "segments(" . "$Xstart, $Ystart , $Xend, $Yend , col = \"$node_colour\", lwd = $width3)\n";
	# $draw_tree_R_commands .= $R_command;
	$current_loop_draw_tree_R_commands .= $R_command;

	my $assign_new_x_offset = $assign_new_x + $plot_tip_label_offset;
	my $Xend_offset = $assign_new_x_offset * cos($y2_proportion);# labels are not right at tip, but a bit further out.

	my $Yend_offset	= ($assign_new_x + $plot_tip_label_offset) * sin($y2_proportion);# Y




	# the following bit, hack some visual representation of values which are on terminal labels
	# in this case, i have the temperature appended. 

	unless($test =~ /^INTERNAL_NODE_/)
		{
	#	print "tip label :$test\n"; # Cachryphora_canadensis_temp5.8

	my $Xstart2 	= $reference_tip_label_x * cos($y2_proportion);# X ax
	my $Ystart2 	= $reference_tip_label_x * sin($y2_proportion);# Y ax

	 # wasted 3 hours of my life finding out how to write this one line:
	my $textangle = atan2 ( $Ystart2 , $Xstart2 ) * (180 / 3.142) ;

	my $pinttextangle = $textangle;
	if($textangle > 90){$pinttextangle = $textangle - 180};
	if($textangle < -90){$pinttextangle = $textangle + 180};
	
	my $R_commandtip;



	if($all_tip_labels_eqidistant == 1)
		{
		 $R_commandtip = "par(srt=0)\ntext($Xstart2, $Ystart2,srt=$pinttextangle , col = \"black\", " . "labels=\"$test\",cex=$tip_tax_label_cex, font=1)\n";# REG
		}else{
		 $R_commandtip = "par(srt=0)\ntext($Xend_offset, $Yend_offset,srt=$pinttextangle , col = \"black\", " . "labels=\"$test\",cex=$tip_tax_label_cex, font=1)\n";# REG
		};

	if($print_tips_circular_tree == 1)
		{$draw_tree_R_commands .= $R_commandtip};


		if($test =~ /_temp(\-*\d+\.\d+)/ && $plot_climate_tree == 1)
			{
			my $temperauter = $1;$temperaure_point_cex = 2;
			my $point_color;

			# simpler:
		#	if($temperauter >= 10){$point_color = "red"}else{$point_color = "blue"};
		#	my $R_command = "points($Xend, $Yend,col=\"$point_color\", pch=16, cex=$temperaure_point_cex)\n";


#	temps range from -32 to +32, so to transform them to be between 0 and 1:
#	(x - -32) / 64
#	colors2<- rgb(colorfunc(0.1) , maxColorValue=255)
#	points(4,2,pch=19,cex=12,col=colors2)

			if($temperauter <= 40)
				{
			my $transformed = ( $temperauter - -32 ) / 64; 
			my $R_command = "colors2<- rgb(colorfunc($transformed) , maxColorValue=255)\n" . 					
			 	"points($Xend, $Yend,col=colors2, pch=16, cex=$temperaure_point_cex)\n";
			$draw_tree_R_commands4 .= $R_command;
				};

			};



			# if want to plot them exactly on tip (good if have branch lengths), then use: $Xend, $Yend
			# otherwise:

			$tip_highlight_X = $highlight_tips_X * cos($y2_proportion);
			$tip_highlight_Y = $highlight_tips_X * sin($y2_proportion);


		if($highlight_tip_string_match == 1 && $test =~ /$tip_string_to_highlight/)
			{
			my $R_command = "points($Xend, $Yend,col=\"red\", pch=16, cex=$highlight_tip_string_cex)\n";
			$draw_tree_R_commands4 .= $R_command;
			};

		if($highlight_tips_from_file == 1 && $highlight_tips{$test} =~ /\w/)
			{
			my $current_tip_color = $highlight_tips{$test};

			if($highlight_tips_type == 1)
			{

			############################################
			print "tip from file found, $test, will highlight\n";
			my $R_command = "points($Xend, $Yend,col=\"$current_tip_color\", pch=16, cex=$highlight_tip_string_cex)\n";
		#	my $R_command = "points($tip_highlight_X, $tip_highlight_Y,col=\"red\", pch=16, cex=$highlight_tip_string_cex)\n";
			$draw_tree_R_commands4 .= $R_command;
			############################################

			}else{
 
			############################################
			my $R_command;
			my $XstartBBB;my $YstartBBB;my $XendBBB;my $YendBBB;

			if($highlight_tips_eqidistant == 1)
				{
				# start position.
				 $XstartBBB 	= $rainbow_X1 * cos($y2_proportion);# X ax
				 $YstartBBB 	= $rainbow_X1 * sin($y2_proportion);# Y ax
				# end position.
				 $XendBBB 	= $rainbow_X2 * cos($y2_proportion);# X
				 $YendBBB 	= $rainbow_X2 * sin($y2_proportion);# Y
				}else{
				 $XstartBBB 	= $assign_new_x * cos($y2_proportion); $YstartBBB 	= $assign_new_x * sin($y2_proportion);# Y ax
				 $XendBBB 	= ($assign_new_x+$colored_tip_bar_length) * cos($y2_proportion); $YendBBB 	= ($assign_new_x+$colored_tip_bar_length) * sin($y2_proportion);# Y
				};

			my $R_command =  "segments(" . "$XstartBBB, $YstartBBB,$XendBBB , $YendBBB , " . "col = \"$current_tip_color\", lwd = $rainbow_line_width)\n";
			$draw_tree_R_commands4  .= $R_command;
			############################################

			}



			};


		my $plot_pie_at_tip = 0; if($known_state{$test} =~ /\w/ || $predicted_state{$test} =~ /\w/){$plot_pie_at_tip = 1}; 
		if($piechart_trait_prob == 1 && $plot_pie_at_tip == 1)
			{# 2021jan, pie charts at tips to depict trait probabilities


			print "ploting pie chart for tip $test\n";
			# $known_state{$sp}; $predicted_state{$sp} = "A"; $state_probabilities_for_pie{$sp}{$pr} = $p;

			# only one state, just draw a circle. also with black boarder to signify prior known,
			if($known_state{$test} =~ /\w/)
				{
				my $trait_state = $known_state{$test};my $plotcolor = $trait_state_colors{$trait_state};
				unless($plotcolor =~ /\w/){die "\nfailed to retrive color for state $trait_state\n"};

				my $R_command = "draw.circle($Xend, $Yend,$pie_radius_border,border=NULL,col=\"black\",lty=1, lwd=1); " . 
					"draw.circle($Xend, $Yend,$pie_radius,border=NULL,col=\"$plotcolor\",lty=1, lwd=1)\n"; # nv=100,
				$draw_tree_R_commands4 .= $R_command;
				};
			
			# multiple states, draw pie chart
			if($predicted_state{$test} =~ /\w/)
				{
				my @color_array;my @prob_array;
				for my $state(@all_states_for_pie_array)
					{
					my $color = $trait_state_colors{$state};
					my $prob = 0; if($state_probabilities_for_pie{$test}{$state} =~ /\d/){$prob = $state_probabilities_for_pie{$test}{$state}};
					push @color_array, $color; push @prob_array, $prob;
					};
				my $colorstring = join '","' , @color_array;
				my $statesstring = join '","' , @all_states_for_pie_array;
				my $probsstring = join ',' , @prob_array;
				my $R_command = "draw.circle($Xend, $Yend,$pie_radius_border,border=NULL,col=\"gray\",lty=1, lwd=0.1); " .
					"states <- c(\"" . $statesstring . "\");probs <- c("  . $probsstring . ");trait_cols <- c(\"" . $colorstring . "\"); " .
 					"t9 <- matrix (NA, nrow=length(states), ncol = 4);t9[ , 1] <- $Xend;t9[ , 2] <- $Yend;t9[ , 3] <- states;t9[ , 4] <- probs;" . 
 					" xyz <- make.xyz(as.numeric(t9[ , 1]),as.numeric(t9[ , 2]),as.numeric(t9[,4]),t9[,3]);" . 
 					" draw.pie(xyz\$x, xyz\$y, xyz\$z, radius = $pie_radius, col=trait_cols) # $test\n";

				$draw_tree_R_commands4b .= $R_command;

 # states <- c("A", "B", "C"); probs <- c(0.1, 0.1, 0.8);trait_cols <- c("red","brown","yellow" )
 # t9 <- matrix (NA, nrow=length(states), ncol = 4)
 # t9[ , 1] <- 0.5;t9[ , 2] <- 0.75;t9[ , 3] <- states;t9[ , 4] <- probs;
 # xyz <- make.xyz(as.numeric(t9[ , 1]),as.numeric(t9[ , 2]),as.numeric(t9[,4]),t9[,3])
 # draw.pie(xyz$x, xyz$y, xyz$z, radius = 0.5, col=trait_cols)


				};
			
			};




		};


	# $draw_tree_R_commands .= $R_command;
	# $default_node_color = "grey"; $node_colour

	if($node_colour eq $default_node_color)
		{
		$draw_tree_R_commands .= $current_loop_draw_tree_R_commands;
		}else{
		$draw_tree_R_commands6 .= $current_loop_draw_tree_R_commands;

		if($label_painted_nodes == 1)
			{
			
			my $R_command = "par(srt=0)\ntext($Xstart, $Ystart, col = \"black\", " . 
			"labels=\"$abreviate_nodeID$internal_taxname\",cex=$painted_branch_label_cex, font=2)\n";
			$draw_tree_R_commands7 .= $R_command;

			};


		};

	};#if($branch_plot_type == 2)




# here highlight and print tax of internal nodes defined by list of decendents.
if($label_these_decendent_defined_nodes_found{$current_node} == 1)
	{
	my $R_command = "points($Xstart, $Ystart, col = \"blue\", pch=16, cex=10)\n" .
		"par(srt=0)\ntext($Xstart, $Ystart, col = \"black\", " . 
		"labels=\"$abreviate_nodeID$internal_taxname\",cex=5, font=2)\n";
	$draw_tree_R_commands3 .= $R_command;
	};


if($count_terminals_node_lable{$current_node} >= $clade_size_limit_for_interal_label && 
	$print_internal_labels_circular_tree == 1)
	{
	my $R_command = "par(srt=0)\ntext($Xstart, $Ystart, col = \"black\", " . 
	"labels=\"$abreviate_nodeID$internal_taxname\",cex=$ref_tree_internal_label_cexB, font=2)\n";# REG
	$draw_tree_R_commands3 .= $R_command;
	};

if( $internal_names_circular_tree{"$abreviate_nodeID"} =~ /\w/)
	{
	my $nodestring49 = $internal_names_circular_tree{"$abreviate_nodeID"};


	my $R_command = "par(srt=0)\ntext($Xstart, $Ystart, col = \"black\", " . 
	"labels=\"$nodestring49\",cex=$ref_tree_internal_label_cexB, font=2)\n";# REG
#	$draw_tree_R_commands3 .= $R_command;

	print_shadow( $Xstart, $Ystart, $nodestring49 ,$rectangular_tree_internal_label_cex);
		
	
	};







if($print_internal_node_labels_on_rectangular_tree == 1)
	{

	print_shadow( $sum_branchlength, $y1_proportion, "$abreviate_nodeID$internal_taxname" ,$rectangular_tree_internal_label_cex);

	my $R_command = "par(srt=0)\ntext($sum_branchlength, $y1_proportion, col = \"black\", " . 
	"labels=\"$abreviate_nodeID$internal_taxname\",cex=$rectangular_tree_internal_label_cex, font=2)\n";# REG
	$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;	
	};



	# stuff inserted for tree drawing

	#######################################################################################
	#
	#######################################################################################







};#sub get_tree_plotting_commands


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub label_user_defined_tips
{
my $tree = shift;

#print "\ntree:$tree\n";# (Squamipalpis_pantoea_HQ921524:0.06063191365842858505,Stenhyp

#	$tip_label_string = $1;$tip_label_file = $2;
#	$add_string_to_this_member{$1} = 1;

$tip_label_string = $tip_label_string . "_";

my @keys = keys %add_string_to_this_member;@keys = sort @keys;

foreach my $label_this(@keys)
	{
#	print "label_this:$label_this\n";
	if($tree =~ s/$label_this\:/$tip_label_string$label_this:/)
		{
		$sucessfully_labelled++;
		}else{
		$not_sucessfully_labelled++;
		
		};	

	};


print "
looked for $#keys user specified tips:
	sucessfully_labelled:$sucessfully_labelled
	not_sucessfully_labelled:$not_sucessfully_labelled
";



return($tree);


};



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub return_tax_string_for_this
{
my $look_for_ID = shift;
my $add_tax_string = ""; # "NA_";

my $completelineage;my $genus;
if(exists($complete_lineage_for_this_species{$look_for_ID}))
	{
	$completelineage = $complete_lineage_for_this_species{$look_for_ID};
	}else{
	$genus = $look_for_ID;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
#	print "no sucess looking for binomial:$binomial, trying genus name:$genus\n";
	$completelineage = $complete_lineage_for_this_species{$genus};
	}

if($completelineage =~ /\w/)
	{
	$add_tax_string = "";my $tax_added = 0;

	foreach my $current_rank (@ranksarray)
		{
		if($completelineage =~ /($current_rank):(\S+)/i)
			{
			my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
			$rank_assignment_counts{$rank}++;

			if($truncate_new_taxonomic_strings == 1)
				{
				$tax =~ s/^(\w{10})\w+/$1/;
				};
			if($uppercase_tax_strings == 1)
				{$tax  =uc($tax)};

			$add_tax_string .= "$tax" . "_";$tax_added++;
			};
		};

		unless($tax_added >= 1)
			{$add_tax_string = ""; # "NA_"; 
			if($warning_print_limit_2 < 250)
				{
			print "warning 3091, maybe inapproriate taxonomic ranks selected, because non found in current lineage\n";
			print "\tcompletelineage:$completelineage\n";
				}elsif($warning_print_limit_2 == 250)
				{
				print " ....... warning print limit reached ..........\n";
				};
			$warning_print_limit_2++;
			};


		unless($add_tax_string =~ /\w/)
			{
		#	foreach my $current_rank (@ranksarray)
		#		{print "\tcurrent_rank:$current_rank\n"};
		#	die "\n\nwhat happened? completelineage:$completelineage\n";
			};

		}else{
		$newick_binomials_for_which_taxonomies_not_found++;	
		$print_warnings++;
		if($print_warnings>=100)
			{
			print "no lineage found for binomial:$binomial or genus:$genus\n";
			$print_warnings = 0;
			};
		
		};


return($add_tax_string);


};


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub print_shadow
{
my $XstartC = $_[0]; my $YstartC = $_[1]; 
my $printlable = $_[2]; my $rectangular_tree_internal_label_cex = $_[3];

my $text_shadow_col = "white"; # azure4
my $pinttextangle = 0;
my $ref_tree_internal_label_cex = $rectangular_tree_internal_label_cex;
my $text_shadow_offset = 0.005;

			# lower left
			my $R_command = "par(srt=0)\ntext($XstartC-$text_shadow_offset, $YstartC-$text_shadow_offset,srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# lower mid
			my $R_command = "par(srt=0)\ntext($XstartC, $YstartC-$text_shadow_offset,srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# lower right
			my $R_command = "par(srt=0)\ntext($XstartC+$text_shadow_offset, $YstartC-$text_shadow_offset,srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# upper left
			my $R_command = "par(srt=0)\ntext($XstartC-$text_shadow_offset, $YstartC+$text_shadow_offset,srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# upper mid
			my $R_command = "par(srt=0)\ntext($XstartC, $YstartC+$text_shadow_offset,srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# upper right
			my $R_command = "par(srt=0)\ntext($XstartC+$text_shadow_offset, $YstartC+$text_shadow_offset,srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;


			# mid left
			my $R_command = "par(srt=0)\ntext($XstartC-$text_shadow_offset, $YstartC+($text_shadow_offset*0.5),srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# mid right
			my $R_command = "par(srt=0)\ntext($XstartC+$text_shadow_offset, $YstartC+($text_shadow_offset*0.5),srt=$pinttextangle , col = \"$text_shadow_col\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;

			# main
			my $R_command = "par(srt=0)\ntext($XstartC, $YstartC,srt=$pinttextangle , col = \"black\", " . 
				"labels=\"$printlable\",cex=$ref_tree_internal_label_cex, font=2)\n";# REG
			$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;$draw_tree_R_commands3 .= $R_command;



#	my $R_command = "par(srt=0)\ntext($sum_branchlength, $y1_proportion, col = \"black\", " . 
#	"labels=\"$abreviate_nodeID$internal_taxname\",cex=$rectangular_tree_internal_label_cex, font=2)\n";# REG
#	$draw_tree_R_commands_RECTANGL_internal_node_labels .= $R_command;	



}

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub plot_ITOL2_9_trees
{
# WP:
$ylim 				= "c(3,9)";# regular :c(-16,16)
$notchsize = 0.081;
$draw_clades_WP = "
rect(1.04, 7.8649, 1.6, 9.162, col=\"gray\",  lwd=0.001, border=\"gray\",lty=NULL,  xpd=FALSE)
rect(1.04,7.378 , 1.6, 7.8649, col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)
rect(1.04,5.270 , 1.6, 7.378, col=\"yellow\",  lwd=0.001, border=\"yellow\",lty=NULL,  xpd=FALSE)
rect(1.04,3.162 , 1.6,5.270 , col=\"pink\",  lwd=0.001, border=\"pink\",lty=NULL,  xpd=FALSE)
text(1.5,8.5,labels=\"Archaeognatha\",cex=1.2, font=2,srt=-45, col = \"white\")
text(1.5,7.6,labels=\"Zygentoma\",cex=1.2, font=2,srt=-45, col = \"white\")
text(1.5,6.3,labels=\"Odonata\",cex=1.2, font=2,srt=-45, col = \"white\")
text(1.5,4.2,labels=\"Ephemeroptera\",cex=1.2, font=2,srt=-45, col = \"white\")





";
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_WP;


#############################################################################################


# POLY
$ylim 				= "c(5.2,6.8)";# regular :c(-16,16)

$collparse_nodes = 0;
if($collparse_nodes == 1)
{
$collapse_nodes{"INTERNAL_NODE_175"} = "Dictyoptera";
$collapse_nodes{"INTERNAL_NODE_81"} = "Acrididae";
$collapse_nodes{"INTERNAL_NODE_28"} = "Tettigoniidae";
$collapse_nodes{"INTERNAL_NODE_45"} = "Pamphagidae";
$collapse_nodes{"INTERNAL_NODE_98"} = "Phasmatidae";
$collapse_nodes{"INTERNAL_NODE_4"} = "Perlidae";
$collapse_nodes{"INTERNAL_NODE_93"} = "Diapheromeridae";
$collapse_nodes{"INTERNAL_NODE_10"} = "Grylloidea";
$collapse_nodes{"INTERNAL_NODE_88"} = "Acridomorpha";
};

$draw_clades_Pol = "
rect(1.04,6.34920634920635 , 1.6,6.50793650793651  , col=\"yellow\",  lwd=0.001, border=\"yellow\",lty=NULL,  xpd=FALSE)
text(1.5, 6.42 ,labels=\"Plecoptera\",cex=0.8, font=2,srt=-33, col = \"white\")
rect(1.04, 6 , 1.6, 6.15873015873016 , col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
text(1.5, 6.065 ,labels=\"Verophasmatodea\",cex=0.8, font=2,srt=-33, col = \"white\")
rect(1.04, 5.80952380952381 , 1.6, 5.96825396825397 , col=\"lightgreen\",  lwd=0.001, border=\"lightgreen\",lty=NULL,  xpd=FALSE)
text(1.5, 5.89 ,labels=\"Caelifera\",cex=0.8, font=2,srt=-33, col = \"white\")
rect(1.04, 5.52380952380952 , 1.6, 5.77777 , col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)
text(1.5, 5.64 , labels=\"Ensifera\",cex=0.8, font=2,srt=-33, col = \"white\")
";

#############################################################################################

# Dicty
$ylim 				= "c(5,7)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.8;	# WP=0.6

# $collapse_nodes{"INTERNAL_NODE_56"} = "Termitidae";
# $collapse_nodes{"INTERNAL_NODE_16"} = "Rhinotermitinae";
# $collapse_nodes{"INTERNAL_NODE_13"} = "Kalotermitidae";
# $collapse_nodes{"INTERNAL_NODE_6"} = "Blaberoidea";

$draw_clades_Dicty = "
rect(1.04, 6.37, 1.6, 6.98630136986301, col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
rect(1.04, 6.49315068493151 , 1.4, 6.65753424657534, col=\"yellow\",  lwd=0.001, border=\"yellow\",lty=NULL,  xpd=FALSE)
rect(1.04, 6.37 , 1.6, 5.17808219178082 , col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)

text(1.35, 6.6 ,labels=\"Mantodea\",cex=1.0, font=2,srt=-45, col = \"white\")
text(1.5, 5.7,labels=\"Isoptera\",cex=1.2, font=2,srt=-45, col = \"white\")
text(1.5, 6.7,labels=\"Blattodea\",cex=1.2, font=2,srt=-45, col = \"white\")

";

#############################################################################################

# PARA

$ylim 				= "c(5,7)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.6;	# WP=0.6

$collparse_nodes_Para = 0;
if($collparse_nodes_Para == 1)
{

$collapse_nodes{"INTERNAL_NODE_29"} = "Aphididae";
$collapse_nodes{"INTERNAL_NODE_117"} = "Pentatomoidea";
$collapse_nodes{"INTERNAL_NODE_90"} = "Miridae";
$collapse_nodes{"INTERNAL_NODE_80"} = "Reduviidae";
$collapse_nodes{"INTERNAL_NODE_39"} = "Aleyrodidae";
$collapse_nodes{"INTERNAL_NODE_59"} = "Cercopoidea";
$collapse_nodes{"INTERNAL_NODE_51"} = "Membracoidea";
$collapse_nodes{"INTERNAL_NODE_33"} = "Psylloidea";
$collapse_nodes{"INTERNAL_NODE_93"} = "Cimicoidea";
$collapse_nodes{"INTERNAL_NODE_20"} = "Fulgoroidea";
$collapse_nodes{"INTERNAL_NODE_6"} = "Philopteridae";
$collapse_nodes{"INTERNAL_NODE_54"} = "Cicadidae";
$collapse_nodes{"INTERNAL_NODE_102"} = "Lygaeoidea";
$collapse_nodes{"INTERNAL_NODE_1"} = "Thripinae";
$collapse_nodes{"INTERNAL_NODE_13"} = "Peloridiidae";
$collapse_nodes{"INTERNAL_NODE_67"} = "Nepinae";
$collapse_nodes{"INTERNAL_NODE_67"} = "Belostomatidae";
$collapse_nodes{"INTERNAL_NODE_107"} = "Coreoidae";
$collapse_nodes{"INTERNAL_NODE_69"} = "Belostomatidae";
};
 
$draw_clades_Para = "
rect(1.04,5.20408163 , 1.6, 6.59183673469388 , col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
rect(1.04, 5.2040816 , 1.4, 5.53061224, col=\"lightgreen\",  lwd=0.001, border=\"lightgreen\",lty=NULL,  xpd=FALSE)
rect(1.04, 5.6530612 , 1.4, 5.97959183 , col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)

text(1.5, 6,labels=\"Hemiptera\",cex=1.8, font=2,srt=-65, col = \"white\")
text(1.30, 5.35 ,labels=\"Pentatomomorpha\",cex=0.7, font=2,srt=-65, col = \"white\")
text(1.30, 5.8 ,labels=\"Nepomorpha\",cex=0.7, font=2,srt=-65, col = \"white\")
";




#############################################################################################

# Hym
$ylim 				= "c(5,7)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.6;	# WP=0.6

$collparse_nodes_Hym = 0;
if($collparse_nodes_Hym == 1)
{
$collapse_nodes{"INTERNAL_NODE_80"} = "Ichneumonoidea";
$collapse_nodes{"INTERNAL_NODE_46"} = "Apoidea";
$collapse_nodes{"INTERNAL_NODE_35"} = "Formicidae";
$collapse_nodes{"INTERNAL_NODE_26"} = "Vespidae";
$collapse_nodes{"INTERNAL_NODE_23"} = "Tenthredinoidea";
$collapse_nodes{"INTERNAL_NODE_9"} = "Scelionidae";
}


#############################################################################################


# Neur
$ylim 				= "c(3,9)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.8;	# WP=0.6

$collparse_nodes_Neur = 0;
if($collparse_nodes_Neur == 1)
{
$collapse_nodes{"INTERNAL_NODE_4"} = "Corydalinae";
$collapse_nodes{"INTERNAL_NODE_0"} = "Chauliodinae";
$collapse_nodes{"INTERNAL_NODE_14"} = "Chrysopidae";
$collapse_nodes{"INTERNAL_NODE_8"} = "Myrmeleontidae";
};

$draw_clades_Neur = "
 rect(1.04,4.666666 , 1.6,6.888888 , col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
 text(1.525, 5.6, labels=\"Neuroptera\",cex=1.6, font=2,srt=-90, col = \"white\")
";



#############################################################################################

# CS
$ylim 				= "c(5.5,6.3)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.9;	# WP=0.6

$collparse_nodes_CS = 0;
if($collparse_nodes_CS == 1)
{
$collapse_nodes{"INTERNAL_NODE_260"} = "Chrysomelidae";
$collapse_nodes{"INTERNAL_NODE_243"} = "Cerambycidae";
$collapse_nodes{"INTERNAL_NODE_231"} = "Curculionidae";
$collapse_nodes{"INTERNAL_NODE_177"} = "Tenebrionidae";
$collapse_nodes{"INTERNAL_NODE_150"} = "Staphylinidae";
$collapse_nodes{"INTERNAL_NODE_93"} = "Leiodidae";
$collapse_nodes{"INTERNAL_NODE_88"} = "Scarabaeidae";
$collapse_nodes{"INTERNAL_NODE_47"} = "Elateridae";
$collapse_nodes{"INTERNAL_NODE_164"} = "Coccinellidae";
$collapse_nodes{"INTERNAL_NODE_15"} = "Caraboidea";
$collapse_nodes{"INTERNAL_NODE_6"} = "Dytiscoidea";
$collapse_nodes{"INTERNAL_NODE_51"} = "Elateroidea";
$collapse_nodes{"INTERNAL_NODE_64"} = "Hydrophilidae";
$collapse_nodes{"INTERNAL_NODE_189"} = "Tenebrionoidea";
$collapse_nodes{"INTERNAL_NODE_58"} = "Bostrichoidea";
$collapse_nodes{"INTERNAL_NODE_60"} = "Histeridae";
$collapse_nodes{"INTERNAL_NODE_74"} = "Scarabaeoidea";
$collapse_nodes{"INTERNAL_NODE_156"} = "Cleroidea";
$collapse_nodes{"INTERNAL_NODE_28"} = "Buprestidae";
$collapse_nodes{"INTERNAL_NODE_24"} = "Byrrhoidea";
$collapse_nodes{"INTERNAL_NODE_268"} = "Cucujiformia";

# half notch size= 0.00995
$draw_clades_CS = "
 rect(1.04,6.22920265780731 , 1.85, 6.28906976744186 , col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)
 rect(1.04,5.76079734219269 , 1.85, 6.22920265780731 , col=\"lightpink\",  lwd=0.001, border=\"lightpink\",lty=NULL,  xpd=FALSE)
 rect(1.04,6.03 , 1.53, 6.13953488372093, col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
 rect(1.04,5.76079734219269 , 1.53, 6.03, col=\"lightgreen\",  lwd=0.001, border=\"lightgreen\",lty=NULL,  xpd=FALSE)
 text(1.575,6.25,labels=\"Strepsiptera\",cex=0.8, font=2,srt=-45, col = \"white\")
 text(1.575,6.00,labels=\"Coleoptera\",cex=1.4, font=2,srt=-45, col = \"white\")
 text(1.40, 5.9 ,labels=\"Polyphaga\",cex=1.0, font=2,srt=-45, col = \"white\")
 text(1.40, 6.09 ,labels=\"Adephaga\",cex=1.0, font=2,srt=-45, col = \"white\")
";

};


#############################################################################################

# MSD
$ylim 				= "c(4.25,7.4)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.4;	# WP=0.6

$collparse_nodes_MSD = 0;
if($collparse_nodes_MSD == 1)
{
$collapse_nodes{"INTERNAL_NODE_73"} = "Calliphoridae";
$collapse_nodes{"INTERNAL_NODE_59"} = "Muscoidea";
$collapse_nodes{"INTERNAL_NODE_42"} = "Acalyptratae";
$collapse_nodes{"INTERNAL_NODE_30"} = "Tabanoidea";
$collapse_nodes{"INTERNAL_NODE_18"} = "Culicidae";
$collapse_nodes{"INTERNAL_NODE_63"} = "Sarcophagidae";
$collapse_nodes{"INTERNAL_NODE_25"} = "Sciaroidea";
$collapse_nodes{"INTERNAL_NODE_48"} = "Tachinidae";
$collapse_nodes{"INTERNAL_NODE_44"} = "Tephritidae";
$collapse_nodes{"INTERNAL_NODE_14"} = "Chironomoidea";
$collapse_nodes{"INTERNAL_NODE_34"} = "Emphidoidea";
$collapse_nodes{"INTERNAL_NODE_7"} = "Trichoceridae";
$collapse_nodes{"INTERNAL_NODE_6"} = "Ptychopteridae";
$collapse_nodes{"INTERNAL_NODE_3"} = "Phlebotominae";

$draw_clades_MSD = "

 rect(1.04, 4.95412844036697 , 1.6, 6.93577981651376 , col=\"lightpink\",  lwd=0.001, border=\"lightpink\",lty=NULL,  xpd=FALSE)
 rect(1.04, 4.95412844036697 , 1.4, 6.22018348623853 , col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
 rect(1.04, 6.6605504587156, 1.4,  6.93577981651376 , col=\"lightgreen\",  lwd=0.001, border=\"lightgreen\",lty=NULL,  xpd=FALSE)
 rect(1.04, 6.38532110091743 , 1.4, 6.60550458715596 , col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)

 text(1.55, 6 ,labels=\"Diptera\",cex=1.8, font=2,srt=-90, col = \"white\")
 text(1.35, 5.5 ,labels=\"Brachycera\",cex=0.7, font=2,srt=-40, col = \"white\")
 text(1.35, 6.8 ,labels=\"Nematocera\",cex=0.7, font=2,srt=-40, col = \"white\")
 text(1.35, 6.5 ,labels=\"Culicomorpha\",cex=0.7, font=2,srt=-40, col = \"white\")
   ";
   
   




};


#############################################################################################


# TL

$ylim 				= "c(5.5,6.5)";# regular :c(-16,16)
$rectangular_plot_tip_label_cex = 0.6;	# WP=0.6

$collparse_nodes_TL = 1;
if($collparse_nodes_TL == 1)
{
 $collapse_nodes{"INTERNAL_NODE_210"} = "Nymphalidae";
 $collapse_nodes{"INTERNAL_NODE_123"} = "Noctuoidea";
 $collapse_nodes{"INTERNAL_NODE_95"} = "Bombycoidea";
 $collapse_nodes{"INTERNAL_NODE_19"} = "Tortricidae";
 $collapse_nodes{"INTERNAL_NODE_50"} = "Hesperidae";
 $collapse_nodes{"INTERNAL_NODE_29"} = "Gelechioidea";
 $collapse_nodes{"INTERNAL_NODE_213"} = "Papilionoidea";
 $collapse_nodes{"INTERNAL_NODE_84"} = "Pyraloidea";
 $collapse_nodes{"INTERNAL_NODE_6"} = "Hepialidae";
 $collapse_nodes{"INTERNAL_NODE_59"} = "Geometridae";
 $collapse_nodes{"INTERNAL_NODE_61"} = "Thyrididae";
 $collapse_nodes{"INTERNAL_NODE_7"} = "Lithocolletinae";
 $collapse_nodes{"INTERNAL_NODE_8"} = "Zygaenidae";


# rect(1.04, , 1.6, , col=\"\",  lwd=0.001, border=\"\",lty=NULL,  xpd=FALSE)
# text(1.5,,labels=\"\",cex=1.2, font=2,srt=-45, col = \"white\")


$draw_clades_TL = "

rect(1.04, 6.345 , 1.6, 6.43388429752066 , col=\"lightgray\",  lwd=0.001, border=\"lightgray\",lty=NULL,  xpd=FALSE)
rect(1.04, 5.615702 , 1.6, 6.345 , col=\"lightpink\",  lwd=0.001, border=\"lightpink\",lty=NULL,  xpd=FALSE)
rect(1.04,5.61570247933884 , 1.5, 6.260330578, col=\"lightblue\",  lwd=0.001, border=\"lightblue\",lty=NULL,  xpd=FALSE)
rect(1.04, 5.615702 , 1.4, 5.838842 , col=\"lightgreen\",  lwd=0.001, border=\"lightgreen\",lty=NULL,  xpd=FALSE)

text(1.54, 6.375 ,labels=\"Trichoptera\",cex=1.0, font=2,srt=-45, col = \"white\")
text(1.54, 6.0 ,labels=\"Lepidoptera\",cex=1.0, font=2,srt=-45, col = \"white\")
text(1.4, 6 ,labels=\"Ditrysia\",cex=1.0, font=2,srt=-45, col = \"white\")
text(1.3, 5.7 ,labels=\"Obtectomera\",cex=1.0, font=2,srt=-45, col = \"white\")


6.33471074380165

";



};

#############################################################################################



# $draw_tree_R_commands_RECTANGL3 = $draw_clades_WP;
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_Pol;
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_Dicty;
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_Poly;
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_Neur;
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_CS;
# $draw_tree_R_commands_RECTANGL3 = $draw_clades_TL;

# rect(1.04, , 1.6, , col=\"\",  lwd=0.001, border=\"\",lty=NULL,  xpd=FALSE)
# text(1.5,,labels=\"\",cex=1.2, font=2,srt=-45, col = \"white\")


};# sub plot_ITOL2_9_trees


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub remove_polytomies
{
my $treefile = shift;

print "
sub remove_polytomies
";


my $tree1= "";
open(IN, $treefile) || die "\n\nerror 2927 cant open file named:$treefile\n\n";
while (my $line = <IN>)
	{if($line =~ /(.+\(.+)/){$tree1 = $1}};
close(IN);
$tree1 =~ s/ //g;$tree_length = length($tree1);

print "\nnewick string of length $tree_length has been read from file:$treefile\n\n";


	my $count_scientific_notation_branchlengths_removed=0;
	my $count_regular_format_branchlengths_removed=0;
	my $count_zero_length_branchlengths_removed=0;


	$remove_branchlengths_while_removing_polytomies = 0;
	if($remove_branchlengths_while_removing_polytomies == 1)
		{
	print "\n3600. user chosen to remove branch-lengths in input tree ...\n";

	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	while($tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/)
		{$count_scientific_notation_branchlengths_removed++}; 
	# remove regular branchlengths: 0.02048
	while($tree1 =~ s/\:\-*\d+\.\d+//)
		{$count_regular_format_branchlengths_removed++}; 
	# remove 0 length branchlengths
	while($tree1 =~ s/\:\d+//)
		{$count_zero_length_branchlengths_removed++};

	print "\t$count_scientific_notation_branchlengths_removed count_scientific_notation_branchlengths_removed\n";
	print "\t$count_regular_format_branchlengths_removed count_regular_format_branchlengths_removed\n";
	print "\t$count_zero_length_branchlengths_removed count_zero_length_branchlengths_removed\n\n";

	if($tree1 =~ /..\:../){die "\nhuh\n"};
		}else{
		print "\nwarning ... not removing brnahclengths while removing polytomies, might not work!\n";
		};

	print "\n3620. user chosen to remove branch-support in input tree ...\n";
	while($tree1 =~ s/(\))\d\.\d+/$1/)
		{$count_proportional_node_support_removed++}; 

	print "
	count_proportional_node_support_removed:$count_proportional_node_support_removed
";


my $newick_string 	= $tree1;

open(LOG8, ">parse_newick_string_log");

#############################################################################################################

print "parsing newick string ....\n";

while ($newick_string =~ /\(([^\(\)]+)\)/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
	{
	my $node = $1; my $nodeID = "INTERNAL_NODE_$interal_node"; my $length_newick_string = length($newick_string);
	my @child_nodes = split /\,/ , $node;
	my $pinrtstring = $node;$pinrtstring =~ s/(.{100}).+/$1/;
print LOG8 "$nodeID\t$#child_nodes\t$pinrtstring\n";

	if(scalar @child_nodes == 1)# an error, internal node with no split
		{

		if($node =~ /^INTERNAL_NODE_\d+$|INTERNAL_NODE_\d+\:\d+\.\d+$/)
			{
			$newick_string =~ s/\(([^\(\)]+)\)/$1/;
			print "why only one child?\n";
			}else{
			print "some other .... length_newick_string:$length_newick_string ID:$nodeID\n";
			if($length_newick_string <= 500){print "\tnode:$node\n";};
			$interal_node++;
			$root_node = $nodeID;

			};
	#	if($newick_string =~ /(.{20}\([^\(\)]+\).{20})/)
	#		{
	#		my $tring = $1;print "\n$tring\n";
	#		};


		}
	elsif(scalar @child_nodes == 2)
		{
	#	print "bifurcation $#child_nodes\n";

		$child_counts{$nodeID} = $#child_nodes;

		for $i(0 .. $#child_nodes)
			{
			$nodes{$nodeID}{$i} = $child_nodes[$i];$nodes{$child_nodes[$i]}{parent} = $nodeID;
			unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++}
			};

		$newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/;

		$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
		$interal_node++;
		}else{
		print "POLYTOMIE, child nodes:$#child_nodes, IN:$interal_node\n";

		while($node =~ s/([^\,]+)\,([^\,]+)/INTERNAL_NODE_$interal_node/)
			{
			my $member1 = $1; my $member2 = $2;			
			print "\tmember1:$member1 member2:$member2, repl:INTERNAL_NODE_$interal_node\n";
			$child_counts{$nodeID} = 1;
				$nodes{$nodeID}{0} = $member1;$nodes{$member1}{parent} = $nodeID;
				$nodes{$nodeID}{1} = $member2;$nodes{$member2}{parent} = $nodeID;


			if($node =~ /([^\,]+)\,([^\,]+)/)
				{
				print "\t\tnothing\n";
				}else{
				$newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/;
				print "\t\tnewick string replacment with INTERNAL_NODE_$interal_node\n";
				};
			$root_node = $nodeID;
			$interal_node++;$nodeID = "INTERNAL_NODE_$interal_node";
			};

			
		};



# $newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/

	}#while ($newick_string =~


close LOG8;

#############################################################################################################


my @terminal_array = keys %terminal_labels;
$count_terminals = scalar @terminal_array;

print "
newick string stored, 
	count of polytomies:$polytomies
	count of internal nodes:$interal_node
	count terminals:$count_terminals
	branchlengths read:$bls_read
	newick_string:$newick_string
";




unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}

$new_newick_string99 = "($root_node)";



##################################################
traverse_backbone_tree_and_print_it($root_node);#
##################################################





open(OUTPR, ">$outprefix.poly_rm") || die "\nerror 2502\n";
print OUTPR "$new_newick_string99\n";
close OUTPR;





print "\npolytomeis rm (output is named $outprefix.poly_rm), quitting.\n";
exit;
};


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub traverse_backbone_tree_and_print_it
{
my $next = $_[0];
my $parent_node = $nodes{$next}{parent};#my $parent_node = $nodes{$next}{parent};

my $child_nodes = $child_counts{$next};#my $child_nodes = $child_counts{$next};
unless($child_nodes =~ /\d/){$child_nodes = 1};
my @next0 = ();


for $i(0 .. $child_nodes)
	{
	my $push_node = $nodes{$next}{$i};
	push @next0 , $push_node;	
	}
my $join_the_child_nodes = join ',', @next0;
$swap_string = "($join_the_child_nodes)";
$new_newick_string99 =~ s/$next(\W)/$swap_string$1/;

$nodes_travsd++;if($nodes_travsd =~ /0000$/){
print "nodes_travsd:$nodes_travsd next:$next child_nodes count:$child_nodes childs:@next0\n";
};

my @next1 = @next0;
	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		if($test =~ /INTERNAL_NODE/)
			{
		#####################################
		traverse_backbone_tree_and_print_it($test );
		#####################################
			};
		}


return();

}




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub paint_clades_circular_tree
{



# OLD:

# coloptera test:
# $paint_clades{'INTERNAL_NODE_24975'} = 'pink';
# $paint_clades{'INTERNAL_NODE_3511'} = 'blue';
# $paint_clades{'INTERNAL_NODE_5899'} = 'brown';
# $paint_clades{'INTERNAL_NODE_12263'} = 'red';
# $paint_clades{'INTERNAL_NODE_3303'} = 'green';
# $paint_clades{'INTERNAL_NODE_17429'} = 'orange';
# $paint_clades{'INTERNAL_NODE_15529'} = 'purple';
# coleop actual:
# $paint_clades{'INTERNAL_NODE_215847'} = 'pink';# 215847_Hymenoptera
# $paint_clades{'INTERNAL_NODE_7240'} = 'blue';# 7240_Neuropterida
# $paint_clades{'INTERNAL_NODE_34014'} = 'brown';# 34014_Coleoptera
# $paint_clades{'INTERNAL_NODE_229312'} = 'green';# 229312_Hemiptera
# $paint_clades{'INTERNAL_NODE_90164'} = 'purple';# 90164_Diptera
# $paint_clades{'INTERNAL_NODE_170823'} = 'yellow';# 170823_Lepidoptera
# $paint_clades{'INTERNAL_NODE_6450'} = 'gold';# 6450_Orthopteriodea
# coleop discarede
# $paint_clades{'INTERNAL_NODE_201888'} = 'cyan';# 201888_Microgastrinae
# $paint_clades{'INTERNAL_NODE_10143'} = 'coral';# 10143_Adephaga
# $paint_clades{'INTERNAL_NODE_12546'} = 'aquamarine';# 12546_Scarabaeoidea
# $paint_clades{'INTERNAL_NODE_56131'} = 'chocolate';# 56131_Empidoidea
# $paint_clades{'INTERNAL_NODE_41645'} = 'cornsilk4';# 41645_Syrphoidea
# $paint_clades{'INTERNAL_NODE_85398'} = 'darkgreen';# 85398_Culicomorpha
# $paint_clades{'INTERNAL_NODE_90151'} = 'darkolivegreen3';# 90151_Chironomidae
# $paint_clades{'INTERNAL_NODE_28629'} = 'red';# 28629_Curciloinoidea
# $paint_clades{'INTERNAL_NODE_185111'} = 'orange';# 185111_Apoidea
# 29448 Cucujiformia
# 215847_Hymenoptera 7240_Neuropterida 34014_Coleoptera 28629_Curciloinoidea
# 229312_Hemiptera 185111_Apoidea 90164_Diptera 170823_Lepidoptera
# 6450_Orthopteriodea 201888_Microgastrinae 10143_Adephaga 12546_Scarabaeoidea
# 56131_Empidoidea 41645_Syrphoidea 85398_Culicomorpha 90151_Chironomidae
# 1799_Orthoptera 213544_Hymenoptera 176707_Formicid 7841_Adephaga 31711_Coleoptera 87862_Diptera
# 168521_Lepidoptera
# Psocoptera





# MAR 2017 ITOL2 final:
# $paint_clades{'INTERNAL_NODE_170823'} = 'pink';	# 170823 Neolepidoptera
# $paint_clades{'INTERNAL_NODE_159096'} = 'blue';	# 159096 Noctuoidea
# $paint_clades{'INTERNAL_NODE_130026'} = 'brown';	# 130026 Papilionoidea
# $paint_clades{'INTERNAL_NODE_90165'} = 'green';	# 90165 Diptera
# $paint_clades{'INTERNAL_NODE_90155'} = 'purple';	# 90155 Culicomorpha
# $paint_clades{'INTERNAL_NODE_58128'} = 'yellow';	# 58128 Muscomorpha  ????????????
# $paint_clades{'INTERNAL_NODE_48851'} = 'gold';	# 48852 Calyptratae   51!
# $paint_clades{'INTERNAL_NODE_34015'} = 'cyan';	# 34015 Coleoptera
# $paint_clades{'INTERNAL_NODE_24975'} = 'coral';	# 22132 Chrysomeloidea; 24975!
# $paint_clades{'INTERNAL_NODE_33988'} = 'aquamarine';	# 33988 Elateriformia
# $paint_clades{'INTERNAL_NODE_10061'} = 'chocolate';	# 10061 Adephaga
# $paint_clades{'INTERNAL_NODE_215847'} = 'cornsilk4';	# 215847 Hymenoptera
# $paint_clades{'INTERNAL_NODE_202828'} = 'darkgreen';	# 202828 Ichneumonoidea
# $paint_clades{'INTERNAL_NODE_214075'} = 'darkolivegreen3';	# 214075 Apocrita
# $paint_clades{'INTERNAL_NODE_185111'} = 'red';		# 185111 Apoidea
# $paint_clades{'INTERNAL_NODE_229313'} = 'orange';		# 229313 Hemiptera
# $paint_clades{'INTERNAL_NODE_2288'} = 'lightblue';	# 2288 Palaeoptera
# $paint_clades{'INTERNAL_NODE_77399'} = 'lightgreen';	# Cecidomyiidae




# nothings:
# $paint_clades{'INTERNAL_NODE_'} = 'pink';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'blue';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'brown';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'green';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'purple';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'yellow';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'gold';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'cyan';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'coral';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'aquamarine';	# 
# $paint_clades{'INTERNAL_NODE_'} = 'chocolate';	# 
# $paint_clades{'INTERNAL_NODE_'} = 'cornsilk4';	# 
# $paint_clades{'INTERNAL_NODE_'} = 'darkgreen';	# 
# $paint_clades{'INTERNAL_NODE_'} = 'darkolivegreen3';	# 
# $paint_clades{'INTERNAL_NODE_'} = 'red';		# 
# $paint_clades{'INTERNAL_NODE_'} = 'orange';		# 

$paint_clades{'INTERNAL_NODE_304'} = 'darkgreen';		# 
 $paint_clades{'INTERNAL_NODE_416'} = 'purple';# Geometridae
# $paint_clades{'INTERNAL_NODE_288'} = 'red';	# Crambidae
$paint_clades{'INTERNAL_NODE_'} = 'chocolate';	
$paint_clades{'INTERNAL_NODE_55'} = 'red'; # pyraloid
$paint_clades{'INTERNAL_NODE_421'} = 'cornsilk4';	# bombycoid

# 416:Geometridae Crambidae 288 Noctuoidea 803 Notodontidae



# my @paint_colors = ('pink' ,  'blue' ,  'brown', 'green',  'purple',  'yellow',  'gold',  'cyan', 
#		'coral',  'aquamarine',  'chocolate', 'cornsilk4', 'darkgreen',  'darkolivegreen3', 'red','orange');
# my @paint_states = ('Koppen_A', 'Koppen_B', 'Koppen_C', 'Koppen_D', 'Koppen_E' );
my @paint_Y 	= ('650',      '600',        '550',    '500',      '450','400');

@paint_states = ('Tropical', 'Arid', 'Temperate', 'Cold', 'Polar' );
my @paint_colors = ('green' ,  'red' ,  'cyan',  'brown',  'blue');



if($paint_clades_from_file == 1)
{


open(PAINT_NODES, $paint_nodes_file) || die "\nerror 3293.\n";

open(LOG22, ">$outprefix.process_newick_paintcladeLOG") || die "\nerror 3330\n";


print "opened new ouptut $outprefix.process_newick_paintcladeLOG\n";

my $paint_color_index = 0;
while (my $line = <PAINT_NODES>)
	{
	my $current_color = $paint_colors[$paint_color_index];
	my $current_staet = $paint_states[$paint_color_index];
	$line =~ s/\n//;$line =~ s/\r//; # print "paint_color_index:$paint_color_index current_color:$current_color line:$line\n";
	my @split = split /\s+/ , $line;


#	if($paint_color_index == 4)
	# if you want print all types on one tree:
	if($paint_color_index =~ /\d/)
		{

		foreach my $member(@split)
			{$paint_clades{$member} = $current_color;$clades_defined_for_coloring++};
		print LOG22 "paint_color_index:$paint_color_index current_color:$current_color $#split\n";
		$draw_tree_R_commands6 .=  "text(-550, $paint_Y[$paint_color_index], col = \"$current_color\", labels=\"$current_staet\",cex=10, font=2)\n";
		};



	$paint_color_index++;
	};
close PAINT_NODES;

if($clades_defined_for_coloring >= 1)
	{
	print "user specified \$paint_clades_from_file == 1\n";
	print "count clades defined:$clades_defined_for_coloring\n";
	print "count different colors to be used:$paint_color_index\n";
	}else{
	die "user specified \$paint_clades_from_file == 1, howewver no clades defined from file. quitting.\n";
	};

close LOG22;

};




};# sub paint_clades_circular_tree








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub replace_scientific_notation_in_branchlengths
{
print "
user chosen to replace_scientific_notation_in_branchlengths ....
";






my $tree1= "";
open(IN, $intree) || die "\n\nerror 3236 cant open file: $intree\n\n";
while (my $line = <IN>)
	{if($line =~ /(.+\(.+)/){$tree1 = $1}};
close(IN);
$tree1 =~ s/ //g;$tree_length = length($tree1);

print "\tnewick string of length $tree_length has been read from file:$intree\n\n";

# regex used for couting them elswhere:
# 	while($tree2 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/)


$scientific_notation_branchlengths_replaced = 0;
#while ($tree1 =~ s/\:(\d\.\d+)[eE](\-\d+)/"\:" .  ($1 * ( 10 ** $2 ))/e)
while ($tree1 =~ /\:(\d\.\d+[eE]\-\d+)/)
	{
	my $valueA = $1; my $recalculated = sprintf("%.5f", $valueA);
	$tree1 =~ s/\:\d\.\d+[eE]\-\d+/:$recalculated/;
#	my $valueA = $1; my $valueB = $2; my $recalculated = $valueA * ( 10 ** $valueB );


	$scientific_notation_branchlengths_replaced++;
	if($scientific_notation_branchlengths_replaced =~ /000$/)
		{
		print "scientific_notation_branchlengths_replaced:$scientific_notation_branchlengths_replaced\n";
		print "\twas:$valueA, changed to $recalculated\n";

		};
	};	

print "
scientific_notation_branchlengths_replaced:$scientific_notation_branchlengths_replaced
";

open(OUT_SNBLRM, ">$outprefix.snbl_rm") || die "\nerror 2502\n";
print OUT_SNBLRM "$tree1\n";
close OUT_SNBLRM;







};


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub traverse_subtree
{

my $current_node = $_[0];my $from_parent = $_[1]; my $sum_bl_subtree = $_[2]; 
my $count_connections = $child_counts{$current_node};

# defined child nodes, these are all connecting nodes except that from direction of origin in traversal.
my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	unless($connecting_node =~ /[\w\d]/)
		{die "\n\nerror 330. NODE ($current_node), from_parent:$from_parent " , 
			" count connecting nodes:$count_connections no:$all_connections connecting_node:$connecting_node\n\n"};

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};

my @replace_array = @next1;
my %branchlengths_for_current_children = ();

if($bls_read >= 2)
	{
	my @node_with_bls;
	foreach my $child(@next1)
		{
		$default_branchlength = 0.00001;
		my $bl = $default_branchlength;
		if(exists($branchlengths{$current_node}{$child}))
			{$bl = $branchlengths{$current_node}{$child}}else{die "\ncant find branchlength for $current_node to $child \n"};
		push @node_with_bls, "$child:$bl";
		if(exists($nodes_to_tip_from_current{$child}))
			{
		#	print "have step count for $child\n";
			}else{
		#	print "\nwhy no step count for $child\n"
			};
	# $nodes_to_tip_from_current{$current_node} = $max_steps_to_tip;


		};
	@replace_array = @node_with_bls;
	};



my $join_the_child_nodes = join ',', @replace_array;
$swap_string = "($join_the_child_nodes)";
$subtree_newick_string =~ s/$current_node(\W)/$swap_string$1/;
# $new_newick_string2 =~ s/$current_node(\W)/$swap_string$current_node$1/;

my @next2 = @next1;




for my $index(0 .. $#next2)
	{
	my $test = $next2[$index];
	my $terminals_from_child = $count_terminals_node_lable{$test};

	my $the_branch_length = "NA";
	if(exists($branchlengths{$current_node}{$test}))
		{
		$the_branch_length = $branchlengths{$current_node}{$test};
		}else{
		# hash out if this breaks something:
		die "\nno bl\n";
		};
	my $next_branch_length = $sum_bl_subtree+$the_branch_length;


	if($test =~ /^INTERNAL_NODE_/)
		{
		my $termainsl_from_node = $store_terminals_derived_from_this_node{$test};
		if($termainsl_from_node =~ /\w/)
			{print NODE_DESCRIPTIONS "AddNode $test	$termainsl_from_node\n"};

		# AddNode NODE_2	Psyllobora_vigintiduopunctata_temp5.0 Psyllobora_borealis_temp11.2 Psyllobora_sp_BOLD_AAU2688_temp6.4 Psyllobora_sp_BOLD_ACN4256_temp11.2

		#####################################
		traverse_subtree($test , $current_node , $next_branch_length);#	# recurse
		#####################################
		}else{
		$list_terminals_in_subtree .= "$test ";

		if($next_branch_length >= $max_distance_to_tip_for_subtree){$max_distance_to_tip_for_subtree = $next_branch_length};
		if($next_branch_length <= $min_distance_to_tip_for_subtree){$min_distance_to_tip_for_subtree = $next_branch_length};

		};

	}; # for my $index(0 .. $#next2)



return();

};






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub combine_node_labels
{

#	$combine_node_label_tree = $1; $second_tree_reroot_node = $2; $combine_node_labels = 1;


print "
sub combine_node_labels
	reading node lables from second tree in file:$combine_node_label_tree
";

my $treefile = $combine_node_label_tree;

my $tree11 = "";
open(IN, $treefile) || die "\n\nerror 561 cant open file named:$treefile\n\n";
while (my $line = <IN>)
	{if($line =~ /(.+\(.+)/){$tree11 = $1}};
close(IN);
$tree11 =~ s/ //g;$tree_length11 = length($tree11);

print "\tnewick string of len $tree_length11 read from:$treefile\n\n";


$remove_branchlengths_from_second_tree = 1;
if($remove_branchlengths_from_second_tree == 1)
	{
	print "\n4235. user chosen to remove branch-lengths in input tree ...\n";
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	while($tree11 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/)
		{$count_scientific_notation_branchlengths_removed11++}; 
	# remove regular branchlengths: 0.02048
	while($tree11 =~ s/\:\-*\d+\.\d+//)
		{$count_regular_format_branchlengths_removed11++}; 
	# remove 0 length branchlengths
	while($tree11 =~ s/\:\d+//)
		{$count_zero_length_branchlengths_removed11++};
	print "\t$count_scientific_notation_branchlengths_removed11 count_scientific_notation_branchlengths_removed\n";
	print "\t$count_regular_format_branchlengths_removed11 count_regular_format_branchlengths_removed\n";
	print "\t$count_zero_length_branchlengths_removed11 count_zero_length_branchlengths_removed\n\n";
	if($tree11 =~ /..\:../){die "\nhuh\n"};
	}else{
	print "user chosen to retain branch-lengths in input tree ...\n";
	};


my $newick_string11 	= $tree11;



#############################################################################################################

# third generation Newick reader, mainly to overcome limitation of reading direction. 

print "\nreading Newick tree ....\n";

print "newick_string11:$newick_string11\n";

while ($newick_string11 =~ s/\(([^\(\)]+)\)([0-9\.\_]*)/INTERNAL_NODE_$interal_node11/)
	{# no point reading the adjacent branchlength, 
		# since the identity of the node to which it connects, is not known at this stage,
		# so read lengths to child nodes.

	my $node = $1; my $support = $2;my $nodeID = "INTERNAL_NODE_$interal_node11"; #print "nodeID:$nodeID node:$node\n";
#	print "nodeID:$nodeID support:$support\n";

	my @child_nodes = split /\,/ , $node;
	$child_counts11{$nodeID} = $#child_nodes;#	print "\nnodeID:$nodeID \@child_nodes:@child_nodes\n";
	$node_support11{$nodeID} = $support;#	print "\nnodeID:$nodeID \@child_nodes:@child_nodes\n";
	for $i(0 .. $#child_nodes)
		{
		my $current_child = $child_nodes[$i];
		my $current_child_branchlength = "NA";# length from node to child
		if($current_child =~ s/\:([0-9\.\-E]+)$//) # allows for scientific notation lengths sometimes given:7.59E-4
			{
			$current_child_branchlength = $1;$bls_read++;

			}elsif($current_child =~ /\:/)
			{
			die "\nerror 4785. cant parse branchlength from:$current_child, quitting\n\n"
			};
		
		if($current_child =~ /INTERNAL_NODE_\d+\D/){die "\n4063. newick parse error:$current_child. quitting.\n"};

	#	print "\ti:$i current_child:$current_child branchlength:$current_child_branchlength\n";


		my $current_child_number_connections;		
		if($child_counts{$current_child} =~ /\d/)
			{$current_child_number_connections = $child_counts{$current_child} + 1
			}else{$current_child_number_connections = 0
			};

		# record the child nodes for current node:
		$nodes11{$nodeID}{$i} 			= $current_child;
		$child_counts11{$current_child} = $current_child_number_connections;

		# record length of branch between the current two nodes,
		# store in both directions
		$branchlengths11{$nodeID}{$current_child} = $current_child_branchlength;
		$branchlengths11{$current_child}{$nodeID} = $current_child_branchlength;

		$nodes11{$current_child}{$current_child_number_connections} = $nodeID;
	#	print "\t\tcurrent_child_number_connections:$current_child_number_connections\n";
	#	print "\ti:$i $current_child\n";

		unless($current_child =~ /^INTERNAL_NODE_/)
			{
			$terminal_labels11{$current_child}++;
			if($terminal_labels11{$current_child}>= 2)
				{
				$duplicated_terminal_labels++;
				print "\nwarning, duplicated terminal label:$current_child\n";
				};
			};

		};


	if(length($newick_string) <= 100)
		{
	#	print "newick_string:$newick_string11\n";
		};	


	$root_node11 = $nodeID;
	$interal_node11++;

	};#while ($newick_string =~

#############################################################################################################


print "newick_string11:$newick_string11\n";


my @terminal_array11 = keys %terminal_labels11; @terminal_array11 = sort @terminal_array11;
$count_terminals11 = scalar @terminal_array11;

print "
second tree read, $count_terminals11 terminals
";



if($second_tree_reroot_node =~ /INTERN/)
	{
	$traversal_start_node11 = $second_tree_reroot_node;
#	$count_terminals_node_lable{$traversal_start_node} = $count_terminals;
	}else{
	$traversal_start_node11 = $root_node11; 
	};

print "traverse second tree from $traversal_start_node11\n";



#######################################################################
traverse_second_tree($traversal_start_node11 , "New_Root");#
#######################################################################






};



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub traverse_second_tree
{
my $current_node = $_[0];my $from_parent = $_[1];
my $count_connections = $child_counts11{$current_node};


# defined child nodes, these are all connecting nodes except that from direction of origin in traversal.
my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes11{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};

if($current_node =~ /INTERNAL_NODE/)
	{
	$terminals_belonging_to_current_node_second_tree = "";
	$terminals_belonging_to_current_node_which_have_tax_info = "";

	###############################################################
	get_terminals_from_this_node_second_tree($current_node , $from_parent , 0);#
	###############################################################

	$terminals_belonging_to_current_node_second_tree =~ s/(\t)\t+/$1/;
	$terminals_belonging_to_current_node_second_tree =~ s/^\t+//;$terminals_belonging_to_current_node_second_tree =~ s/\t+$//;
	my @count_terms_array = split /\t/ , $terminals_belonging_to_current_node_second_tree;
	my $count_termnials = scalar @count_terms_array;


	my @sorted_terminals_array = sort @count_terms_array;
	my $join_sorted_line = join ' ' , @sorted_terminals_array;
	my $sec_tree_suport = $node_support11{$current_node};
	$store_second_tree_support{$join_sorted_line} = $sec_tree_suport;

#	if(length($join_sorted_line) <= 60)
		{
	#	print "storin suport $sec_tree_suport node $current_node, termnials:$join_sorted_line\n";
		};

	}else{

	}; # if($current_node =~ /INTERNAL_NODE/)


for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		traverse_second_tree($test , $current_node );#	# recurse
		#####################################
		}
	}

return();

}

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub get_terminals_from_this_node_second_tree
{
my $current_node = $_[0];my $from_parent = $_[1]; my $steps_to_tip = $_[2];$steps_to_tip++;
my $count_connections = $child_counts11{$current_node};
my @next1 = ();my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes11{$current_node}{$all_connections};
	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};
	};


for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];
	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node_second_tree($test , $current_node , $steps_to_tip);#	# recurse
		#####################################
		}elsif($test =~ /[\w\d]/)
		{
		$terminals_belonging_to_current_node_second_tree .= "$test\t";
		};
	}
return();

};#sub get_terminals_from_this_node_second_tree



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




# code from ITOL 1 plotting script, prints most frequent tax for user specfied number of terminal labels



sub print_terminal_tax_labels
{

my @y_positions = keys %store_terminals_at_each_y_position;
@y_positions = sort { $a <=> $b } @y_positions;
my $len_y  = scalar @y_positions;
my $lowest = $y_positions[0];my $highest = $y_positions[$#y_positions];
my $range = $highest-$lowest;
my $bin_size = $range/$num_labs_to_print;
print "
sub print_terminal_labs
	num_labs_to_print:$num_labs_to_print
	number of y positions: $len_y
	lowest:$lowest 
	highest:$highest
	range:$range
	bin_size:$bin_size
";

my $current_position = $lowest;
for my $j(1 .. $num_labs_to_print)
	{
	if($j =~ /00$/){print "$j of $num_labs_to_print\n"};

	my $current_end = $current_position + $bin_size;	#	print "lable no:$j $current_position to $current_end\n";
	my @current_terminals_for_this_y = ();my %count_each_gen = (); my %count_each_taxon = (); my %count_each_label = ();

	foreach my $y(@y_positions)	
		{
		my $terminal = $store_terminals_at_each_y_position{$y};
		if($y >= $current_position && $y <= $current_end)
			{
			push @current_terminals_for_this_y , $terminal; # print "\ty:$y terminal:$terminal\n";
			};
		};


	my $query_assigned_to_this_tip = 0;
	foreach my $term(@current_terminals_for_this_y)
		{
	#	print "term:$term ";
		if($term =~ /^([A-Z][a-z]+)/)# inadeqate, these might not neccessarily be genus
			{
			my $tentative_genus = $1;
			my $test_taxonomy = $complete_lineage_for_this_species{$tentative_genus};#	print "tentative_genus:$tentative_genus test_taxonomy:$test_taxonomy\n";
			if($test_taxonomy =~ / genus\:/)
				{my $geeynus = $tentative_genus;$count_each_gen{$geeynus}++};
			$count_each_taxon{$tentative_genus}++;
			};
		$count_each_label{$term}++;
		};

	my @gen = keys %count_each_gen;my $b6;my $b7; # only those with lineage info
	my @tax2 = keys %count_each_taxon;my $b66;my $b77; # all
	my @labels9 = keys %count_each_label;  my $b666;my $most_freq_label; # all

#	print "stored $#gen genera\n";

	foreach my $g(@gen)
		{
		if($b6 <= $count_each_gen{$g}){$b6 = $count_each_gen{$g};$b7 = $g};
		};
	# b7 will contain most frequent genus name (of those with lineage);  b77 most freq of all. 
	foreach my $g(@tax2)
		{if($b66 <= $count_each_taxon{$g}){$b66 = $count_each_taxon{$g};$b77 = $g}};

	foreach my $g(@labels9)
		{if($b666 <= $count_each_label{$g}){$b666 = $count_each_label{$g};$most_freq_label = $g}};


	$most_freq_label =~ s/([A-Z][a-z]+)_[a-z].+/$1/;

	my $print_col = "black"; # 	if(exists($which_color_is_genus{$b7})){$print_col = $which_color_is_genus{$b7}}else{};
	my $famname = $b7; my $order_name; my $append_tip_string = "";

#	print "\nb7:$b7 b77:$b77\n";

	my $test_taxonomy = $complete_lineage_for_this_species{$b7};
	if($b7 =~ /\w/)
		{
		}else{
		if($b77 =~ /\w/)
			{
			$test_taxonomy = $complete_lineage_for_this_species{$b77}; $b7 = $b77
			}else{
			$b7 = $most_freq_label;
			};
		};

	if($test_taxonomy =~ / order\:(\w\w\w\w)\w+\s/)
		{$order_name = uc($1)};
	if($test_taxonomy =~ / family\:(\w+)\s/)
		{
		$append_tip_string = " (" .  $1 . ")";$famname = $1;
	#	$append_tip_string = " (" . $order_name. " " . $1 . ")";$famname = $1;

		}elsif($test_taxonomy =~ / superfamily\:(\w+)\s/)
		{
		$append_tip_string = " (" . $order_name . " " . $1. ")";$famname = $1;
		}else{
		$append_tip_string = " (NA)";
		};


	if($append_terminal_label_with_fam_name == 1)
		{
	$b7 .= $append_tip_string;
		};

	my $current_mid = $current_position + ($bin_size*0.5);	
	my $y1_proportion = ($current_mid / $count_terminals) + 1; $y1_proportion *= $multipleier;
	my $Xstart 	= $reference_tip_label_x * cos($y1_proportion);# X ax
	my $Ystart 	= $reference_tip_label_x * sin($y1_proportion);# Y ax
	 # wasted 3 hours of my life finding out how to write this one line:
	my $textangle = atan2 ( $Ystart , $Xstart ) * (180 / 3.142) ;

	my $pinttextangle = $textangle;
	if($textangle > 90){$pinttextangle = $textangle - 180};
	if($textangle < -90){$pinttextangle = $textangle + 180};
	
	my $R_command;
	 $R_command = "par(srt=0)\ntext($Xstart, $Ystart,srt=$pinttextangle , col = \"$print_col\", " . 
			"labels=\"$b7\",cex=$tip_tax_label_cex, font=1)\n";# REG

#	 $R_command = "par(srt=0)\ntext($Xstart, $Ystart,srt=$pinttextangle , col = \"$print_col\", " . 
#			"labels=\"$b7\",cex=$tip_tax_label_cex, font=1)\n";# REG


	unless(exists($printed_this_tip_alrady{$b7}))
		{
	#	if($b7 =~ /\w\s\(/)
	#		{
			$draw_tree_R_commands_tip_labels  .= $R_command;
	#		};	
		};

	$dont_print_duplicate_tip_labels_on_reference_tree = 0;

	if($dont_print_duplicate_tip_labels_on_reference_tree == 1)
		{$printed_this_tip_alrady{$b7} = 1}



	$current_position+=$bin_size;
	}



}; # sub print_terminal_tax_labels




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub print_rainbow
{


%gene_synonyms = (
	'Dermaptera' 			=> 'red',
	'Blattodea' 		=> 'green',
	'Odonata' 		=> 'blue',
	'Megaloptera' 		=> 'brown',
	'Trichoptera' 		=> 'orange',
	'Hymenoptera' 		=> 'grey',
	'Hemiptera' 		=> 'pink',
	'Mecoptera' 		=> 'purple',
	'Siphonaptera' 		=> 'antiquewhite',
	'Phasmatodea' 		=> 'aquamarine',
	'Thysanoptera' 		=> 'azure',
	'Plecoptera' 		=> 'beige',
	'Grylloblattodea' 		=> 'bisque',
	'Mantodea' 		=> 'chocolate',
	'Neuroptera' 		=> 'coral',
	'Raphidioptera' 		=> 'cornsilk',
	'Psocoptera' 		=> 'cyan',
	'Zoraptera' 		=> 'darkcyan',
	'Mantophasmatodea' 		=> 'darkred',
	'Ephemeroptera' 		=> 'red4',
	'Coleoptera' 		=> 'darkblue',
	'Strepsiptera' 		=> 'lightblue',
	'Phthiraptera' 		=> 'darkgreen',
	'Embioptera' 		=> 'lightgreen',
	'Orthoptera' 		=> 'antiquewhite3',
	'Lepidoptera' 		=> 'aquamarine3',
	'Archaeognatha' 		=> 'brown4',
	'Isoptera' 		=> 'darkorange',
	'Diptera' 		=> 'pink4',
	'Entomobryomorpha' 	=> 'darkblue'

	);

 %gene_synonyms = (
	'Geometridae' 			=> 'red',
	'Crambidae' 		=> 'green',
	'Noctuoidea' 		=> 'blue',
	'Sphingidae'	=> 'aquamarine3',
	'Limacodidae' 		=> 'pink4'
	);
  






my @y_positions = keys %store_terminals_at_each_y_position;
@y_positions = sort { $a <=> $b } @y_positions;
my $len_y  = scalar @y_positions;
my $lowest = $y_positions[0];my $highest = $y_positions[$#y_positions];
my $range = $highest-$lowest;
my $bin_size = $range/$num_labs_to_print_rainbow;
print "
sub print rainbow
";

# this is very inefficient, would be nice to re-write


my $current_position = $lowest;
for my $j(1 .. $num_labs_to_print_rainbow)
	{
	if($j =~ /00$/){print "$j of $num_labs_to_print_rainbow\n"};
	my $current_end = $current_position + $bin_size;	#	print "lable no:$j $current_position to $current_end\n";
	my @current_terminals_for_this_y = ();my %count_each_gen = (); my %count_each_order = ();my %count_each_taxon = ();

	foreach my $y(@y_positions)	
		{
		my $terminal = $store_terminals_at_each_y_position{$y};
		if($y >= $current_position && $y <= $current_end)
			{
			push @current_terminals_for_this_y , $terminal; # print "\ty:$y terminal:$terminal\n";
			};
		};


	my $query_assigned_to_this_tip = 0;
	foreach my $term(@current_terminals_for_this_y)
		{
		# print "term:$term ";
		if($term =~ /^([A-Z][a-z]+)/)# inadeqate, these might not neccessarily be genus
			{
			my $tentative_genus = $1;
			my $test_taxonomy = $complete_lineage_for_this_species{$tentative_genus};#	print "tentative_genus:$tentative_genus test_taxonomy:$test_taxonomy\n";
			if($test_taxonomy =~ / genus\:/)
				{my $geeynus = $tentative_genus;$count_each_gen{$geeynus}++}

			if($tentative_genus =~ /^[A-Z][a-z]+tera$/)
				{
				$count_each_order{$tentative_genus}++;
				};
$count_each_taxon{$tentative_genus}++;
			};
		};

	my @gen = keys %count_each_gen;my $b6;my $b7;
	my @tax2 = keys %count_each_taxon;my $b66;my $b77;

	foreach my $g(@gen)
		{
		# for rainbow, only get those with order designation
		if($b6 <= $count_each_gen{$g} && $complete_lineage_for_this_species{$g} =~ /\w/)
			{$b6 = $count_each_gen{$g};$b7 = $g}
		};
	# b7 will contain most frequent genus name.
	foreach my $g(@tax2)
		{if($b66 <= $count_each_taxon{$g}){$b66 = $count_each_taxon{$g};$b77 = $g}};



	my $print_col = "black"; # 	if(exists($which_color_is_genus{$b7})){$print_col = $which_color_is_genus{$b7}}else{};
	my $test_taxonomy = $complete_lineage_for_this_species{$b7};
	unless($b7 =~ /\w/){$test_taxonomy = $complete_lineage_for_this_species{$b77}; $b7 = $b77};

	my $famname = $b7; my $order_name; my $append_tip_string = "";
	$bar_color = "grey";

	print "b7:$b7 test_taxonomy:$test_taxonomy\n";

	# try get from formal taxonomic lineage, 
	if($test_taxonomy =~ / order\:(\w+)\s|suborder\:(\w+)\s/)
		{
		$order_name = $1;
		
		my $fam_name;if($test_taxonomy =~ / family\:(\w+)\s/){$fam_name = $1};
		my $superfam_name;if($test_taxonomy =~ / superfamily\:(\w+)\s/){$superfam_name = $1};
		if(exists($gene_synonyms{$order_name}))
			{
			$bar_color = $gene_synonyms{$order_name}
			}elsif(exists($gene_synonyms{$superfam_name}))
			{
			$bar_color = $gene_synonyms{$superfam_name};
			}elsif(exists($gene_synonyms{$fam_name}))
			{
			$bar_color =  $gene_synonyms{$fam_name};
			};
	#	print "order_name:$order_name fam_name:$fam_name bar_color:$bar_color\n";

		}elsif(exists($gene_synonyms{$b7}))
		{
		# 'genus' name might itself be an order name, so try this
		$bar_color = $gene_synonyms{$b7};
		}; # elsif($test_taxonomy =~ / family\:(\w+)\s/)# if not 
	#	{
	#	$bar_color = $gene_synonyms{$b77};
	#	};




	my $current_mid = $current_position + ($bin_size*0.5);	
	my $y1_proportion = ($current_mid / $count_terminals) + 1; $y1_proportion *= $multipleier;
#	my $Xstart 	= $reference_tip_label_x * cos($y1_proportion);# X ax
#	my $Ystart 	= $reference_tip_label_x * sin($y1_proportion);# Y ax



	
	my $R_command;

	# start position.
	my $XstartB 	= $rainbow_X1 * cos($y1_proportion);# X ax
	my $YstartB 	= $rainbow_X1 * sin($y1_proportion);# Y ax
	# end position.
	my $XendB 	= $rainbow_X2 * cos($y1_proportion);# X
	my $YendB 	= $rainbow_X2 * sin($y1_proportion);# Y
	my $R_command =  "segments(" . "$XstartB, $YstartB,$XendB , $YendB , " . 
			"col = \"$bar_color\", lwd = $rainbow_line_width)\n";

	$draw_tree_R_commands_tip_labels  .= $R_command;




	$current_position+=$bin_size;
	}


};

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub proccess_tree_differently
{
my $tree = shift;

unless($tree =~ /\(\w/){die "\nerror, sub process_tree called with no tree.\n"}

# print "tree:$tree\n";



my @terminal_array = keys %terminal_labels; @terminal_array = sort @terminal_array;
$count_terminals = scalar @terminal_array;


foreach my $binomial(@terminal_array)
	{




	my $add_tax_string = "NA_";
	my $look_for_ID = $binomial;

	my $completelineage;my $genus;
	if(exists($complete_lineage_for_this_species{$look_for_ID}))
		{
		$completelineage = $complete_lineage_for_this_species{$look_for_ID};
		}else{
		$genus = $look_for_ID;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
	#	print "no sucess looking for binomial:$binomial, trying genus name:$genus\n";
		$completelineage = $complete_lineage_for_this_species{$genus};
		}

# $binomial $genus

	if($completelineage =~ /\w/)
		{
		$add_tax_string = "";my $tax_added = 0;

		foreach my $current_rank (@ranksarray)
			{
			if($completelineage =~ /($current_rank):(\S+)/i)
				{
				my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
				$rank_assignment_counts{$rank}++;

				if($truncate_new_taxonomic_strings == 1)
					{
					$tax =~ s/^(\w\w\w\w\w\w)\w+/$1/;
					};
				if($uppercase_tax_strings == 1)
					{$tax  =uc($tax)};

				$add_tax_string .= "$tax" . "_";$tax_added++;
				};
			};

		unless($tax_added >= 1)
			{$add_tax_string = "NA_"; 
			if($warning_print_limit_1 < 50)
				{
			print "warning 1759, maybe inapproriate taxonomic ranks selected, because non found in current lineage\n";
			print "\tcompletelineage:$completelineage\n";
				}elsif($warning_print_limit_1 == 50)
				{print "warning 1759 print limit reached\n"};
			$warning_print_limit_1++;
			};

		#$add_tax_string =~ s/_$//;
		$lineages_assigned++;

		$tree =~ s/($binomial[\:\(\)\,])/$add_tax_string$1/;

		unless($add_tax_string =~ /\w/)
			{
		#	foreach my $current_rank (@ranksarray)
		#		{print "\tcurrent_rank:$current_rank\n"};
		#	die "\n\nwhat happened? completelineage:$completelineage\n";
			};

		}else{ # 	if($completelineage =~ /\w/)
		$newick_binomials_for_which_taxonomies_not_found++;	
		print "no lineage found for binomial:$binomial or genus:$genus\n";
		};



};


return($tree);












};

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub color_lineages_from_file
{



open(CLFILE, $color_lineages_file) || die "\nerror 5148. cant find file $color_lineages_file\n";
print "looking for terminals in file $color_lineages_file\n";

my $paint_color_index = 0;
while (my $line = <CLFILE>)
	{
	$line =~ s/^>//;$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^([A-Z][a-z]+).*\t([a-z]+)$/)
		{
		$color_lineage_of_terminal{$1} = $2;
		};

	}
close CLFILE;

my @count_terminals_to_color_lineage = keys %color_lineage_of_terminal;
my @treetermanls = keys %terminal_labels;

print "count_terminals_to_color_lineage:$#count_terminals_to_color_lineage\n";



foreach my $color_lineage_of_tterminal(@count_terminals_to_color_lineage)
	{
#	print "color_lineage_of_representtive of taxon:$color_lineage_of_tterminal\n";	
	my $color_to_which = $color_lineage_of_terminal{$color_lineage_of_tterminal};
	unless($color_to_which =~ /\w/){die "\nerror 5192\n"};

	my $found_tree_terminal_to_color_lineage_of = 0;
	foreach my $tree_term(@treetermanls)
		{
		if($tree_term =~ /^$color_lineage_of_tterminal[_]/)
			{
			$found_tree_terminal_to_color_lineage_of++;
			if($found_tree_terminal_to_color_lineage_of == 1)
				{$color_lineage_of_tree_terminal{$tree_term} = $color_to_which};
		#	print "\tfound thing to color lineage:$tree_term\n";
			};
		};

	
	};


# list of things to color liinega might be genera, and will be muliplt species of each, have now found a single one to color linegge.
# now go through each terminal and record all bracnh labels up to root


my @color_lineage_of_tree_terminals = keys %color_lineage_of_tree_terminal;

foreach my $termina1 (@color_lineage_of_tree_terminals)
	{
	my $which_color = $color_lineage_of_tree_terminal{$termina1};unless($which_color =~ /\w/){die "\nerorr 454565\n"};
#	print "color_lineage_of_terminal:$termina1\n";	
	$color_these_branches{$termina1} = $which_color;


	if($parent_node_IDs{$termina1} =~ /./)
		{
	#	my $parent_ID = $parent_node_IDs{$termina1};
	#	print "\tparent_ID:$parent_ID\n";
		my $parent_ID = $termina1;

		my $X = 0;		
		while ( $X == 0)
			{
			$color_these_branches{$parent_ID} = $which_color;
			$parent_ID = $parent_node_IDs{$parent_ID};
			if($parent_ID =~ /./)
				{
				
				}else{$X = 1};
			};

		}else{
		print "\twarning, cant find parent of terminal:$termina1\n";
		};
		
	
	};



my @color_these_branch = keys %color_these_branches;

print "color_these_branch:$#color_these_branch\n";


};



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






