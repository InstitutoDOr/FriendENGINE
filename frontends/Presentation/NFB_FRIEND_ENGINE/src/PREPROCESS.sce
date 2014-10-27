
#scenario_type = trials;
scenario_type = fMRI_emulation;
#scenario_type = fMRI;

#Parameters required for fMRI scenario type
scan_period = 2000;
pulses_per_scan = 1;      
pulse_code = 53;
sequence_interrupt = false;

#General parameters
default_background_color = 0, 0, 0;
default_picture_duration = next_picture;
default_font_size = 50;

#----------------------------------------------------
begin;
#----------------------------------------------------

picture{ bitmap { filename = "fix.jpg"; };x=0;y=0; } fix;

picture { text { caption = "Initializing..."; } t; } x=0; y=0; } pic;

picture { filename= } fix;

trial {
	picture fix;
	time = 0;
} trialRest;


#STARTING PCL
begin_pcl;

pic.present();

include "communication.pcl";

string sess_id = createSession();



