
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
default_font_size = 30;

#----------------------------------------------------
begin;
#----------------------------------------------------

picture{ bitmap { filename = "fix.jpg"; } bit_fix; x=0;y=0; } fix;

picture { 
	text { caption = "Initializing..."; text_align=align_left; } t; x=-300; y=300; 
} pic;

picture { 
	bitmap{ filename = "fix.jpg"; } bit; x=0; y=0; 
} picFeedback;

trial {
	picture fix;
	time = 0;
} trialRest;


#STARTING PCL
begin_pcl;

bool flip = true;

## load feedback pictures
array<bitmap> feedback_bitmap[0];
int num_feedback = 47;
int k = 1;
loop until k > num_feedback
begin 
	bitmap b = new bitmap;
	string fname = printf( k, "%03i.jpg" );
	b.set_filename( fname );
	b.load();
	if flip then
		b.flip_horizontally();
	end;
	
	feedback_bitmap.add( b );
	k = k+1;
end;

## declare global variable sess_id
string sess_id;

pic.present();

include "communication.pcl";

# open once the main thread
open_mainThread();

t.set_caption( "Creating session ..." );
t.redraw();
pic.present();

sess_id = createSession();

t.set_caption( "Sending plugin information ..." );
t.redraw();
pic.present();

sendPluginInformation();

# start preprocessing step
t.set_caption( "Preprocessing ..." );
t.redraw();
pic.present();

main_request( "NBPREPROC" );

string preproc_ok = "";
loop until preproc_ok == "OK"
begin
	preproc_ok = getCommandState( "PREPROC" );
end;
logfile.add_event_entry("PREPROC_FINISHED");

# start feedback
t.set_caption( "Feedback starting ..." );
t.redraw();
pic.present();

main_request( "NBFEEDBACK" );

t.set_caption( "Waiting for scanning ..." );
t.redraw();
pic.present();

int volume_delay = pulse_manager.main_pulse_count();

double feedback = 0.0;
int volume = 0;
int volume_value_received = 0;
int volumeMax = 296;
string feedback_ok = "";
loop until (volume > volumeMax)
begin
	
	volume = pulse_manager.main_pulse_count() - volume_delay;

	feedback_ok = "";
	
	if volume > 1 && volume-1 > volume_value_received then
		feedback = getFeedbackValue( volume-1 , feedback_ok );
		if feedback_ok != "0" then
		   volume_value_received = volume-1;
		end;
	end;
	
	# if not finished yet, try the anterior value
	if feedback_ok == "0" && volume > 2 && volume-2 > volume_value_received then
		feedback = getFeedbackValue( volume-2, feedback_ok );
		if feedback_ok != "0" then
		   volume_value_received = volume-2;
		end;
	end;
	

	if feedback_ok == "1" then
		picFeedback.set_part( 1, bit_fix );
		picFeedback.present();
	elseif feedback_ok != "0" && feedback_ok != "" then
		if feedback < 0.0 then feedback = 0.0; end;
		if feedback > 1.0 then feedback = 1.0; end;
				
		int ind = int( feedback * double(num_feedback-1) ) + 1; 
		picFeedback.set_part( 1, feedback_bitmap[ind] );
		picFeedback.present();
	end;


end;

main_request( "GLM" );
main_request( "FEATURESELECTION" );
main_request( "ENDSESSION" );

