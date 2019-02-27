
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

picture { 
	text { caption = "Initializing..."; text_align=align_left; } t1; x=-300; y=300; 
} picInit;

picture { 
	bitmap{ filename = "1.jpg"; } bit; x=0; y=0; 
} picFeedback;

bitmap { filename = "fix.jpg"; } bdefault;

picture{
text{ caption = "..."; } t;
x = 0; y = 0;
}pic;

bitmap { filename = "suprime.jpg"; }bsuprime;
bitmap { filename = "experimenta.jpg"; }bexperimenta;
bitmap { filename = "observa.jpg"; }bobserva;
bitmap { filename = "valora.jpg"; }bvalora;


#STARTING PCL
begin_pcl;

int numberOfBlurImages = 10;

bool flip = false;

array<string> codes[35] = {
"experimenta",
"experimenta1",
"valora",
"suprime",
"suprime5",
"valora",
"default",

"experimenta",
"experimenta3",
"valora",
"suprime",
"suprime4",
"valora",
"default",

"experimenta",
"experimenta3",
"valora",
"suprime",
"suprime3",
"valora",
"default",

"experimenta",
"experimenta3",
"valora",
"suprime",
"suprime2",
"valora",
"default",

"experimenta",
"experimenta3",
"valora",
"suprime",
"suprime1",
"valora",
"default"
};


## load feedback pictures
array<bitmap> experimenta_bitmaps[0];
int num_feedback = 20;
int k = 1;
loop until k > num_feedback
begin 
	bitmap b = new bitmap;
	string fname = printf( k, "exp%i.jpg" );
	b.set_filename( fname );
	b.load();
	if flip then
		b.flip_horizontally();
	end;
	
	experimenta_bitmaps.add( b );
	k = k+1;
end;

array<bitmap> suprime_bitmaps[0];
int j;
num_feedback = 20;
k = 1;
int bitmap_index = 0;
loop until k > num_feedback
begin 
	j = 1;
	loop until j > numberOfBlurImages
	begin
	   bitmap b = new bitmap;
	   bitmap_index = k * 100 + (j-1);
	   string fname = printf( bitmap_index, "sup%i.jpg" );
	   b.set_filename( fname );
	   b.load();
	   if flip then
		   b.flip_horizontally();
	   end;
	   suprime_bitmaps.add( b );
	   j = j +1;
	end;
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
	wait_interval(1000);
end;
logfile.add_event_entry("PREPROC_FINISHED");

# start feedback
t.set_caption( "Feedback starting ..." );
t.redraw();
pic.present();

loop until pulse_manager.main_pulse_count() > 0 begin end;

main_request( "NBFEEDBACK" );

t.set_caption( "Waiting for scanning ..." );
t.redraw();
pic.present();
	
array<int> durations[7] = { 4000, 40000, 4000, 4000, 40000, 4000, 10000};

double feedback = 0.0;
int volume = 0;
int volume_value_received = 0;
int ci = 0;
int bi = 1;
int subindex = 1;
int bitmap_stop_time = 0;
int queryDelay = 700;
int lastQuery = clock.time();
string classnumber = "";

int volume_delay = pulse_manager.main_pulse_count();
int start_time = clock.time();

loop until false
begin
	
	volume = pulse_manager.main_pulse_count() - volume_delay + 1;

	# mudar bitmaps conforme protocolo
	if clock.time() - start_time > bitmap_stop_time  then
		ci = ci + 1;
		
		if ci > codes.count() then break; end;
		
		if codes[ci] == "suprime" then
			picFeedback.set_part(1, bsuprime );
		elseif codes[ci] == "valora" then
			picFeedback.set_part(1,bvalora );
		elseif codes[ci] == "experimenta" then
			picFeedback.set_part(1,bexperimenta );
		elseif codes[ci] == "default" then
			picFeedback.set_part(1,bdefault );
		else
			subindex = int( codes[ci].substring(codes[ci].count(), 1));
 
			string condition = codes[ci].substring(1,codes[ci].count()-1);

			if condition == "suprime" then
				bi = subindex * numberOfBlurImages + 1;
			   picFeedback.set_part(1, suprime_bitmaps[bi] );
			   logfile.add_event_entry("presenting_" + condition + "_" + string(bi));
			elseif condition  == "experimenta" then
				bi = subindex;
			   picFeedback.set_part(1, experimenta_bitmaps[bi]);
			   logfile.add_event_entry("presenting_" + condition + "_" + string(bi));
			end;
		end;
		
		picFeedback.present();
		logfile.add_event_entry(codes[ci]);
		bitmap_stop_time = bitmap_stop_time + durations[mod(ci-1,7)+1];
		logfile.add_event_entry("new_stop_" + string(bitmap_stop_time));
	end;

	classnumber = "";
	
	if clock.time()-lastQuery > queryDelay then
	   if volume > 0 && volume > volume_value_received then
		   feedback = getFeedbackValue( volume , classnumber );
		   if classnumber != "0" then
		      volume_value_received = volume;
		   end;
	   end;
	
	   # if not finished yet, try the anterior value
	   if volume > 1 && volume-1 > volume_value_received then
		   feedback = getFeedbackValue( volume-1 , classnumber );
		   if classnumber != "0" then
		      volume_value_received = volume-1;
		   end;
	   end;
	
	   /*
	   # if not finished yet, try the anterior value
	   if classnumber == "0" && volume > 2 && volume-2 > volume_value_received then
		   feedback = getFeedbackValue( volume-2, classnumber );
		   if classnumber != "0" then
		      volume_value_received = volume-2;
		   end;
	   end;
	   */
	   lastQuery = clock.time();
	end;

	if classnumber == "4" then
		if feedback < 0.0 then feedback = 0.0; end;
		if feedback > 1.0 then feedback = 1.0; end;
		
		int index = subindex * numberOfBlurImages + int(feedback * 9) + 1;
      picFeedback.set_part(1, suprime_bitmaps[index] );		
		picFeedback.present();
   end;
end;

main_request( "ENDSESSION" );

