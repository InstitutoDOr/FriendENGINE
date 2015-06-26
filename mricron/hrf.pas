unit hrf;

interface
uses
    define_types, utypes,metagraph;
const
	 kHRFdur = 24000; //ms for 'full' HRF - window size for HRF
function CreateHRF (lTRsec: double; var lKernelBins: integer; lDefaultsStatsFmriT: integer; var lHRFra, lTDra: doublep): boolean;

function ConvolveTimeCourse(var lTimeCourse: PMatrix; var lKernel: doublep; var l4DTrace: T4DTrace;
lCond,lCondOut,lnVol,lKernelBins,lDefaultsStatsFmriT,lDefaultsStatsFmriT0: integer;
         lTRSec: single; lSliceTime: boolean): boolean;


implementation
uses math {power}, ugamma {gamma},sysutils,dialogs;

const
     kHRFkernelSec = 32;
     //kDefaultsStatsFmriT = 16; //each TR is supersampled at 16x resolution


//from SPM's hrf.m
function spm_Gpdf(x,h,l: double): double;
//emulates spm_Gpdf
begin
	result := power(l,h)*power(x,(h-1))* exp(-l*x);
	result := result / gamma(h);
end;

function fHRF (u,dt: double): double;
//emulates spm_hrf.m
const
	//TR = 1;
	p1= 6; //delay of response
	p2=16;//delay of undershoot (relative to onset)
	p3=1; //dispersion of response
	p4=1; //dispersion of undershoot
	p5=6;  //ratio of response to undershoot
	p7=kHRFkernelSec;//length of kernel (seconds)
begin
    if u <= 0 then
       result := 0
    else
        result   := spm_Gpdf(u,p1/p3,dt/p3) - spm_Gpdf(u,p2/p4,dt/p4)/p5;
end;

function CreateHRF (lTRsec: double; var lKernelBins: integer; lDefaultsStatsFmriT: integer; var lHRFra, lTDra: doublep): boolean;
//NOTE: if this returns TRUE, you MUST freemem lHRFra, lTDra
//returns lHRFra and lTDra with lBins of data - equal to 32sec convolution kernel for
//hemodynamic response (HRF) and the HRF's temporal derivative
var
   lDT,lSum,l1sec: double;
   lI: integer;
begin
     result := false;
     if lDefaultsStatsFmriT < 1 then exit;
     lDT := (lTRsec / lDefaultsStatsFmriT); //DeltaTime - width of each sample in sec
     lKernelBins := round ( kHRFkernelSec / lDT);
     if lKernelBins < 1 then
        exit;
     getmem(lHRFra,lKernelBins*sizeof(double));
     //generate whole HRF kernel
     for lI := 1 to lKernelBins do
            lHRFra^[lI] := fHRF (lI-1,lDT);
     //find sum
     lSum := 0;
     for lI := 1 to lKernelBins do
            lSum := lSum + lHRFra^[lI];
     //normalize - so sum = 1
     for lI := 1 to lKernelBins do
         lHRFra^[lI] := lHRFra^[lI]/lsum;
     //next temporal derivative
     getmem(ltdra,lKernelBins*sizeof(double));
      l1sec := 1/lDT;
     for lI := 1 to lKernelBins do
            ltdra^[lI] := fHRF((lI-1)-l1sec,lDT); //tdHRF (lI-1,lDT);
     //find sum
     lSum := 0;
     for lI := 1 to lKernelBins do
            lSum := lSum + ltdra^[lI];
     //normalize - so sum = 1
     for lI := 1 to lKernelBins do
         ltdra^[lI] := ltdra^[lI]/lsum;
     //temporal derivative is difference between normalized TD and normalized HRF
     for lI := 1 to lKernelBins do
         ltdra^[lI] := lHRFra^[lI]- ltdra^[lI];
     result := true;
end;

function Convolve(var lTimeCoursePrecise,lKernel: doublep; lEventBin,lnVolPrecise,lKernelBins: integer): boolean;
var
   lVol,lStart,lEnd: integer;
begin
    result := false;
    if (lEventBin > lnVolPrecise) then exit; //event too late to influence timecourse
    if ((lEventBin+lKernelBins)< 1) then exit;//event too early to influence timecourse
    lStart := lEventBin;
    if lStart < 1 then
       lStart := 1;
    lEnd := (lEventBin+lKernelBins-1);
    if lEnd > lnVolPrecise then
       lEnd := lnVolPrecise;
    //lOffset := lEventBin;
    for lVol := lStart to lEnd do begin
        lTimeCoursePrecise^[lVol] := lTimeCoursePrecise^[lVol] + lKernel^[lVol -lEventBin+1];
    end;
    result := true;
end;

function ConvolveTimeCourse(var lTimeCourse: PMatrix; var lKernel: doublep; var l4DTrace: T4DTrace; lCond,lCondOut, lnVol,lKernelBins,lDefaultsStatsFmriT,lDefaultsStatsFmriT0: integer;
         lTRSec: single; lSliceTime: boolean): boolean;
var
   lnVolPrecise,lEvent,lVol,lVolx,lEventBin,lEventEnd: integer;
   lDT: double;
   lTimeCoursePrecise: doublep;//supersampled by kDefaultsStatsFmriT
   lAllEvents: boolean;
begin
     result := false;
     if (l4DTrace.Conditions[lCond].Events < 1) or (lnVol < 1) or (lTRSec <= 0) then exit;
     lnVolPrecise := lnVol * lDefaultsStatsFmriT;
     getmem(lTimeCoursePrecise,lnVolPrecise * sizeof(double));
     for lVol := 1 to lnVolPrecise do
         lTimeCoursePrecise^[lVol] := 0;
     lDT := (lTRsec / lDefaultsStatsFmriT); //DeltaTime - width of each sample in sec
     //spm_fmri_design
     //X is supersampled at 16 times (fMRI_T) the number of volumes - with  (32 bin offset)
     //k   = SPM.nscan(s);
     //X = X([0:(k - 1)]*fMRI_T + fMRI_T0 + 32,:);
     for lEvent := 1 to l4DTrace.Conditions[lCond].Events do begin

         lEventBin := round((l4DTrace.Conditions[lCond].EventRA^[lEvent])/lDT);
         //incorrect: same dur will have different number of bins due to rounding:
         //lEventEnd := round((l4DTrace.Conditions[lCond].EventRA^[lEvent]+l4DTrace.Conditions[lCond].DurRA^[lEvent])/lDT);
         //correct: all stimuli of same duration will have identical number of bins
         lEventEnd := lEventBin+round(l4DTrace.Conditions[lCond].DurRA^[lEvent]/lDT);
         //if lEvent = 1 then fx(lEventBin,lEventEnd,l4DTrace.Conditions[lCond].DurRA^[lEvent]);
         repeat
               if (lEventBin > 0) and (lEventBin <= lnVolPrecise) then
                  Convolve(lTimeCoursePrecise,lKernel,lEventBin,lnVolPrecise,lKernelBins);
               inc(lEventBin);
         until lEventBin > lEventEnd;
     end; //for each event
     //output - scaled by reciprocal of DT: e.g. if TR=2, DT=0.125, Scale = 8
     //if TR=2.2, DT=0.1375 Scale = 7.2727
     //this linear scaling does not change any effects - it simply clones SPM2
     lAllEvents := true;
     for lEvent := 1 to l4DTrace.Conditions[lCond].Events do
         if l4DTrace.Conditions[lCond].DurRA^[lEvent] > lDT then
            lAllEvents := false;
     if lAllEvents then
         lDT := 1/lDT
     else
        lDT := 1;
     lVolx := lDefaultsStatsFmriT0;
     for lVol := 1 to lnVol do begin
         if (lVolx > 0) and (lVolx < lnVolPrecise)  then
            lTimeCourse^[lCondOut]^[lVol] := lDT * lTimeCoursePrecise^[lVolx];
         inc(lVolx,lDefaultsStatsFmriT);
     end;
     
     freemem(lTimeCoursePrecise);
     result := true;
end;//func ConvolveTimeCourse

end.
