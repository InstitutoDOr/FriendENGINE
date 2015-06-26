program mricron;

{$mode objfpc}{$H+}
uses
  {$IFDEF UNIX}
  cthreads,
  {$ENDIF}
     Interfaces, // this includes the LCL widgetset
  Forms, nifti_img_view, nifti_hdr_view,
  about, Text, ReadInt, histoform, autoroi, ROIfilt, render,
   MultiSlice, CropEdges, bet, mni,
  voismooth, prefs, perisettings, graphx, cutout, ReadFloat, landmarks, 
batchstatselect, nii_label;
{$IFNDEF UNIX}
 {$IFDEF FPC}
   {$R manifest.res}
{$ELSE}
  {$R *.res}//windows icon
 {$ENDIF}
{$ENDIF}



{$IFDEF WINDOWS}{$R mricron.rc}{$ENDIF}

begin
  Application.Title:='MRIcron';
  Application.Initialize;
  Application.CreateForm(TImgForm, ImgForm);
  Application.CreateForm(THdrForm, HdrForm);
  Application.CreateForm(TAnatForm, AnatForm);
  Application.CreateForm(TAboutForm, AboutForm);
  Application.CreateForm(TTextForm, TextForm);
  Application.CreateForm(TReadIntForm, ReadIntForm);
  Application.CreateForm(TAutoROIForm, AutoROIForm);
  Application.CreateForm(THistogramForm, HistogramForm);
  Application.CreateForm(TFilterROIform, FilterROIform);
  Application.CreateForm(TMultiSliceForm, MultiSliceForm);
  Application.CreateForm(TRenderForm, RenderForm);
  Application.CreateForm(TCropEdgeForm, CropEdgeForm);
  Application.CreateForm(TBETForm, BETForm);
  Application.CreateForm(TMNIForm, MNIForm);
  Application.CreateForm(Tvoismoothform, voismoothform);
  Application.CreateForm(TPrefForm, PrefForm);
  Application.CreateForm(TPSForm, PSForm);
  Application.CreateForm(TGraph4DForm, Graph4DForm);
  Application.CreateForm(TCutoutForm, CutoutForm);
  Application.CreateForm(TReadFloatForm, ReadFloatForm);
  Application.Run;
end.
