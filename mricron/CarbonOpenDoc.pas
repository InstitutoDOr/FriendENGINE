unit CarbonOpenDoc;


interface
{$H+}
uses
  //FPCMacOSAll,
  MacOSAll, CarbonProc;
  
//type
//  TFormOpenFileMethod = procedure(const FileName : string) of object;
//procedure InitOpenDocHandler(MethodToUse : TFormOpenFileMethod);
procedure InitOpenDocHandler;


implementation
uses
 nifti_img_view;
//var
//  OpenFileMethod : TFormOpenFileMethod;

function OpenDocEventHandler(var theAppleEvent: AppleEvent;
                             var reply: AppleEvent;
                                 handlerRefcon: SInt32): OSErr; stdcall;
var
  DocList: AEDescList;
  FileCount: Integer;
  FileIdx: Integer;
  Keyword: AEKeyword;
  FileDesc: AEDesc;
  FileRef: FSRef;
  FileURL: CFURLRef;
  FileCFStr: CFStringRef;
begin
  if OSError(AEGetParamDesc(theAppleEvent, keyDirectObject, typeAEList, DocList),
             'OpenDocEventHandler', '', 'AEGetParamDesc') then
    Exit;

  try
    if OSError(AECountItems(DocList, FileCount),
               'OpenDocEventHandler', '', 'AECountItems') then
      Exit;

    for FileIdx := 1 to FileCount do
      begin

      if OSError(AEGetNthDesc(DocList, FileIdx, typeFSRef, @Keyword, FileDesc), 
                 'OpenDocEventHandler', '', 'AEGetNthDesc') then
        Exit;

      if OSError(AEGetDescData(FileDesc, @FileRef, SizeOf(FSRef)), 
                 'OpenDocEventHandler', '', 'AEGetDescData') then 
        Exit;

      if OSError(AEDisposeDesc(FileDesc), 
                 'OpenDocEventHandler', '', 'AEDisposeDesc') then 
        Exit;

      FileURL := CFURLCreateFromFSRef(kCFAllocatorDefault, FileRef);
      FileCFStr := CFURLCopyFileSystemPath(FileURL, kCFURLPOSIXPathStyle);
      ImgForm.FormOpenFileMethod(CFStringToStr(FileCFStr));
      //ImgForm.OpenFileMethod(CFStringToStr(FileCFStr));

      FreeCFString(FileURL);
      FreeCFString(FileCFStr);
      end;
  finally
    AEDisposeDesc(DocList);
  end;

end;  {OpenDocEventHandler}


procedure InitOpenDocHandler {(MethodToUse : TFormOpenFileMethod)};
begin
  //OpenFileMethod := MethodToUse;
  
  AEInstallEventHandler(kCoreEventClass, kAEOpenDocuments,
                        NewAEEventHandlerUPP(
                         AEEventHandlerProcPtr(Pointer(@OpenDocEventHandler))),
                        0, False);

end;  {InitOpenDocHandler}


end.

