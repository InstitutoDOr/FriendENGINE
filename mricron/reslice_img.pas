unit reslice_img;
//12 April 2009 - added lTrilinearSmooth option to allow nearest neighbor interpolation
interface
uses
{$ifndef fpc}{windows,} {$endif}
GraphicsMathLibrary,nifti_hdr, nifti_types;
function Reslice_Img_To_Unaligned (var lTargHdr: TNIfTIhdr; var lSrcHdr: TMRIcroHdr; lTrilinearSmoothIn: boolean): boolean;
function Hdr2InvMat (lHdr: TNiftiHdr; var lOK: boolean): TMatrix;
procedure Voxel2mm(var X,Y,Z: single; var lHdr: TNIfTIHdr);
procedure mm2Voxel (var X,Y,Z: single; var lInvMat: TMatrix);
implementation



uses  dialogs, define_types;


function Hdr2Mat (lHdr:  TNIFTIhdr): TMatrix;
begin
  Result := Matrix3D (
  lHdr.srow_x[0],lHdr.srow_x[1],lHdr.srow_x[2],lHdr.srow_x[3],      // 3D "graphics" matrix
  lHdr.srow_y[0],lHdr.srow_y[1],lHdr.srow_y[2],lHdr.srow_y[3],      // 3D "graphics" matrix
  lHdr.srow_z[0],lHdr.srow_z[1],lHdr.srow_z[2],lHdr.srow_z[3],      // 3D "graphics" matrix
						   0,0,0,1);
end;

(*procedure ReportMatrix (lM:TMatrix);
const
	kCR = chr (13);
begin
	showmessage(RealToStr(lM.matrix[1,1],6)+','+RealToStr(lM.matrix[1,2],6)+','+RealToStr(lM.matrix[1,3],6)+','+RealToStr(lM.matrix[1,4],6)+kCR+
		RealToStr(lM.matrix[2,1],6)+','+RealToStr(lM.matrix[2,2],6)+','+RealToStr(lM.matrix[2,3],6)+','+RealToStr(lM.matrix[2,4],6)+kCR+
		RealToStr(lM.matrix[3,1],6)+','+RealToStr(lM.matrix[3,2],6)+','+RealToStr(lM.matrix[3,3],6)+','+RealToStr(lM.matrix[3,4],6)+kCR
    +RealToStr(lM.matrix[4,1],6)+','+RealToStr(lM.matrix[4,2],6)+','+RealToStr(lM.matrix[4,3],6)+','+RealToStr(lM.matrix[4,4],6)
	  );
end;     *)

(*
procedure  SPMmat(var lDestMat: TMatrix);
//SPM matrices are indexed from 1
//This function is only useful for direct comparisons with SPM
var
  lTemp,lVS: TMatrix;
begin
  lVS := Matrix3D (1,0,0,-1,
    0,1,0,-1,
    0,0,1,-1, 0,0,0,1);//VoxelShift
  lTemp := lDestMat;
  lDestMat := MultiplyMatrices(lTemp,lVS);
end;*)

procedure  Coord(var lV: TVector; var lMat: TMatrix);
//transform X Y Z by matrix
var
  lXi,lYi,lZi: single;
begin
  lXi := lV.x; lYi := lV.y; lZi := lV.z;
  lV.x := (lXi*lMat.matrix[1][1]+lYi*lMat.matrix[1][2]+lZi*lMat.matrix[1][3]+lMat.matrix[1][4]);
  lV.y := (lXi*lMat.matrix[2][1]+lYi*lMat.matrix[2][2]+lZi*lMat.matrix[2][3]+lMat.matrix[2][4]);
  lV.z := (lXi*lMat.matrix[3][1]+lYi*lMat.matrix[3][2]+lZi*lMat.matrix[3][3]+lMat.matrix[3][4]);
end;

procedure  Transposemat(var lMat: TMatrix);
var
  lTemp: TMatrix;
  i,j: integer;
begin
  lTemp := lMat;
  for i := 1 to lMat.size do
    for j := 1 to lMat.size do
      lMat.matrix[i,j] := lTemp.matrix[j,i];
end;

function gaussj(VAR a: TMatrix): boolean;//Invert a Matrix - see Numerical Recipes
label
  666;
VAR
   big,dum,pivinv: real;
   n,i,icol,irow,j,k,l,ll: integer;
   indxc,indxr,ipiv: array [1..4] of integer;
BEGIN
   result := true;
   icol := 1;//not used - avoids compiler warning
   irow := 1;//not used - avoids compiler warning
   n := a.size;
   FOR j := 1 TO n DO BEGIN
      ipiv[j] := 0
   END;
   FOR i := 1 TO n DO BEGIN
      big := 0.0;
      FOR j := 1 TO n DO BEGIN
         IF (ipiv[j] <> 1) THEN BEGIN
            FOR k := 1 TO n DO BEGIN
               IF (ipiv[k] = 0) THEN BEGIN
                  IF (abs(a.matrix[j,k]) >= big) THEN BEGIN
                     big := abs(a.matrix[j,k]);
                     irow := j;
                     icol := k
                  END
               END ELSE IF (ipiv[k] > 1) THEN BEGIN
                  goto 666;
               END
            END
         END
      END;
      ipiv[icol] := ipiv[icol]+1;
      IF (irow <> icol) THEN BEGIN
         FOR l := 1 TO n DO BEGIN
            dum := a.matrix[irow,l];
            a.matrix[irow,l] := a.matrix[icol,l];
            a.matrix[icol,l] := dum
         END;
      END;
      indxr[i] := irow;
      indxc[i] := icol;
      IF (a.matrix[icol,icol] = 0.0) THEN
        goto 666;
      pivinv := 1.0/a.matrix[icol,icol];
      a.matrix[icol,icol] := 1.0;
      FOR l := 1 TO n DO BEGIN
         a.matrix[icol,l] := a.matrix[icol,l]*pivinv
      END;
      FOR ll := 1 TO n DO BEGIN
         IF (ll <> icol) THEN BEGIN
            dum := a.matrix[ll,icol];
            a.matrix[ll,icol] := 0.0;
            FOR l := 1 TO n DO BEGIN
               a.matrix[ll,l] := a.matrix[ll,l]-a.matrix[icol,l]*dum
            END;
         END
      END
   END;
   FOR l := n DOWNTO 1 DO BEGIN
      IF (indxr[l] <> indxc[l]) THEN BEGIN
         FOR k := 1 TO n DO BEGIN
            dum := a.matrix[k,indxr[l]];
            a.matrix[k,indxr[l]] := a.matrix[k,indxc[l]];
            a.matrix[k,indxc[l]] := dum
         END
      END
   END;
   exit;
   666: //only get here if there is an error
   Showmessage('error in reslice_img - singular matrix. Spatial orientation is ambiguous.');
   a := Eye3D;
   result := false;
END;

procedure SubVec (var lVx: TVector; lV0: TVector);
begin
  lVx.x := lVx.x - lV0.x;
  lVx.y := lVx.y - lV0.y;
  lVx.z := lVx.z - lV0.z;
end;

(*procedure mm2Voxel (var X,Y,Z: single; var lInvMat: TMatrix);
//returns voxels indexed from 1 not 0!
var
   lV: TVector;
   lSrcMatInv,lSrcMat: TMatrix;
begin
     lV := Vector3D (X,Y,Z);
     lV := Transform (lV,lInvMat);
     X := lV.x+1;
     Y := lV.y+1;
     Z := lV.z+1;
end;*)

procedure mm2Voxel (var X,Y,Z: single; var lInvMat: TMatrix);
//returns voxels indexed from 1 not 0!
var
   lV: TVector;
   lSrcMatInv,lSrcMat: TMatrix;
begin
     lV := Vector3D (X,Y,Z);
     Coord (lV,lInvMat);
     X := lV.x+1;
     Y := lV.y+1;
     Z := lV.z+1;
end;

procedure Voxel2mm(var X,Y,Z: single; var lHdr: TNIfTIHdr);
var
   lV: TVector;
   lMat: TMatrix;
begin
     //lV := Vector3D (X-1,Y-1,Z-1);
     lV := Vector3D (X-1,Y-1,Z-1);
     lMat := Hdr2Mat(lHdr);
     Coord(lV,lMat);
     X := lV.x;
     Y := lV.y;
     Z := lV.z;
end;

function Voxel2Voxel (var lDestHdr,lSrcHdr: TNIFTIhdr): TMatrix;
//returns matrix for transforming voxels from one image to the other image
//results are in VOXELS not mm
var
   lV0,lVx,lVy,lVz: TVector;
   lDestMat,lSrcMatInv,lSrcMat: TMatrix;

begin
     //Step 1 - compute source coordinates in mm for 4 voxels
     //the first vector is at 0,0,0, with the
     //subsequent voxels being left, up or anterior
     lDestMat := Hdr2Mat(lDestHdr);
     //SPMmat(lDestMat);
     lV0 := Vector3D  (0,0,0);
     lVx := Vector3D  (1,0,0);
     lVy := Vector3D  (0,1,0);
     lVz := Vector3D  (0,0,1);
     Coord(lV0,lDestMat);
     Coord(lVx,lDestMat);
     Coord(lVy,lDestMat);
     Coord(lVz,lDestMat);
     lSrcMat := Hdr2Mat(lSrcHdr);
     //SPMmat(lSrcMat);
     lSrcMatInv := lSrcMat;
     gaussj(lSrcMatInv);
     //the vectors should be rows not columns....
     //therefore we transpose the matrix
     Transposemat(lSrcMatInv);
     //the 'transform' multiplies the vector by the matrix
     lV0 := Transform (lV0,lSrcMatInv);
     lVx := Transform (lVx,lSrcMatInv);
     lVy := Transform (lVy,lSrcMatInv);
     lVz := Transform (lVz,lSrcMatInv);
     //subtract each vector from the origin
     // this reveals the voxel-space influence for each dimension
     SubVec(lVx,lV0);
     SubVec(lVy,lV0);
     SubVec(lVz,lV0);
     result := Matrix3D(lVx.x,lVy.x,lVz.x,lV0.x,
      lVx.y,lVy.y,lVz.y,lV0.y,
      lVx.z,lVy.z,lVz.z,lV0.z, 0,0,0,1);
end;

procedure CopyHdrMat(var lTarg,lDest: TNIfTIHdr);
//destination has dimensions and rotations of destination
var
   lI: integer;
begin
     //destination will have dimensions of target
   lDest.dim[0] := 3; //3D
   for lI := 1 to 3 do
       lDest.dim[lI] := lTarg.dim[lI];
   lDest.dim[4] := 1; //3D
   //destination will have pixdim of target
   for lI := 0 to 7 do
       lDest.pixdim[lI] := lTarg.pixdim[lI];
   lDest.xyzt_units := lTarg.xyzt_units; //e.g. mm and sec
   lDest.qform_code := lTarg.qform_code;
   lDest.sform_code := lTarg.sform_code;
   lDest.quatern_b := lTarg.quatern_b;
   lDest.quatern_c := lTarg.quatern_c;
   lDest.quatern_d := lTarg.quatern_d;
   lDest.qoffset_x := lTarg.qoffset_x;
   lDest.qoffset_y := lTarg.qoffset_y;
   lDest.qoffset_z := lTarg.qoffset_z;
   for lI := 0 to 3 do begin
       lDest.srow_x[lI] := lTarg.srow_x[lI];
       lDest.srow_y[lI] := lTarg.srow_y[lI];
       lDest.srow_z[lI] := lTarg.srow_z[lI];
   end;
end;

function OneToOne(lM:TMatrix): boolean;
var
  lC,lR: integer;
begin
  result := false;
  for lC := 1 to 3 do
    for lR := 1 to 3 do
      if (lM.matrix[lC,lR] <> 0) and ((abs(lM.matrix[lC,lR])- 1) > 0.00001) then
        exit;
  result := true;
end;

 function Reslice_Img_To_Unaligned (var lTargHdr: TNIfTIhdr; var lSrcHdr: TMRIcroHdr; lTrilinearSmoothIn: boolean): boolean;
var
  lXrM1,lYrM1,lZrM1,lZx,lZy,lZz,lYx,lYy,lYz,lXreal,lYreal,lZreal: single;
 lXo,lYo,lZo,lMinY,lMaxY,lMinZ,lMaxZ,
 lPos,lXs,lYs,lZs,lXYs,lXYZs,lX,lY,lZ,lOutVolItems,
 lXi,lYi,lZi: integer;
 lDestHdr: TNIFTIhdr;
 lMat: TMatrix;
 lTrilinearSmooth,lOverlap: boolean;
 lXx,lXy,lXz: Singlep0;
 l32fs,l32f : SingleP;
 l16is,l16i : SmallIntP;
 l8i,l8is,lSrcBuffer,lBuffUnaligned,lBuffAligned,lBuffOutUnaligned: bytep;
begin
  lTrilinearSmooth := lTrilinearSmoothIn;
     result := false;
     lOverlap := false;
     lDestHdr := lSrcHdr.NIfTIHdr; //destination has the comments and voxel BPP of source
     CopyHdrMat(lTargHdr,lDestHdr);//destination has dimensions and rotations of destination
     lXs := lSrcHdr.NIfTIHdr.Dim[1];
     lYs := lSrcHdr.NIfTIHdr.Dim[2];
     lZs := lSrcHdr.NIfTIHdr.Dim[3];

     lXYs:=lXs*lYs; //slicesz
     lXYZs := lXYs*lZs;
     lX := lDestHdr.Dim[1];
     lY := lDestHdr.Dim[2];
     lZ := lDestHdr.Dim[3];
     lOutVolItems :=lX*lY*lZ;
    if lSrcHdr.ImgBufferBPP = 4 then begin
	    l32fs := SingleP(lSrcHdr.ImgBuffer);
	    GetMem(lBuffOutUnaligned,(lOutVolItems*sizeof(single))+16);
	    {$IFDEF FPC}
	    l32f := align(lBuffOutUnaligned,16);
   	    {$ELSE}
	    l32f := SingleP($fffffff0 and (integer(lBuffOutUnaligned)+15));
            {$ENDIF}
	    for lPos := 1 to lOutVolItems do
		    l32f^[lPos] := 0; //set all to zero
    end else if lSrcHdr.ImgBufferBPP = 2 then begin
	    l16is := SmallIntP(lSrcHdr.ImgBuffer);
	    GetMem(lBuffOutUnaligned,(lOutVolItems*sizeof(smallint))+16);
	       {$IFDEF FPC}
	    l16i := align(lBuffOutUnaligned,16);
   {$ELSE}
	    l16i := SmallIntP($fffffff0 and (integer(lBuffOutUnaligned)+15));
   {$ENDIF}
	    for lPos := 1 to lOutVolItems do
		    l16i^[lPos] := 0; //set all to zero
    end else if lSrcHdr.ImgBufferBPP = 1 then begin
      l8is := ByteP(lSrcHdr.ImgBuffer);
	    GetMem(l8i,lOutVolItems);
	    Fillchar(l8i^,lOutVolItems,0); //set all to zero
    end;
    lMat := Voxel2Voxel (lTargHdr,lSrcHdr.NIfTIHdr);
    //lDestHdr := lSrcHdr; //destination has the comments and voxel BPP of source
    //CopyHdrMat(lTargHdr,lDestHdr);//destination has dimensions and rotations of destination
     //now we can apply the transforms...
     //build lookup table - speed up inner loop
     getmem(lXx, lX*sizeof(single));
     getmem(lXy, lX*sizeof(single));
     getmem(lXz, lX*sizeof(single));
     for lXi := 0 to (lX-1) do begin
      lXx^[lXi] := lXi*lMat.matrix[1][1];
      lXy^[lXi] := lXi*lMat.matrix[2][1];
      lXz^[lXi] := lXi*lMat.matrix[3][1];
     end;
     lPos := 0;
     if (lTrilinearSmooth) and (OneToOne(lMat)) then begin
        lTrilinearSmooth := false;
     end;
if lTrilinearSmooth then begin//compute trilinear interpolation
//compute trilinear interpolation
     for lZi := 0 to (lZ-1) do begin
         //these values are the same for all voxels in the slice
         // compute once per slice
         lZx := lZi*lMat.matrix[1][3];
         lZy := lZi*lMat.matrix[2][3];
         lZz := lZi*lMat.matrix[3][3];
         for lYi := 0 to (lY-1) do begin
             //these values change once per row
             // compute once per row
             lYx :=  lYi*lMat.matrix[1][2];
             lYy :=  lYi*lMat.matrix[2][2];
             lYz :=  lYi*lMat.matrix[3][2];
             for lXi := 0 to (lX-1) do begin
                 //compute each column
                 inc(lPos);
                 lXreal := (lXx^[lXi]+lYx+lZx+lMat.matrix[1][4]);
                 lYreal := (lXy^[lXi]+lYy+lZy+lMat.matrix[2][4]);
                 lZreal := (lXz^[lXi]+lYz+lZz+lMat.matrix[3][4]);
                 //need to test Xreal as -0.01 truncates to zero
                 if (lXreal >= 0) and (lYreal >= 0) and (lZreal >= 0) and
                     (lXreal < (lXs -1)) and (lYreal < (lYs -1) ) and (lZreal <= (lZs -1))   //June09 lZReal <= instead of <
                  then begin
                    //compute the contribution for each of the 8 source voxels
                    //nearest to the target
                                     lOverlap := true;
                                     lXo := trunc(lXreal);
			              lYo := trunc(lYreal);
			              lZo := trunc(lZreal);
			              lXreal := lXreal-lXo;
			              lYreal := lYreal-lYo;
			              lZreal := lZreal-lZo;
                    lXrM1 := 1-lXreal;
			              lYrM1 := 1-lYreal;
			              lZrM1 := 1-lZreal;
			              lMinY := lYo*lXs;
			              lMinZ := lZo*lXYs;
			              lMaxY := lMinY+lXs;
                    inc(lXo);//images incremented from 1 not 0
                    //Check if sample is perfectly in the Z-plane.
                    //This requires only 8 samples, so its faster
                    //in addition, for very thin volumes, it allows us to sample to the edge
                    if lZReal = 0 then begin // perfectly in plane, only sample 4 voxels near each other
                      case lSrcHdr.ImgBufferBPP  of
                        1 : l8i^[lPos] :=
                          round ( ( (lXrM1*lYrM1)*l8is^[lXo+lMinY+lMinZ])+((lXreal*lYrM1)*l8is^[lXo+1+lMinY+lMinZ])+((lXrM1*lYreal)*l8is^[lXo+lMaxY+lMinZ])+((lXreal*lYreal)*l8is^[lXo+1+lMaxY+lMinZ]));
 	                      2: l16i^[lPos] :=
                           round (( (lXrM1*lYrM1)*l16is^[lXo+lMinY+lMinZ])+((lXreal*lYrM1)*l16is^[lXo+1+lMinY+lMinZ])+((lXrM1*lYreal)*l16is^[lXo+lMaxY+lMinZ])+((lXreal*lYreal)*l16is^[lXo+1+lMaxY+lMinZ]));
	                      4: l32f^[lPos] :=
		 	                      ( (lXrM1*lYrM1)*l32fs^[lXo+lMinY+lMinZ])+((lXreal*lYrM1)*l32fs^[lXo+1+lMinY+lMinZ])+((lXrM1*lYreal)*l32fs^[lXo+lMaxY+lMinZ])+((lXreal*lYreal)*l32fs^[lXo+1+lMaxY+lMinZ]);
                      end; //case
                    end else begin  //not perfectly in plane... we need 8 samples...
                      lMaxZ := lMinZ+lXYs;
                      case lSrcHdr.ImgBufferBPP  of
                          1 : l8i^[lPos] :=
                           round ({all min} ( (lXrM1*lYrM1*lZrM1)*l8is^[lXo+lMinY+lMinZ])
                           {x+1}+((lXreal*lYrM1*lZrM1)*l8is^[lXo+1+lMinY+lMinZ])
                           {y+1}+((lXrM1*lYreal*lZrM1)*l8is^[lXo+lMaxY+lMinZ])
                           {z+1}+((lXrM1*lYrM1*lZreal)*l8is^[lXo+lMinY+lMaxZ])
                           {x+1,y+1}+((lXreal*lYreal*lZrM1)*l8is^[lXo+1+lMaxY+lMinZ])
                           {x+1,z+1}+((lXreal*lYrM1*lZreal)*l8is^[lXo+1+lMinY+lMaxZ])
                           {y+1,z+1}+((lXrM1*lYreal*lZreal)*l8is^[lXo+lMaxY+lMaxZ])
                           {x+1,y+1,z+1}+((lXreal*lYreal*lZreal)*l8is^[lXo+1+lMaxY+lMaxZ]) );
                          2:l16i^[lPos] :=
                           round ({all min} ( (lXrM1*lYrM1*lZrM1)*l16is^[lXo+lMinY+lMinZ])
                           {x+1}+((lXreal*lYrM1*lZrM1)*l16is^[lXo+1+lMinY+lMinZ])
                           {y+1}+((lXrM1*lYreal*lZrM1)*l16is^[lXo+lMaxY+lMinZ])
                           {z+1}+((lXrM1*lYrM1*lZreal)*l16is^[lXo+lMinY+lMaxZ])
                           {x+1,y+1}+((lXreal*lYreal*lZrM1)*l16is^[lXo+1+lMaxY+lMinZ])
                           {x+1,z+1}+((lXreal*lYrM1*lZreal)*l16is^[lXo+1+lMinY+lMaxZ])
                           {y+1,z+1}+((lXrM1*lYreal*lZreal)*l16is^[lXo+lMaxY+lMaxZ])
                           {x+1,y+1,z+1}+((lXreal*lYreal*lZreal)*l16is^[lXo+1+lMaxY+lMaxZ]) );
	                        4: l32f^[lPos] :=
                           {all min} ( (lXrM1*lYrM1*lZrM1)*l32fs^[lXo+lMinY+lMinZ])
                           {x+1}+((lXreal*lYrM1*lZrM1)*l32fs^[lXo+1+lMinY+lMinZ])
                           {y+1}+((lXrM1*lYreal*lZrM1)*l32fs^[lXo+lMaxY+lMinZ])
                           {z+1}+((lXrM1*lYrM1*lZreal)*l32fs^[lXo+lMinY+lMaxZ])
                           {x+1,y+1}+((lXreal*lYreal*lZrM1)*l32fs^[lXo+1+lMaxY+lMinZ])
                           {x+1,z+1}+((lXreal*lYrM1*lZreal)*l32fs^[lXo+1+lMinY+lMaxZ])
                           {y+1,z+1}+((lXrM1*lYreal*lZreal)*l32fs^[lXo+lMaxY+lMaxZ])
                           {x+1,y+1,z+1}+((lXreal*lYreal*lZreal)*l32fs^[lXo+1+lMaxY+lMaxZ]) ;
                      end; //case
                    end; //not perfectly in plane
                 end; //if voxel is in source image's bounding box
             end;//z
         end;//y
     end;//z
end else begin //if trilinear, else nearest neighbor
//nearest neighbor - added 12 April 2009
     for lZi := 0 to (lZ-1) do begin
         //these values are the same for all voxels in the slice
         // compute once per slice
         lZx := lZi*lMat.matrix[1][3];
         lZy := lZi*lMat.matrix[2][3];
         lZz := lZi*lMat.matrix[3][3];
         for lYi := 0 to (lY-1) do begin
             //these values change once per row
             // compute once per row
             lYx :=  lYi*lMat.matrix[1][2];
             lYy :=  lYi*lMat.matrix[2][2];
             lYz :=  lYi*lMat.matrix[3][2];
             for lXi := 0 to (lX-1) do begin
                 //compute each column
                 inc(lPos);
                 lXo := round(lXx^[lXi]+lYx+lZx+lMat.matrix[1][4]);
                 lYo := round(lXy^[lXi]+lYy+lZy+lMat.matrix[2][4]);
                 lZo := round(lXz^[lXi]+lYz+lZz+lMat.matrix[3][4]);
                 //need to test Xreal as -0.01 truncates to zero
                 if (lXo >= 0) and (lYo >= 0{1}) and (lZo >= 0{1}) and
                  (lXo < (lXs)) and (lYo < (lYs) ) and (lZo < (lZs))
                  //2012 removed -1 for nearest neighbor (lXo < (lXs -1)) and (lYo < (lYs -1) ) and (lZo < (lZs))
                  then begin
                    lOverlap := true;
                    inc(lXo);//images incremented from 1 not 0
		    lYo := lYo*lXs;
		    lZo := lZo*lXYs;
                    case lSrcHdr.ImgBufferBPP  of
                         1 : l8i^[lPos] :=l8is^[lXo+lYo+lZo];
 	                 2:  l16i^[lPos] :=l16is^[lXo+lYo+lZo];
	                 4: l32f^[lPos] :=l32fs^[lXo+lYo+lZo] ;
                    end; //case
                 end; //if voxel is in source image's bounding box
             end;//z
         end;//y
     end;//z
//end nearest neighbor
end;

     //release lookup tables
     freemem(lXx);
     freemem(lXy);
     freemem(lXz);
     //check to see if image is empty...
    if not lOverlap then
         Showmessage('No overlap between overlay and background - these images do not appear coregistered.');

    if lSrcHdr.ImgBufferBPP = 4 then begin
	    FreeMem(lSrcHdr.ImgBufferUnaligned);
	    GetMem(lSrcHdr.ImgBufferUnaligned ,(lOutVolItems*sizeof(Single)) + 16);
   {$IFDEF FPC}
	    lSrcHdr.ImgBuffer := align(lSrcHdr.ImgBufferUnaligned,16);
   {$ELSE}
	    lSrcHdr.ImgBuffer := ByteP($fffffff0 and (integer(lSrcHdr.ImgBufferUnaligned)+15));
   {$ENDIF}
	    lSrcHdr.ImgBufferItems := lOutVolItems;
	    move(l32f^,lSrcHdr.ImgBuffer^,(lOutVolItems*sizeof(Single)));
      FreeMem(lBuffOutUnaligned);
    end else if lSrcHdr.ImgBufferBPP = 2 then begin
	    FreeMem(lSrcHdr.ImgBufferUnaligned);
	    GetMem(lSrcHdr.ImgBufferUnaligned ,(lOutVolItems*sizeof(SmallInt)) + 16);
   {$IFDEF FPC}
	    lSrcHdr.ImgBuffer := align(lSrcHdr.ImgBufferUnaligned,16);
   {$ELSE}
	    lSrcHdr.ImgBuffer := ByteP($fffffff0 and (integer(lSrcHdr.ImgBufferUnaligned)+15));
   {$ENDIF}

	    lSrcHdr.ImgBufferItems := lOutVolItems;
	    //CopyMemory(Pointer(lSrcHdr.ImgBuffer),Pointer(l16i),(lOutVolItems*sizeof(SmallInt)));
	    move(l16i^,lSrcHdr.ImgBuffer^,(lOutVolItems*sizeof(SmallInt)));
	    FreeMem(lBuffOutUnaligned);
    end else if lSrcHdr.ImgBufferBPP = 1 then begin
	    FreeMem(lSrcHdr.ImgBufferUnaligned);
	    GetMem(lSrcHdr.ImgBufferUnaligned ,lOutVolItems + 16);
   {$IFDEF FPC}
	    lSrcHdr.ImgBuffer := align(lSrcHdr.ImgBufferUnaligned,16);
   {$ELSE}
	    lSrcHdr.ImgBuffer := ByteP($fffffff0 and (integer(lSrcHdr.ImgBufferUnaligned)+15));
   {$ENDIF}
	    lSrcHdr.ImgBufferItems := lOutVolItems;
	    //CopyMemory(Pointer(lSrcHdr.ImgBuffer),Pointer(l8i),lOutVolItems);
	    move(l8i^,lSrcHdr.ImgBuffer^,lOutVolItems);
	    FreeMem(l8i);
    end;
    lSrcHdr.NIfTIHdr := lDestHdr; //header inherits coordinates of target
end;


function Hdr2InvMat (lHdr: TNiftiHdr; var lOK: boolean): TMatrix;
var
  lSrcMat,lSrcMatInv: TMatrix;
begin
  lSrcMat := Hdr2Mat( lHdr);
  lSrcMatInv := lSrcMat;
  lOK := gaussj(lSrcMatInv);
  //the vectors should be rows not columns....
  //therefore we transpose the matrix
  //use this if you use transform instead of coord
  //Transposemat(lSrcMatInv);
  result := lSrcMatInv;
end;

end.