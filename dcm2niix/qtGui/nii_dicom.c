//#define MY_DEBUG

#include "nifti1.h"
#include "nii_dicom.h"
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> //toupper
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <float.h>
#include <stdint.h>
#include "nifti1_io_core.h"
#ifdef myUseCOut
#include <iostream>
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifdef MY_DEBUG
float deFuzz(float v) {
    if (fabs(v) < 0.00001)
        return 0;
    else
        return v;
    
}

void reportMat33(char *str, mat33 A) {
    printf("%s = [%g %g %g ; %g %g %g; %g %g %g ]\n",str,
           deFuzz(A.m[0][0]),deFuzz(A.m[0][1]),deFuzz(A.m[0][2]),
           deFuzz(A.m[1][0]),deFuzz(A.m[1][1]),deFuzz(A.m[1][2]),
           deFuzz(A.m[2][0]),deFuzz(A.m[2][1]),deFuzz(A.m[2][2]));
}

void reportMat44(char *str, mat44 A) {
    printf("%s = [%g %g %g %g; %g %g %g %g; %g %g %g %g; 0 0 0 1]\n",str,
           deFuzz(A.m[0][0]),deFuzz(A.m[0][1]),deFuzz(A.m[0][2]),deFuzz(A.m[0][3]),
           deFuzz(A.m[1][0]),deFuzz(A.m[1][1]),deFuzz(A.m[1][2]),deFuzz(A.m[1][3]),
           deFuzz(A.m[2][0]),deFuzz(A.m[2][1]),deFuzz(A.m[2][2]),deFuzz(A.m[2][3]));
}
#endif


int verify_slice_dir (struct TDICOMdata d, struct TDICOMdata d2, struct nifti_1_header *h, mat44 *R){
    //returns slice direction: 1=sag,2=coronal,3=axial, -= flipped
    if (h->dim[3] < 2) return 0; //don't care direction for single slice
    int iSL = 1; //find Z-slice direction: row with highest magnitude of 3rd column
    if ( (fabs(R->m[1][2]) >= fabs(R->m[0][2]))
        && (fabs(R->m[1][2]) >= fabs(R->m[2][2]))) iSL = 2; //
    if ( (fabs(R->m[2][2]) >= fabs(R->m[0][2]))
        && (fabs(R->m[2][2]) >= fabs(R->m[1][2]))) iSL = 3; //axial acquisition
    float pos = NAN;
    if ( !isnan(d2.patientPosition[iSL]) ) { //patient position fields exist
        pos = d2.patientPosition[iSL];
        if (isSameFloat(pos, d.patientPosition[iSL])) pos = NAN;
#ifdef MY_DEBUG
        if (!isnan(pos)) printf("position determined using lastFile %f\n",pos);
#endif
    }
    if (isnan(pos) &&( !isnan(d.patientPositionLast[iSL]) ) ) { //patient position fields exist
        pos = d.patientPositionLast[iSL];
        if (isSameFloat(pos, d.patientPosition[iSL])) pos = NAN;
#ifdef MY_DEBUG
        if (!isnan(pos)) printf("position determined using last (4d) %f\n",pos);
#endif
    }
    if (isnan(pos) && ( !isnan(d.stackOffcentre[iSL])) )
        pos = d.stackOffcentre[iSL];
    if (isnan(pos) && ( !isnan(d.lastScanLoc)) )
        pos = d.lastScanLoc;
    vec4 x;
    x.v[0] = 0; x.v[1] = 0; x.v[2]=h->dim[3]-1; x.v[3] = 1;
    vec4 pos1v = nifti_vect44mat44_mul(x, *R);
    float pos1 = pos1v.v[iSL-1];//-1 as C indexed from 0
    bool flip = false;
    if (!isnan(pos)) // we have real SliceLocation for last slice or volume center
        flip = (pos > R->m[iSL-1][3]) != (pos1 > R->m[iSL-1][3]); // same direction?, note C indices from 0
    else {// we do some guess work and warn user
#ifdef myUseCOut
     	std::cout<<"WARNING: Unable to determine slice direction: please check whether slices are flipped" <<std::endl;
#else
        printf("WARNING: Unable to determine slice direction: please check whether slices are flipped\n");
#endif
    }
    if (flip) {
        for (int i = 0; i < 4; i++)
            R->m[i][2] = -R->m[i][2];
    }
#ifdef MY_DEBUG
    printf("verify slice dir %d %d %d\n",h->dim[1],h->dim[2],h->dim[3]);
    reportMat44("Rout",*R);
    printf("iSL = %d\n",iSL);
    printf(" pos1 = %f\n",pos1);
#endif
    if (flip)
        return -iSL;
    else
        return iSL;
} //verify_slice_dir()

void setQSForm(struct nifti_1_header *h, mat44 Q44) {
    h->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    h->srow_x[0] = Q44.m[0][0];
    h->srow_x[1] = Q44.m[0][1];
    h->srow_x[2] = Q44.m[0][2];
    h->srow_x[3] = Q44.m[0][3];
    h->srow_y[0] = Q44.m[1][0];
    h->srow_y[1] = Q44.m[1][1];
    h->srow_y[2] = Q44.m[1][2];
    h->srow_y[3] = Q44.m[1][3];
    h->srow_z[0] = Q44.m[2][0];
    h->srow_z[1] = Q44.m[2][1];
    h->srow_z[2] = Q44.m[2][2];
    h->srow_z[3] = Q44.m[2][3];
    float dumdx, dumdy, dumdz;
    nifti_mat44_to_quatern( Q44 , &h->quatern_b, &h->quatern_c, &h->quatern_d,&h->qoffset_x, &h->qoffset_y, &h->qoffset_z, &dumdx, &dumdy, &dumdz,&h->pixdim[0]) ;
    h->qform_code = NIFTI_XFORM_SCANNER_ANAT;
} //setQSForm()



int headerDcm2NiiSForm(struct TDICOMdata d, struct TDICOMdata d2,  struct nifti_1_header *h) { //fill header s and q form
    //see http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1_io.c
    //returns sliceDirection: 0=unknown,1=sag,2=coro,3=axial,-=reversed slices
    //
    int sliceDir = 0;
    if (h->dim[3] < 2) return sliceDir; //don't care direction for single slice
    h->sform_code = NIFTI_XFORM_UNKNOWN;
    h->qform_code = NIFTI_XFORM_UNKNOWN;
    bool isOK = false;
    for (int i = 1; i <= 6; i++)
        if (d.orient[i] != 0.0) isOK = true;
    if (!isOK) return sliceDir;
    mat44 Q44 = nifti_dicom2mat(d.orient, d.patientPosition, d.xyzMM);
    if (d.CSA.mosaicSlices > 1) {
        double nRowCol = ceil(sqrt(d.CSA.mosaicSlices));
        double lFactorX = (d.xyzDim[1] -(d.xyzDim[1]/nRowCol)   )/2.0;
        double lFactorY = (d.xyzDim[2] -(d.xyzDim[2]/nRowCol)   )/2.0;
        Q44.m[0][3] =(Q44.m[0][0]*lFactorX)+(Q44.m[0][1]*lFactorY)+Q44.m[0][3];
        Q44.m[1][3] =(Q44.m[1][0]*lFactorX)+(Q44.m[1][1]*lFactorY)+Q44.m[1][3];
        Q44.m[2][3] =(Q44.m[2][0]*lFactorX)+(Q44.m[2][1]*lFactorY)+Q44.m[2][3];
        
        /* #ifdef obsolete_mosaic_flip
         double val = d.xyzDim[2]/nRowCol; //obsolete!!!
         //Q44 now equals 'dicom_to_patient' in spm_dicom_convert
         mat44 patient_to_tal, analyze_to_dicom;
         LOAD_MAT44(patient_to_tal, -1.0l,0.0l,0.0l,0.0l, 0.0l,-1.0l,0.0l,0.0l, 0.0l,0.0l,1.0l,0.0l);
         LOAD_MAT44(analyze_to_dicom, 1.0l,0.0l,0.0l,-1.0l, 0.0l,-1.0l,0.0l,val, 0.0l,0.0l,1.0l,-1.0l);
         Q44 = nifti_mat44_mul(patient_to_tal,Q44);
         Q44 = nifti_mat44_mul(Q44,analyze_to_dicom);
         //Q44 now equals 'mat' in spm_dicom_convert
         //subasgn.m in SPM5 translates by one voxel...
         LOAD_MAT44(analyze_to_dicom, 1.0l,0.0l,0.0l,1.0l, 0.0l,1.0l,0.0l,1.0l, 0.0l,0.0l,1.0l,1.0l);
         Q44 = nifti_mat44_mul(Q44,analyze_to_dicom);
         #else */
        for (int c=0; c<2; c++)
            for (int r=0; r<4; r++)
                Q44.m[c][r] = -Q44.m[c][r];
        // #endif
        mat33 Q;
        LOAD_MAT33(Q, d.orient[1], d.orient[4],d.CSA.sliceNormV[1],
                   d.orient[2],d.orient[5],d.CSA.sliceNormV[2],
                   d.orient[3],d.orient[6],d.CSA.sliceNormV[3]);
        if  (nifti_mat33_determ(Q) < 0) { //Siemens sagittal are R>>L, whereas NIfTI is L>>R, we retain Siemens order on disk so ascending is still ascending, but we need to have the spatial transform reflect this.
            mat44 det;
            sliceDir = kSliceOrientMosaicNegativeDeterminant; //we need to handle DTI vectors accordingly
            LOAD_MAT44(det, 1.0l,0.0l,0.0l,0.0l, 0.0l,1.0l,0.0l,0.0l, 0.0l,0.0l,-1.0l,0.0l);
            //patient_to_tal.m[2][3] = 1-d.CSA.MosaicSlices;
            Q44 = nifti_mat44_mul(Q44,det);
        }
    } else { //not a mosaic
        sliceDir = verify_slice_dir(d, d2, h, &Q44);
        for (int c=0; c<4; c++)// LPS to nifti RAS, xform matrix before reorient
            for (int r=0; r<2; r++) //swap rows 1 & 2
                Q44.m[r][c] = - Q44.m[r][c];
#ifdef MY_DEBUG
        reportMat44("final",Q44);
#endif
    }
    setQSForm(h,Q44);
    return sliceDir;
} //headerDcm2NiiSForm()

int headerDcm2Nii2(struct TDICOMdata d, struct TDICOMdata d2, struct nifti_1_header *h) { //final pass after de-mosaic
    char txt[1024] = {""};
    if (h->slice_code == NIFTI_SLICE_UNKNOWN) h->slice_code = d.CSA.sliceOrder;
    if (h->slice_code == NIFTI_SLICE_UNKNOWN) h->slice_code = d2.CSA.sliceOrder; //sometimes the first slice order is screwed up https://github.com/eauerbach/CMRR-MB/issues/29
    sprintf(txt, "TE=%.2g;Time=%.3f", d.TE,d.acquisitionTime);// d.dateTime);
    if (d.CSA.phaseEncodingDirectionPositive > 0) {
        char dtxt[1024] = {""};
        sprintf(dtxt, ";phase+=%d", d.CSA.phaseEncodingDirectionPositive);
        strcat(txt,dtxt);
    }
    if ((d.CSA.bandwidthPerPixelPhaseEncode > 0) && ((d.phaseEncodingRC =='C') || (d.phaseEncodingRC =='R'))) {
        float dwellTime = 0;
        if (d.phaseEncodingRC =='C')
            dwellTime =  1000/d.CSA.bandwidthPerPixelPhaseEncode/h->dim[2];
        else
            dwellTime =  1000/d.CSA.bandwidthPerPixelPhaseEncode/h->dim[1];
        char dtxt[1024] = {""};
        sprintf(dtxt, ";dwell=%.3f", dwellTime);
        strcat(txt,dtxt);
    }
    //printf(" Description= '%s'  [length=%lu]\n",txt, strlen(txt));
    //x strlcpy(h->descrip,txt,80);
    snprintf(h->descrip,80, "%s",txt);
    if (strlen(d.imageComments) > 0)
        snprintf(h->aux_file,24,"%s",d.imageComments);
    return headerDcm2NiiSForm(d,d2, h);
} //headerDcm2Nii2()

int dcmStrLen (int len) {
    if (len < kDICOMStr)
        return len+1;
    else
        return kDICOMStr;
} //dcmStrLen()

struct TDICOMdata clear_dicom_data() {
    struct TDICOMdata d;
    
    d.locationsInAcquisition = 0;
    for (int i=0; i < 4; i++) {
        for (int n=0; n < kMaxDTIv; n++)
            d.CSA.dtiV[n][i] = 0;
        d.patientPosition[i] = NAN;
        //d.patientPosition2nd[i] = NAN; //used to distinguish XYZT vs XYTZ for Philips 4D
        d.patientPositionLast[i] = NAN; //used to compute slice direction for Philips 4D
        d.stackOffcentre[i] = NAN;
        d.angulation[i] = 0.0f;
        d.xyzMM[i] = 1;
    }
    d.CSA.numDti = 0;
    for (int i=0; i < 5; i++) d.xyzDim[i] = 1;
        for (int i = 0; i < 7; i++) d.orient[i] = 0.0f;
            d.patientPositionSequentialRepeats = 0;//d.isHasMixed = false;
            d.isHasPhase = false;
            d.isHasMagnitude = false;
            d.sliceOrient = kSliceOrientUnknown;
            
            //strcpy(d.sliceOrient,"\n");
            strcpy(d.patientName, "John_Doe");
            strcpy(d.patientID, "ID123");
            strcpy(d.imageComments, "imgComments");
            strcpy(d.studyDate, "1/1/1977");
            strcpy(d.studyTime, "11:11:11");
            d.dateTime = (double)19770703150928.0;
            d.acquisitionTime = 0.0f;
            strcpy(d.protocolName, "MPRAGE");
            d.manufacturer = kMANUFACTURER_UNKNOWN;
            d.isPlanarRGB = false;
            d.lastScanLoc = NAN;
            d.TR = 0;
            d.TE = 0;
            //d.locationsInAcquisition = 0;
            d.numberOfDynamicScans = 0;
            d.imageNum = 0;
            d.intenScale = 1;
            d.intenIntercept = 0;
            d.seriesNum = 1;
            d.acquNum = 0;
            d.imageNum = 1;
            d.imageStart = 0;
            d.is3DAcq = false; //e.g. MP-RAGE, SPACE, TFE
            d.bitsAllocated = 16;//bits
            d.bitsStored = 0;
            d.samplesPerPixel = 1;
            d.isValid = false;
            d.isSigned = false; //default is unsigned!
            d.isExplicitVR = true;
            d.isLittleEndian = true; //DICOM initially always little endian
            d.converted2NII = 0;
            d.phaseEncodingRC = '?';
            d.CSA.bandwidthPerPixelPhaseEncode = 0.0;
            d.CSA.mosaicSlices = 0;
            d.CSA.sliceOrder = NIFTI_SLICE_UNKNOWN;
            d.CSA.protocolSliceNumber1 = 0;
            d.CSA.phaseEncodingDirectionPositive = -1; //unknown
            return d;
} //clear_dicom_data()

void dcmStrDigitsOnly(char* lStr) {
    //e.g. change "H11" to " 11"
    size_t len = strlen(lStr);
    if (len < 1) return;
    for (int i = 0; i < len; i++)
        if (!isdigit(lStr[i]) )
            lStr[i] = ' ';
    
}


void dcmStr(int lLength, unsigned char lBuffer[], char* lOut) {
    //char test[] = " 1     2    3    ";
    //lLength = (int)strlen(test);
    if (lLength < 1) return;
    char cString[lLength+1];
    cString[lLength] =0;
    memcpy(cString, (char*)&lBuffer[0], lLength);
    //memcpy(cString, test, lLength);
    for (int i = 0; i < lLength; i++)
        if (cString[i]==' ') cString[i] = '_';
    int len = 1;
    for (int i = 1; i < lLength; i++) { //remove repeated "_"
        if ((cString[i-1]!='_') || (cString[i]!='_')) {
            cString[len] =cString[i];
            len++;
        }
    } //for each item
    if (cString[len-1] == '_') len--;
    //while ((len > 0) && (cString[len]=='_')) len--; //remove trailing '_'
    cString[len] = 0; //null-terminate, strlcpy does this anyway
    len = dcmStrLen(len);
    //strlcpy(lOut,cString,len);
    memcpy(lOut,cString,len-1);
    lOut[len-1] = 0;
    //printf(">>%s<<>>%s<<\n",cString,lOut);
    //lOut[len-1]='z';
    //printf(">>%s<<>>%s<<\n",cString,lOut);
    //strlcpy(lOut,cString,len);
} //dcmStr()

float dcmFloat(int lByteLength, unsigned char lBuffer[], bool littleEndian) {//read binary 32-bit float
    //http://stackoverflow.com/questions/2782725/converting-float-values-from-big-endian-to-little-endian
#ifdef __BIG_ENDIAN__
    bool swap = littleEndian;
#else
    bool swap = !littleEndian;
#endif
    float retVal;
    memcpy(&retVal, (char*)&lBuffer[0], 4);
    if (!swap) return retVal;
    char *floatToConvert = ( char* ) & lBuffer;
    char *returnFloat = ( char* ) & retVal;
    //swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];
    //printf("swapped val = %f\n",retVal);
    return retVal;
} //dcmFloat()

int dcmInt (int lByteLength, unsigned char lBuffer[], bool littleEndian) { //read binary 16 or 32 bit integer
    if (littleEndian) {
        if (lByteLength <= 3)
            return  lBuffer[0] | (lBuffer[1]<<8); //shortint vs word?
        return lBuffer[0]+(lBuffer[1]<<8)+(lBuffer[2]<<16)+(lBuffer[3]<<24); //shortint vs word?
    }
    if (lByteLength <= 3)
        return  lBuffer[1] | (lBuffer[0]<<8); //shortint vs word?
    return lBuffer[3]+(lBuffer[2]<<8)+(lBuffer[1]<<16)+(lBuffer[0]<<24); //shortint vs word?
} //dcmInt()

int dcmStrInt (int lByteLength, unsigned char lBuffer[]) {//read float stored as a string
    char cString[lByteLength+1];
    cString[lByteLength] =0;
    memcpy(cString, (char*)&lBuffer[0], lByteLength);
    //printf(" --> *%s* %s%s\n",cString, &lBuffer[0],&lBuffer[1]);
    return atoi(cString);
} //dcmStrInt()

int dcmStrManufacturer (int lByteLength, unsigned char lBuffer[]) {//read float stored as a string
    if (lByteLength < 2) return kMANUFACTURER_UNKNOWN;
    char cString[lByteLength+1];
    cString[lByteLength] =0;
    memcpy(cString, (char*)&lBuffer[0], lByteLength);
    //printf("MANU %s\n",cString);
    if ((toupper(cString[0])== 'S') && (toupper(cString[1])== 'I'))
        return kMANUFACTURER_SIEMENS;
    if ((toupper(cString[0])== 'G') && (toupper(cString[1])== 'E'))
        return kMANUFACTURER_GE;
    if ((toupper(cString[0])== 'P') && (toupper(cString[1])== 'H'))
        return kMANUFACTURER_PHILIPS;
    return kMANUFACTURER_UNKNOWN;
} //dcmStrManufacturer

typedef struct __attribute__((packed)) {
    char name[64]; //null-terminated
    int32_t vm;
    char vr[4]; //  possibly nul-term string
    int32_t syngodt;//  ??
    int32_t nitems;// number of items in CSA
    int32_t xx;// maybe == 77 or 205
} TCSAtag; //Siemens csa tag structure
typedef struct __attribute__((packed)) {
    int32_t xx1, xx2_Len, xx3_77, xx4;
} TCSAitem; //Siemens csa item structure

float csaMultiFloat (unsigned char buff[], int nItems, float Floats[], int *ItemsOK) {
    //warning: lFloats indexed from 1! will fill lFloats[1]..[nFloats]
    //if lnItems == 1, returns first item, if lnItems > 1 returns index of final successful conversion
    TCSAitem itemCSA;
    *ItemsOK = 0;
    if (nItems < 1)  return 0.0f;
    int lPos = 0;
    for (int lI = 1; lI <= nItems; lI++) {
        memcpy(&itemCSA, &buff[lPos], sizeof(itemCSA));
        lPos +=sizeof(itemCSA);
        if (itemCSA.xx2_Len > 0) {
            char cString[itemCSA.xx2_Len];
            memcpy(&cString, &buff[lPos], sizeof(cString));
            lPos += ((itemCSA.xx2_Len +3)/4)*4;
            //printf(" %d item length %d = %s\n",lI, itemCSA.xx2_Len, cString);
            Floats[lI] = atof(cString);
            *ItemsOK = lI; //some sequences have store empty items
        }
    } //for each item
    return Floats[1];
} //csaMultiFloat()

int readCSAImageHeader(unsigned char *buff, int lLength, struct TCSAdata *CSA, bool isVerbose) {
    //see also http://afni.nimh.nih.gov/pub/dist/src/siemens_dicom_csa.c
    //printf("%c%c%c%c\n",buff[0],buff[1],buff[2],buff[3]);
    if (lLength < 36) return EXIT_FAILURE;
    if ((buff[0] != 'S') || (buff[1] != 'V') || (buff[2] != '1') || (buff[3] != '0') ) return EXIT_FAILURE;
    int lPos = 8; //skip 8 bytes of data, 'SV10' plus  2 32-bit values unused1 and unused2
    int lnTag = buff[lPos]+(buff[lPos+1]<<8)+(buff[lPos+2]<<16)+(buff[lPos+3]<<24);
    if (buff[lPos+4] != 77) return EXIT_FAILURE;
    lPos += 8; //skip 8 bytes of data, 32-bit lnTag plus 77 00 00 0
    TCSAtag tagCSA;
    TCSAitem itemCSA;
    int itemsOK;
    float lFloats[7];
    for (int lT = 1; lT <= lnTag; lT++) {
        memcpy(&tagCSA, &buff[lPos], sizeof(tagCSA)); //read tag
        lPos +=sizeof(tagCSA);
        //printf("%d CSA of %s %d\n",lPos, tagCSA.name, tagCSA.nitems);
        if (tagCSA.nitems > 0) {
            if (strcmp(tagCSA.name, "NumberOfImagesInMosaic") == 0)
                CSA->mosaicSlices = round(csaMultiFloat (&buff[lPos], 1,lFloats, &itemsOK));
            else if (strcmp(tagCSA.name, "B_value") == 0) {
                CSA->dtiV[0][0] = csaMultiFloat (&buff[lPos], 1,lFloats, &itemsOK);
                CSA->numDti = 1; //triggered by b-value, as B0 images do not have DiffusionGradientDirection tag
            }
            else if ((strcmp(tagCSA.name, "DiffusionGradientDirection") == 0) && (tagCSA.nitems > 2)){
                CSA->dtiV[0][1] = csaMultiFloat (&buff[lPos], 3,lFloats, &itemsOK);
                CSA->dtiV[0][2] = lFloats[2];
                CSA->dtiV[0][3] = lFloats[3];
                if (isVerbose)
                    printf("DiffusionGradientDirection %f %f %f\n",lFloats[1],lFloats[2],lFloats[3]);
            } else if ((strcmp(tagCSA.name, "SliceNormalVector") == 0) && (tagCSA.nitems > 2)){
                CSA->sliceNormV[1] = csaMultiFloat (&buff[lPos], 3,lFloats, &itemsOK);
                CSA->sliceNormV[2] = lFloats[2];
                CSA->sliceNormV[3] = lFloats[3];
                if (isVerbose)
                    printf("SliceNormalVector %f %f %f\n",CSA->sliceNormV[1],CSA->sliceNormV[2],CSA->sliceNormV[3]);
            } else if (strcmp(tagCSA.name, "SliceMeasurementDuration") == 0)
                CSA->sliceMeasurementDuration = csaMultiFloat (&buff[lPos], 3,lFloats, &itemsOK);
            else if (strcmp(tagCSA.name, "BandwidthPerPixelPhaseEncode") == 0)
                CSA->bandwidthPerPixelPhaseEncode = csaMultiFloat (&buff[lPos], 3,lFloats, &itemsOK);
            else if ((strcmp(tagCSA.name, "MosaicRefAcqTimes") == 0) && (tagCSA.nitems > 3)  ){
                float sliceTimes[tagCSA.nitems+1];
                csaMultiFloat (&buff[lPos], tagCSA.nitems,sliceTimes, &itemsOK);
                float minTimeIndex, minTimeValue, timeValue1;
                CSA->multiBandFactor = 1;
                timeValue1 = sliceTimes[1];
                minTimeIndex = 1;
                minTimeValue= sliceTimes[1];
                for (int z = 2; z <= itemsOK; z++) { //find index and value of fastest time
                    if (sliceTimes[z] < minTimeValue) {
                        minTimeValue = sliceTimes[z];
                        minTimeIndex = z;
                    }
                    if (sliceTimes[z] == timeValue1) CSA->multiBandFactor++;
                }
                if (minTimeIndex == 2)
                    CSA->sliceOrder = NIFTI_SLICE_ALT_INC2;// e.g. 3,1,4,2
                else if (minTimeIndex == (itemsOK-1))
                    CSA->sliceOrder = NIFTI_SLICE_ALT_DEC2;// e.g. 4,3,2,1
                else if ((minTimeIndex == 1) && (sliceTimes[2] < sliceTimes[3]))
                    CSA->sliceOrder = NIFTI_SLICE_SEQ_INC;
                else if ((minTimeIndex == 1) && (sliceTimes[2] > sliceTimes[3]))
                    CSA->sliceOrder = NIFTI_SLICE_ALT_INC;
                else if ((minTimeIndex == itemsOK) && (sliceTimes[itemsOK-1] < sliceTimes[itemsOK]))
                    CSA->sliceOrder = NIFTI_SLICE_SEQ_DEC;
                else if ((minTimeIndex == itemsOK) && (sliceTimes[itemsOK-1] > sliceTimes[itemsOK-2]))
                    CSA->sliceOrder = NIFTI_SLICE_ALT_DEC;
                else {
                    /*NSMutableArray *sliceTimesNS = [NSMutableArray arrayWithCapacity:tagCSA.nitems];
                     for (int z = 1; z <= itemsOK; z++)
                     [sliceTimesNS addObject:[NSNumber numberWithFloat:sliceTimes[z]]];
                     NSLog(@" Warning: unable to determine slice order for %lu slice mosaic: %@",(unsigned long)[sliceTimesNS count],sliceTimesNS );
                     */
                    printf("Warning: unable to determine slice order from CSA tag MosaicRefAcqTimes\n");
                }
            } else if (strcmp(tagCSA.name, "ProtocolSliceNumber") == 0)
                CSA->protocolSliceNumber1 = round (csaMultiFloat (&buff[lPos], 1,lFloats, &itemsOK));
            else if (strcmp(tagCSA.name, "PhaseEncodingDirectionPositive") == 0)
                CSA->phaseEncodingDirectionPositive = round (csaMultiFloat (&buff[lPos], 1,lFloats, &itemsOK));
            /*if (strcmp(tagCSA.name, "SlicePosition_PCS") == 0) {
             float f =  (csaMultiFloat (&buff[lPos], 1,lFloats, &itemsOK));
             for (int k = 1; k <= tagCSA.nitems; k++)
             printf("PCS %d = %f\n",k, lFloats[k]);
             }
             else {
             printf("unused CSA tag %s with %d items\n",tagCSA.name, tagCSA.nitems);
             }*/
            for (int lI = 1; lI <= tagCSA.nitems; lI++) {
                memcpy(&itemCSA, &buff[lPos], sizeof(itemCSA));
                lPos +=sizeof(itemCSA);
                lPos += ((itemCSA.xx2_Len +3)/4)*4;
            }
        } //if at least 1 item
    }// for lT 1..lnTag
    return EXIT_SUCCESS;
} // readCSAImageHeader()

void dcmMultiFloat (int lByteLength, char lBuffer[], int lnFloats, float *lFloats) {
    //warning: lFloats indexed from 1! will fill lFloats[1]..[nFloats]
    if ((lnFloats < 1) || (lByteLength < 1)) return;
    char cString[lByteLength+1];
    memcpy(cString, (char*)&lBuffer[0], lByteLength);
    cString[lByteLength] = 0; //null terminate
    char *temp=( char *)malloc(lByteLength+1);
    int f = 0,lStart = 0;
    bool isOK = false;
    for (int i = 0; i <= lByteLength; i++) {
        if ((lBuffer[i] >= '0') && (lBuffer[i] <= '9')) isOK = true;
        if ((isOK) && ((i == (lByteLength)) || (lBuffer[i] == '/')  || (lBuffer[i] == ' ')  || (lBuffer[i] == '\\') )){
            //x strlcpy(temp,&cString[lStart],i-lStart+1);
            snprintf(temp,i-lStart+1,"%s",&cString[lStart]);
            //printf("dcmMultiFloat %s\n",temp);
            if (f < lnFloats) {
                f ++;
                lFloats[f] = atof(temp);
                //printf("%d == %f\n", f, atof(temp));
            } //if f <= nFloats
            lStart = i+1;
        } //if isOK
    }  //for i to length
    free(temp);
} //dcmMultiFloat()

float dcmStrFloat (int lByteLength, unsigned char lBuffer[]) { //read float stored as a string
    char cString[lByteLength+1];
    memcpy(cString, (char*)&lBuffer[0], lByteLength);
    cString[lByteLength] = 0; //null terminate
    return atof(cString);
} //dcmStrFloat()

/*
 void dcmStrReport (int lByteLength, char lBuffer[]) {//print value at location
 char cString[lByteLength+1];
 memcpy(cString, (char*)&lBuffer[0], lByteLength);
 cString[lByteLength] = 0; //null terminate
 printf("%d dcmStrReport '%s'\n",lByteLength, cString);
 } //dcmStrReport
 */

int headerDcm2Nii(struct TDICOMdata d, struct nifti_1_header *h) {
    //printf("bytes %dx%dx%d %d, %d\n",d.XYZdim[1],d.XYZdim[2],d.XYZdim[3], d.Allocbits_per_pixel, d.samplesPerPixel);
    for (int i = 0; i < 80; i++) h->descrip[i] = 0;
    for (int i = 0; i < 24; i++) h->aux_file[i] = 0;
    for (int i = 0; i < 18; i++) h->db_name[i] = 0;
    for (int i = 0; i < 10; i++) h->data_type[i] = 0;
    for (int i = 0; i < 16; i++) h->intent_name[i] = 0;
    if ((d.bitsAllocated == 8) && (d.samplesPerPixel == 3))
        h->datatype = DT_RGB24;
    else if ((d.bitsAllocated == 8) && (d.samplesPerPixel == 1))
        h->datatype = DT_UINT8;
    else if ((d.bitsAllocated == 16) && (d.samplesPerPixel == 1) && (d.isSigned))
        h->datatype = DT_INT16;
    else if ((d.bitsAllocated == 16) && (d.samplesPerPixel == 1) && (!d.isSigned))
        h->datatype = DT_UINT16;
    else if (d.bitsAllocated == 32)
        h->datatype = DT_INT32;
    else {
#ifdef myUseCOut
     	std::cout<<"Unsupported DICOM bit-depth " <<d.bitsAllocated << " with " << d.samplesPerPixel << "samples per pixel" <<std::endl;
#else
        printf("Unsupported DICOM bit-depth %d with %d samples per pixel\n",d.bitsAllocated,d.samplesPerPixel);
#endif
        
        return EXIT_FAILURE;
    }
    if ((h->datatype == DT_UINT16) && (d.bitsStored > 0) &&(d.bitsStored < 16))
        h->datatype = DT_INT16; // DT_INT16 is more widely supported, same represenation for values 0..32767
    for (int i = 0; i < 8; i++) {
        h->pixdim[i] = 0.0f;
        h->dim[i] = 0;
    }
    h->regular = 114;
    h->scl_inter = d.intenIntercept;
    h->scl_slope = d.intenScale;
    h->cal_max = 0;
    h->cal_min = 0;
    h->magic[0]='n';
    h->magic[1]='+';
    h->magic[2]='1';
    h->magic[3]='\0';
    h->vox_offset = d.imageStart;
    h->bitpix = d.bitsAllocated * d.samplesPerPixel;
    h->pixdim[1] = d.xyzMM[1];
    h->pixdim[2] = d.xyzMM[2];
    h->pixdim[3] = d.xyzMM[3];
    h->pixdim[4] = d.TR/1000; //TR reported in msec, time is in sec
    h->dim[1] = d.xyzDim[1];
    h->dim[2] = d.xyzDim[2];
    h->dim[3] = d.xyzDim[3];
    h->dim[4] = d.xyzDim[4];
    if (h->dim[4] < 2)
        h->dim[0] = 3;
    else
        h->dim[0] = 4;
    
    for (int i = 0; i <= 3; i++) {
        h->srow_x[i] = 0.0f;
        h->srow_y[i] = 0.0f;
        h->srow_z[i] = 0.0f;
    }
    h->srow_x[0] = -1;
    h->srow_y[2] = 1;
    h->srow_z[1] = -1;
    h->srow_x[3] =  (h->dim[1] /2);
    h->srow_y[3] = -(h->dim[3] /2);
    h->srow_z[3] =  (h->dim[2] /2);
    h->qform_code = NIFTI_XFORM_UNKNOWN;
    h->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    h->toffset = 0;
    h->intent_code = NIFTI_INTENT_NONE;
    h->dim_info = 0; //Freq, Phase and Slice all unknown
    h->xyzt_units = NIFTI_UNITS_UNKNOWN;
    h->slice_duration = 0; //avoid +inf/-inf, NaN
    h->intent_p1 = 0;  //avoid +inf/-inf, NaN
    h->intent_p2 = 0;  //avoid +inf/-inf, NaN
    h->intent_p3 = 0;  //avoid +inf/-inf, NaN
    h->pixdim[0] = 1; //QFactor should be 1 or -1
    h->sizeof_hdr = 348; //used to signify header does not need to be byte-swapped
    h->slice_code = d.CSA.sliceOrder;
    headerDcm2Nii2(d, d, h);
    return EXIT_SUCCESS;
} // headerDcm2Nii()

bool isFloatDiff (float a, float b) {
    return (fabs (a - b) > FLT_EPSILON);
} //isFloatDiff()

ivec3 setVec3i(int x, int y, int z)
{
    ivec3 v = {x, y, z};
    return v;
} //setVec3i()

mat33 nifti_mat33_reorder_cols( mat33 m, ivec3 v ) {
    // matlab equivalent ret = m(:, v); where v is 1,2,3 [INDEXED FROM ONE!!!!]
    mat33 ret;
    for (int r=0; r<3; r++) {
        for(int c=0; c<3; c++)
            ret.m[r][c] = m.m[r][v.v[c]-1];
    }
    return ret;
} //nifti_mat33_reorder_cols()

/*bool isExt (char *file_name, const char* ext) {
 char *p_extension;
 if((p_extension = strrchr(file_name,'.')) != NULL )
 if(strcmp(p_extension,ext) == 0) return true;
 return false;
 
 }*/

void changeExt (char *file_name, const char* ext) {
    char *p_extension;
    p_extension = strrchr(file_name, '.');
    if (p_extension)
    {
        strcpy(++p_extension, ext);
    }
} //changeExt()

struct TDICOMdata  nii_readParRec (char * parname) {
    struct TDICOMdata d = clear_dicom_data();
    FILE *fp = fopen(parname, "r");
    if (fp == NULL) return d;
#define LINESZ 2048
#define	kSlice	0
#define	kEcho	1
#define	kDyn	2
#define	kCardiac	3
#define	kImageType	4
#define	kSequence	5
#define	kIndex	6
#define	kBitsPerVoxel	7
#define	kXdim	9
#define	kYdim	10
#define	kRI	11
#define	kRS	12
#define	kSS	13
#define	kAngulationAPs	16 //In V4, offcentre and Angulation labeled as y z x, but actually x y z!
#define	kAngulationFHs	17
#define	kAngulationRLs	18
#define	kPositionAP	19
#define	kPositionFH	20
#define	kPositionRL	21
#define	kThickmm	22
#define	kGapmm	23
#define kSliceOrients 25
#define	kXmm	28
#define	kYmm	29
#define	kTEcho	30
#define	kDynTime	31
#define	kGradientNumber 42
#define	kbval 33
#define	kv1	47
#define	kv2	45
#define	kv3	46
#define	kASL	48
    char buff[LINESZ];
    bool ADCwarning = false;
    for (int n=0; kMaxDTIv < 4; n++)
        d.CSA.dtiV[n][0] = -1; //set to impossible value to detect re-usage by ADC map
    int parVers = 0;
    int nCols = 26;
    int slice = 0;
    //int prevSliceIndex = 0; //index of prior slice: detect if images are not in order
    const int kMaxCols = 49;
    float *cols = (float *)malloc(sizeof(float) * kMaxCols);
    char *p = fgets (buff, LINESZ, fp);
    bool isIntenScaleVaries = false;
    bool isIndexSequential = true;
    while (p) {
        if (strlen(buff) < 1)
            continue;
        if (buff[0] == '#') { //comment
            char Comment[7][50];
            sscanf(buff, "# %s %s\n", Comment[0], Comment[1]);
            if (strcmp(Comment[1], "TRYOUT") == 0) {
                sscanf(buff, "# %s %s %s %s %s %s V%s\n", Comment[0], Comment[1], Comment[2], Comment[3]
                       ,Comment[4], Comment[5],Comment[6]);
                parVers = (atof(Comment[6])*10); //4.2 = 42 etc
                if (parVers < 40) {
#ifdef myUseCOut
                    std::cout<<"This software is unable to convert ancient PAR files: please use legacy dcm2nii" <<std::endl;
#else
                    printf("This software is unable to convert ancient PAR files: please use legacy dcm2nii\n");
#endif
                    return d;
                    //nCols = 26; //e.g. PAR 3.0 has 26 relevant columns
                } else if (parVers < 41)
                    nCols = 32; //e.g PAR 4.0
                else if (parVers < 42)
                    nCols = 47; //e.g. PAR 4.1
                else
                    nCols = kMaxCols; //e.g. PAR 4.2
            }
            p = fgets (buff, LINESZ, fp);//get next line
            continue;
        } //process '#' comment
        if (buff[0] == '.') { //tag
            char Comment[8][50];
            sscanf(buff, ". %s %s %s %s %s %s %s %s\n", Comment[0], Comment[1],Comment[2], Comment[3], Comment[4], Comment[5], Comment[6], Comment[7]);
            
            if ((strcmp(Comment[0], "Acquisition") == 0) && (strcmp(Comment[1], "nr") == 0)) {
                d.acquNum = atoi( Comment[3]);
                d.seriesNum = d.acquNum;
            }
            if ((strcmp(Comment[0], "Repetition") == 0) && (strcmp(Comment[1], "time") == 0))
                d.TR = atof(Comment[4]);
            if ((strcmp(Comment[0], "Patient") == 0) && (strcmp(Comment[1], "name") == 0)) {
                strcpy(d.patientName, Comment[3]);
                strcat(d.patientName, Comment[4]);
                strcat(d.patientName, Comment[5]);
                strcat(d.patientName, Comment[6]);
                strcat(d.patientName, Comment[7]);
                //printf("%s\n",d.patientName);
                
            }
            if ((strcmp(Comment[0], "Protocol") == 0) && (strcmp(Comment[1], "name") == 0)) {
                strcpy(d.protocolName, Comment[3]);
                strcat(d.protocolName, Comment[4]);
                strcat(d.protocolName, Comment[5]);
                strcat(d.protocolName, Comment[6]);
                strcat(d.protocolName, Comment[7]);
                //printf("%s\n",d.protocolName);
            }
            if ((strcmp(Comment[0], "Examination") == 0) && (strcmp(Comment[1], "date/time") == 0)) {
                strcpy(d.studyDate, Comment[3]);
                strcpy(d.studyTime, Comment[5]);
                //to do convert to traditional DICOM style date time
            }
            if ((strcmp(Comment[0], "Off") == 0) && (strcmp(Comment[1], "Centre") == 0)) {
                //Off Centre midslice(ap,fh,rl) [mm]
                d.stackOffcentre[2] = atof(Comment[5]);
                d.stackOffcentre[3] = atof(Comment[6]);
                d.stackOffcentre[1] = atof(Comment[7]);
            }
            if ((strcmp(Comment[0], "Patient") == 0) && (strcmp(Comment[1], "position") == 0)) {
                //Off Centre midslice(ap,fh,rl) [mm]
                d.patientOrient[0] = toupper(Comment[3][0]);
                d.patientOrient[1] = toupper(Comment[4][0]);
                d.patientOrient[2] = toupper(Comment[5][0]);
                d.patientOrient[3] = 0;
            }
            if ((strcmp(Comment[0], "Max.") == 0) && (strcmp(Comment[3], "slices/locations") == 0)) {
                d.xyzDim[3] = atoi(Comment[5]);
            }
            p = fgets (buff, LINESZ, fp);//get next line
            continue;
        } //process '.' tag
        if (strlen(buff) < 24) { //empty line
            p = fgets (buff, LINESZ, fp);//get next line
            continue;
        }
        if (parVers < 20) {
#ifdef myUseCOut
            std::cout<<"Error: PAR files should have 'CLINICAL TRYOUT' line with a version from 2.0-4.2: " << parname<<std::endl;
#else
            printf("Error: PAR files should have 'CLINICAL TRYOUT' line with a version from 2.0-4.2: %s\n", parname);
#endif
            return d;
        }
        for (int i = 0; i < nCols; i++)
            cols[i] = strtof(p, &p); // p+1 skip comma, read a float
        if ((cols[kIndex]) != slice) isIndexSequential = false; //slices 0,1,2.. should have indices 0,1,2,3...
        slice ++;
        if (slice == 1) {
            //for (int i = 0; i < nCols; i++)
            //    cols1[i] = cols[i]; //store first slice to see if dimensions or intensity scale varies between slices
            d.xyzDim[1] = cols[kXdim];
            d.xyzDim[2] = cols[kYdim];
            d.xyzMM[1] = cols[kXmm];
            d.xyzMM[2] = cols[kYmm];
            d.xyzMM[3] = cols[kThickmm] + cols[kGapmm];
            d.patientPosition[1] = cols[kPositionRL];
            d.patientPosition[2] = cols[kPositionAP];
            d.patientPosition[3] = cols[kPositionFH];
            d.angulation[1] = cols[kAngulationRLs];
            d.angulation[2] = cols[kAngulationAPs];
            d.angulation[3] = cols[kAngulationFHs];
            d.sliceOrient = cols[kSliceOrients];
            d.TE = cols[kTEcho];
            d.bitsAllocated = cols[kBitsPerVoxel];
            d.bitsStored = cols[kBitsPerVoxel];
            d.intenIntercept = cols[kRI];
            d.intenScale = cols[kRS];
        } else {
            if ((d.xyzDim[1] != cols[kXdim]) || (d.xyzDim[2] != cols[kYdim]) || (d.bitsAllocated != cols[kBitsPerVoxel]) ) {
#ifdef myUseCOut
                std::cout<<"Error: slice dimensions or bit depth varies "<< parname <<std::endl;
#else
                printf("Error: slice dimensions or bit depth varies %s\n", parname);
#endif
                return d;
            }
            if ((d.patientPositionSequentialRepeats == 0) && ((!isSameFloat(d.patientPosition[1],cols[kPositionRL])) ||
                                                              (!isSameFloat(d.patientPosition[2],cols[kPositionAP])) ||
                                                              (!isSameFloat(d.patientPosition[3],cols[kPositionFH])) ) )//this is the first slice with different position
                d.patientPositionSequentialRepeats = slice-1;
            
            if ((d.intenScale != cols[kRS]) || (d.intenIntercept != cols[kRI]))
                isIntenScaleVaries = true;
        }
        if (cols[kImageType] == 0) d.isHasMagnitude = true;
        if (cols[kImageType] != 0) d.isHasPhase = true;
        if (cols[kGradientNumber] > 0) {
            int dir = cols[kGradientNumber];
            if ((cols[kbval] > 0.0) && (cols[kv1] == 0.0) && (cols[kv1] == 0.0) && (cols[kv1] == 0.0) ) {
                if (d.CSA.dtiV[dir-1][0] >= 0) dir = dir + 1; //Philips often stores an ADC map along with B0 and weighted images, unfortunately they give it the same kGradientNumber as the B0! (seen in PAR V4.2)
                //the logic here is that IF the gradient was previously used we increment the gradient number. This should provide compatibility when Philips fixes this bug
                //it seems like the ADC is always saved as the final volume, so this solution SHOULD be foolproof.
                ADCwarning = true;
            }
            if (dir > d.CSA.numDti) d.CSA.numDti = dir;

            if (cols[dir] <= kMaxDTIv) {
                d.CSA.dtiV[dir-1][0] = cols[kbval];
                d.CSA.dtiV[dir-1][1] = cols[kv1];
                d.CSA.dtiV[dir-1][2] = cols[kv2];
                d.CSA.dtiV[dir-1][3] = cols[kv3];

            } //save DTI direction
        } //if DTI directions
        //printf("%f %f %lu\n",cols[9],cols[kGradientNumber], strlen(buff))
        p = fgets (buff, LINESZ, fp);//get next line
    }
    free (cols);
    fclose (fp);
    d.isValid = true;
    d.isSigned = true;
    d.xyzDim[4] = slice/d.xyzDim[3];
    d.locationsInAcquisition = d.xyzDim[3];
#ifdef myUseCOut
    if (isIntenScaleVaries)
        std::cout<<"Warning: intensity slope/intercept varies between slices! [solution: user dcm2nii instead]" <<std::endl;
    
    if (!isIndexSequential)
        std::cout<<"Warning: slice order not saved to disk sequentially! [solution: user dcm2nii instead]" <<std::endl;
    printf("Warning: slice order not saved to disk sequentially! [solution: user dcm2nii instead]\n");
    std::cout<<"Done reading PAR header version "<< (float)parVers/10<<" with "<< d.CSA.numDti << "DTI directions"<<std::endl;
#else
    if (ADCwarning)
        printf("Warning: PAR/REC dataset includes an ADC map that could disrupt analysis. Please remove volume and ensure vectors are reported correctly\n");
    
    if (isIntenScaleVaries)
       printf("Warning: intensity slope/intercept varies between slices! [solution: user dcm2nii instead]\n");
        if (!isIndexSequential)
            printf("Warning: slice order not saved to disk sequentially! [solution: user dcm2nii instead]\n");
            printf("Done reading PAR header version %.1f, with %d DTI directions\n", (float)parVers/10, d.CSA.numDti);
#endif
            
            //see Xiangrui Li 's dicm2nii (also BSD license)
            // http://www.mathworks.com/matlabcentral/fileexchange/42997-dicom-to-nifti-converter
            // Rotation order and signs are figured out by try and err, not 100% sure
            float d2r = M_PI/180;
            vec3 ca = setVec3(cos(d.angulation[1]*d2r),cos(d.angulation[2]*d2r),cos(d.angulation[3]*d2r));
            vec3 sa = setVec3(sin(d.angulation[1]*d2r),sin(d.angulation[2]*d2r),sin(d.angulation[3]*d2r));
            mat33 rx,ry,rz;
    LOAD_MAT33(rx,1.0f, 0.0f, 0.0f, 0.0f, ca.v[0], -sa.v[0], 0.0f, sa.v[0], ca.v[0]);
    LOAD_MAT33(ry, ca.v[1], 0.0f, sa.v[1], 0.0f, 1.0f, 0.0f, -sa.v[1], 0.0f, ca.v[1]);
    LOAD_MAT33(rz, ca.v[2], -sa.v[2], 0.0f, sa.v[2], ca.v[2], 0.0f, 0.0f, 0.0f, 1.0f);
    mat33 R = nifti_mat33_mul( rx,ry );
    R = nifti_mat33_mul( R,rz);
    ivec3 ixyz = setVec3i(1,2,3);
    if (d.sliceOrient == kSliceOrientSag) {
        ixyz = setVec3i(2,3,1);
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                if (c != 1) R.m[r][c] = -R.m[r][c]; //invert first and final columns
    }else if (d.sliceOrient == kSliceOrientCor) {
        ixyz = setVec3i(1,3,2);
        for (int r = 0; r < 3; r++)
            R.m[r][2] = -R.m[r][2]; //invert rows of final column
    }
    R = nifti_mat33_reorder_cols(R,ixyz); //dicom rotation matrix
    d.orient[1] = R.m[0][0]; d.orient[2] = R.m[1][0]; d.orient[3] = R.m[2][0];
    d.orient[4] = R.m[0][1]; d.orient[5] = R.m[1][1]; d.orient[6] = R.m[2][1];
    mat33 diag;
    LOAD_MAT33(diag, d.xyzMM[1],0.0f,0.0f,  0.0f,d.xyzMM[2],0.0f,  0.0f,0.0f, d.xyzMM[3]);
    R= nifti_mat33_mul( R, diag );
    mat44 R44;
    LOAD_MAT44(R44, R.m[0][0],R.m[0][1],R.m[0][2],d.stackOffcentre[1],
               R.m[1][0],R.m[1][1],R.m[1][2],d.stackOffcentre[2],
               R.m[2][0],R.m[2][1],R.m[2][2],d.stackOffcentre[3]);
    vec3 x;
    if (parVers > 40) //guess
        x = setVec3(((float)d.xyzDim[1]-1)/2,((float)d.xyzDim[2]-1)/2,((float)d.xyzDim[3]-1)/2);
        else
            x = setVec3((float)d.xyzDim[1]/2,(float)d.xyzDim[2]/2,((float)d.xyzDim[3]-1)/2);
            mat44 eye;
    LOAD_MAT44(eye, 1.0f,0.0f,0.0f,x.v[0],
               0.0f,1.0f,0.0f,x.v[1],
               0.0f,0.0f,1.0f,x.v[2]);
    eye= nifti_mat44_inverse( eye ); //we wish to compute R/eye, so compute invEye and calculate R*invEye
    R44= nifti_mat44_mul( R44 , eye );
    vec4 y;
    y.v[0]=0.0f; y.v[1]=0.0f; y.v[2]=d.xyzDim[3]-1; y.v[3]=1.0f;
    y= nifti_vect44mat44_mul(y, R44 );
    int iOri = 2; //for axial, slices are 3rd dimenson (indexed from 0) (k)
    if (d.sliceOrient == kSliceOrientSag) iOri = 0; //for sagittal, slices are 1st dimension (i)
        if (d.sliceOrient == kSliceOrientCor) iOri = 1; //for coronal, slices are 2nd dimension (j)
            if  (( (y.v[iOri]-R44.m[iOri][3])>0 ) == ( (y.v[iOri]-d.stackOffcentre[iOri+1])>0 ) ) {
                d.patientPosition[1] = R44.m[0][3];
                d.patientPosition[2] = R44.m[1][3];
                d.patientPosition[3] = R44.m[2][3];
                d.patientPositionLast[1] = y.v[0];
                d.patientPositionLast[2] = y.v[1];
                d.patientPositionLast[3] = y.v[2];
            }else {
                //d.patientPosition
                d.patientPosition[1] = y.v[0];
                d.patientPosition[2] = y.v[1];
                d.patientPosition[3] = y.v[2];
                d.patientPositionLast[1] = R44.m[0][3];
                d.patientPositionLast[2] = R44.m[1][3];
                d.patientPositionLast[3] = R44.m[2][3];
            }
    //finish up
    changeExt (parname, "REC");
    d.locationsInAcquisition = d.xyzDim[3];
    d.manufacturer = kMANUFACTURER_PHILIPS;
    d.imageStart = 0;
    return d;
} //nii_readParRec()

size_t nii_ImgBytes(struct nifti_1_header hdr) {
    //unsigned long imgsz = nii_ImgBytes(hdr);
    size_t imgsz = hdr.bitpix/8;
    for (int i = 1; i < 8; i++)
        if (hdr.dim[i]  > 1)
            imgsz = imgsz * hdr.dim[i];
    return imgsz;
} //nii_ImgBytes()

unsigned char * nii_demosaic(unsigned char* inImg, struct nifti_1_header *hdr, int nMosaicSlices, int ProtocolSliceNumber1) {
    //demosaic http://nipy.org/nibabel/dicom/dicom_mosaic.html
    if (nMosaicSlices < 2) return inImg;
    //Byte inImg[ [img length] ];
    //[img getBytes:&inImg length:[img length]];
    int nRowCol = ceil(sqrt(nMosaicSlices));
    int colBytes = hdr->dim[1]/nRowCol * hdr->bitpix/8;
    int lineBytes = hdr->dim[1] * hdr->bitpix/8;
    int rowBytes = hdr->dim[1] * hdr->dim[2]/nRowCol * hdr->bitpix/8;
    int col = 0;
    int row = 0;
    int lOutPos = 0;
    hdr->dim[1] = hdr->dim[1]/nRowCol;
    hdr->dim[2] = hdr->dim[2]/nRowCol;
    hdr->dim[3] = nMosaicSlices;
    size_t imgsz = nii_ImgBytes(*hdr);
    unsigned char *outImg = (unsigned char *)malloc(imgsz);
    for (int m=1; m <= nMosaicSlices; m++) {
        int lPos = (row * rowBytes) + (col * colBytes);
        for (int y = 0; y < hdr->dim[2]; y++) {
            memcpy(&outImg[lOutPos], &inImg[lPos], colBytes); // dest, src, bytes
            lPos += lineBytes;
            lOutPos +=colBytes;
        }
        col ++;
        if (col >= nRowCol) {
            row ++;
            col = 0;
        } //start new column
    } //for m = each mosaic slice
    if (ProtocolSliceNumber1 > 1) {
#ifdef myUseCOut
     	std::cout<<"WARNING: CSA 'ProtocolSliceNumber' SUGGESTS REVERSED SLICE ORDER: SPATIAL AND DTI COORDINATES UNTESTED" <<std::endl;
#else
        printf("WARNING: CSA 'ProtocolSliceNumber' SUGGESTS REVERSED SLICE ORDER: SPATIAL AND DTI COORDINATES UNTESTED\n");
#endif
    }
    /*if ((ProtocolSliceNumber1 > 1) && (hdr->dim[3] > 1)) { //exceptionally rare: reverse order of slices - now handled in matrix...
     int sliceBytes = hdr->dim[1] * hdr->dim[2] * hdr->bitpix/8;
     memcpy(&inImg[0], &outImg[0],sliceBytes*hdr->dim[3]); //copy data with reversed order dest, src, bytes
     int lOutPos = sliceBytes * (hdr->dim[3]-1);
     int lPos = 0;
     for (int m=0; m < nMosaicSlices; m++) {
     memcpy( &outImg[lOutPos], &inImg[lPos], sliceBytes);
     lPos += sliceBytes;
     lOutPos -= sliceBytes;
     }
     }*/
    free(inImg);
    return outImg;
    //return [NSData dataWithBytes:&outImg length:outlen];
} //nii_demosaic()

unsigned char * nii_flipImgY(unsigned char* bImg, struct nifti_1_header *hdr){
    //DICOM row order opposite from NIfTI
    int dim3to7 = 1;
    for (int i = 3; i < 8; i++)
        if (hdr->dim[i] > 1) dim3to7 = dim3to7 * hdr->dim[i];
    size_t lineBytes = hdr->dim[1] * hdr->bitpix/8;
    if ((hdr->datatype == DT_RGB24) && (hdr->bitpix == 24)) {
        lineBytes = hdr->dim[1];
        dim3to7 = dim3to7 * 3;
    } //rgb data saved planar (RRR..RGGGG..GBBB..B
    
    //printf("smoking gun %d %d\n",hdr->dim[3],hdr->dim[4]); return bImg;
    
    unsigned char line[lineBytes];
    size_t sliceBytes = hdr->dim[2] * lineBytes;
    int halfY = hdr->dim[2] / 2; //note truncated toward zero, so halfY=2 regardless of 4 or 5 columns
    for (int sl = 0; sl < dim3to7; sl++) { //for each 2D slice
        size_t slBottom = (size_t)sl*sliceBytes;
        size_t slTop = (((size_t)sl+1)*sliceBytes)-lineBytes;
        for (int y = 0; y < halfY; y++) {
            //swap order of lines
            memcpy(&line, &bImg[slBottom], lineBytes);
            memcpy(&bImg[slBottom], &bImg[slTop], lineBytes);
            memcpy(&bImg[slTop], &line, lineBytes);
            slTop -= lineBytes;
            slBottom += lineBytes;
        } //for y
    } //for each slice
    return bImg;
} //nii_flipImgY()

unsigned char * nii_flipImgZ(unsigned char* bImg, struct nifti_1_header *hdr){
    //DICOM row order opposite from NIfTI
    int halfZ = hdr->dim[3] / 2; //note truncated toward zero, so halfY=2 regardless of 4 or 5 columns
    if (halfZ < 1) return bImg;
    int dim4to7 = 1;
    for (int i = 4; i < 8; i++)
        if (hdr->dim[i] > 1) dim4to7 = dim4to7 * hdr->dim[i];
    int sliceBytes = hdr->dim[1] * hdr->dim[2] * hdr->bitpix/8;
    long long volBytes = sliceBytes * hdr->dim[3];
    unsigned char slice[sliceBytes];
    
    for (int vol = 0; vol < dim4to7; vol++) { //for each 2D slice
        long long slBottom = vol*volBytes;
        long long slTop = ((vol+1)*volBytes)-sliceBytes;
        for (int z = 0; z < halfZ; z++) {
            //swap order of lines
            memcpy(&slice, &bImg[slBottom], sliceBytes);
            memcpy(&bImg[slBottom], &bImg[slTop], sliceBytes);
            memcpy(&bImg[slTop], &slice, sliceBytes);
            slTop -= sliceBytes;
            slBottom += sliceBytes;
        } //for Z
    } //for each volume
    return bImg;
} //nii_flipImgZ()

unsigned char * nii_flipZ(unsigned char* bImg, struct nifti_1_header *h){
    //flip slice order
    if (h->dim[3] < 2) return bImg;
    mat33 s;
    mat44 Q44;
    LOAD_MAT33(s,h->srow_x[0],h->srow_x[1],h->srow_x[2], h->srow_y[0],h->srow_y[1],h->srow_y[2],
               h->srow_z[0],h->srow_z[1],h->srow_z[2]);
    LOAD_MAT44(Q44,h->srow_x[0],h->srow_x[1],h->srow_x[2],h->srow_x[3],
               h->srow_y[0],h->srow_y[1],h->srow_y[2],h->srow_y[3],
               h->srow_z[0],h->srow_z[1],h->srow_z[2],h->srow_z[3]);
    vec4 v= setVec4(0,0,h->dim[3]-1);
    v = nifti_vect44mat44_mul(v, Q44); //after flip this voxel will be the origin
    mat33 mFlipZ;
    LOAD_MAT33(mFlipZ,1.0f, 0.0f, 0.0f, 0.0f,1.0f,0.0f, 0.0f,0.0f,-1.0f);
    s= nifti_mat33_mul( s , mFlipZ );
    LOAD_MAT44(Q44, s.m[0][0],s.m[0][1],s.m[0][2],v.v[0],
               s.m[1][0],s.m[1][1],s.m[1][2],v.v[1],
               s.m[2][0],s.m[2][1],s.m[2][2],v.v[2]);
    //printf(" ----------> %f %f %f\n",v.v[0],v.v[1],v.v[2]);
    setQSForm(h,Q44);
    //printf("nii_flipImgY dims %dx%dx%d %d \n",h->dim[1],h->dim[2], dim3to7,h->bitpix/8);
    return nii_flipImgZ(bImg,h);
}//nii_flipZ()


unsigned char * nii_flipY(unsigned char* bImg, struct nifti_1_header *h){
    mat33 s;
    mat44 Q44;
    LOAD_MAT33(s,h->srow_x[0],h->srow_x[1],h->srow_x[2], h->srow_y[0],h->srow_y[1],h->srow_y[2],
               h->srow_z[0],h->srow_z[1],h->srow_z[2]);
    LOAD_MAT44(Q44,h->srow_x[0],h->srow_x[1],h->srow_x[2],h->srow_x[3],
               h->srow_y[0],h->srow_y[1],h->srow_y[2],h->srow_y[3],
               h->srow_z[0],h->srow_z[1],h->srow_z[2],h->srow_z[3]);
    vec4 v= setVec4(0,h->dim[2]-1,0);
    v = nifti_vect44mat44_mul(v, Q44); //after flip this voxel will be the origin
    mat33 mFlipY;
    LOAD_MAT33(mFlipY,1.0f, 0.0f, 0.0f, 0.0f,-1.0f,0.0f, 0.0f,0.0f,1.0f);
    
    s= nifti_mat33_mul( s , mFlipY );
    LOAD_MAT44(Q44, s.m[0][0],s.m[0][1],s.m[0][2],v.v[0],
               s.m[1][0],s.m[1][1],s.m[1][2],v.v[1],
               s.m[2][0],s.m[2][1],s.m[2][2],v.v[2]);
    setQSForm(h,Q44);
    //printf("nii_flipImgY dims %dx%dx%d %d \n",h->dim[1],h->dim[2], dim3to7,h->bitpix/8);
    return nii_flipImgY(bImg,h);
}//nii_flipY()

unsigned char * nii_loadImgCore(char* imgname, struct nifti_1_header hdr) {
    size_t imgsz = nii_ImgBytes(hdr);
    FILE *file = fopen(imgname , "rb");
	if (!file) {
         printf("Error: unable to open %s\n", imgname);
         return NULL;
    }
	fseek(file, 0, SEEK_END);
	long long fileLen=ftell(file);
    if (fileLen < (imgsz+hdr.vox_offset)) return NULL;
	fseek(file, hdr.vox_offset, SEEK_SET);
    unsigned char *bImg = (unsigned char *)malloc(imgsz);
	fread(bImg, imgsz, 1, file);
	fclose(file);
    return bImg;
} //nii_loadImg()

unsigned char * nii_rgb2Planar(unsigned char* bImg, struct nifti_1_header *hdr, int isPlanar) {
    //DICOM data saved in triples RGBRGBRGB, NIfTI RGB saved in planes RRR..RGGG..GBBBB..B
    if (hdr->datatype != DT_RGB24) return bImg;
    if (isPlanar == 1) return bImg;//return nii_bgr2rgb(bImg,hdr);
    int dim3to7 = 1;
    for (int i = 3; i < 8; i++)
        if (hdr->dim[i] > 1) dim3to7 = dim3to7 * hdr->dim[i];
    int sliceBytes24 = hdr->dim[1]*hdr->dim[2] * hdr->bitpix/8;
    int sliceBytes8 = hdr->dim[1]*hdr->dim[2];
    unsigned char  slice24[ sliceBytes24 ];
    int sliceOffsetR = 0;
    for (int sl = 0; sl < dim3to7; sl++) { //for each 2D slice
        memcpy(&slice24, &bImg[sliceOffsetR], sliceBytes24);
        int sliceOffsetG = sliceOffsetR + sliceBytes8;
        int sliceOffsetB = sliceOffsetR + 2*sliceBytes8;
        int i = 0;
        int j = 0;
        for (int rgb = 0; rgb < sliceBytes8; rgb++) {
            bImg[sliceOffsetR+j] =slice24[i];
            i++;
            bImg[sliceOffsetG+j] =slice24[i];
            i++;
            bImg[sliceOffsetB+j] =slice24[i];
            i++;
            j++;
        }
        sliceOffsetR += sliceBytes24;
    } //for each slice
    return bImg;
} //nii_rgb2Planar()

unsigned char * nii_iVaries(unsigned char *img, struct nifti_1_header *hdr){
    //each DICOM image can have its own intesity scaling, whereas NIfTI requires the same scaling for all images in a file
    //WARNING: do this BEFORE nii_check16bitUnsigned!!!!
    //if (hdr->datatype != DT_INT16) return img;
    int dim3to7 = 1;
    for (int i = 3; i < 8; i++)
        if (hdr->dim[i] > 1) dim3to7 = dim3to7 * hdr->dim[i];
    int nVox = hdr->dim[1]*hdr->dim[2]* dim3to7;
    if (nVox < 1) return img;
    float * img32=(float*)malloc(nVox*sizeof(float));
    if (hdr->datatype == DT_UINT8) {
        uint8_t * img8i = (uint8_t*) img;
        for (int i=0; i < nVox; i++)
            img32[i] = img8i[i];
    } else if (hdr->datatype == DT_UINT16) {
        uint16_t * img16ui = (uint16_t*) img;
        for (int i=0; i < nVox; i++)
            img32[i] = img16ui[i];
    } else if (hdr->datatype == DT_INT16) {
        int16_t * img16i = (int16_t*) img;
        for (int i=0; i < nVox; i++)
            img32[i] = img16i[i];
    } else if (hdr->datatype == DT_INT32) {
        int32_t * img32i = (int32_t*) img;
        for (int i=0; i < nVox; i++)
            img32[i] = img32i[i];
    }
    free (img); //release previous image
    for (int i=0; i < nVox; i++)
        img32[i] = (img32[i]* hdr->scl_slope)+hdr->scl_inter;
    hdr->scl_slope = 1;
    hdr->scl_inter = 0;
    hdr->datatype = DT_FLOAT;
    hdr->bitpix = 32;
    return (unsigned char*) img32;
} //nii_iVaries()

unsigned char * nii_XYTZ_XYZT(unsigned char* bImg, struct nifti_1_header *hdr, int seqRepeats) {
    //Philips can save time as 3rd dimensions, NIFTI requires time is 4th dimension
    int dim4to7 = 1;
    for (int i = 4; i < 8; i++)
        if (hdr->dim[i] > 1) dim4to7 = dim4to7 * hdr->dim[i];
    if ((hdr->dim[3] < 2) || (dim4to7 < 2)) return bImg;
#ifdef myUseCOut
    std::cout<<"Converting XYTZ to XYZT with "<<hdr->dim[3]<<" slices (Z) and "<<  dim4to7<< "volumes" <<std::endl;
    if ((dim4to7 % seqRepeats) != 0) {
        std::cout<<"Error: patient position repeats "<<seqRepeats<<"  times, but this does not evenly divide number of volumes: "<<  dim4to7 <<std::endl;
        seqRepeats = 1;
    }
    
#else
    printf("Converting XYTZ to XYZT with %d slices (Z) and %d volumes (T).\n",hdr->dim[3], dim4to7);
    if ((dim4to7 % seqRepeats) != 0) {
        printf("Error: patient position repeats %d times, but this does not evenly divide number of volumes (%d)", seqRepeats,dim4to7);
        seqRepeats = 1;
    }
    
#endif
    uint64_t typeRepeats = dim4to7 / seqRepeats;
    uint64_t sliceBytes = hdr->dim[1]*hdr->dim[2]*hdr->bitpix/8;
    uint64_t seqBytes = sliceBytes * seqRepeats;
    uint64_t typeBytes = seqBytes * hdr->dim[3];
    
    uint64_t imgSz = nii_ImgBytes(*hdr);
    //this uses a lot of RAM, someday this could be done in place...
    unsigned char *outImg = (unsigned char *)malloc(imgSz);
    //memcpy(&tempImg[0], &bImg[0], imgSz);
    uint64_t origPos = 0;
    uint64_t Pos = 0; //
    
    for (int t = 0; t < typeRepeats; t++) { //for each volume
        for (int s = 0; s < seqRepeats; s++) {
            origPos = (t*typeBytes) +s*sliceBytes;
            for (int z = 0; z < hdr->dim[3]; z++) { //for each slice
                memcpy( &outImg[Pos],&bImg[origPos], sliceBytes);
                Pos += sliceBytes;
                origPos += seqBytes;
            }
        }//for s
    }
    free(bImg);
    return outImg;
} //nii_ImgBytes()

unsigned char * nii_byteswap(unsigned char *img, struct nifti_1_header *hdr){
    if (hdr->bitpix < 9) return img;
    uint64_t nvox = nii_ImgBytes(*hdr) / (hdr->bitpix/8);
    void *ar = (void*) img;
    if (hdr->bitpix == 16) nifti_swap_2bytes( nvox , ar );
    if (hdr->bitpix == 32) nifti_swap_4bytes( nvox , ar );
    if (hdr->bitpix == 64) nifti_swap_8bytes( nvox , ar );
    return img;
} //nii_byteswap()

unsigned char * nii_loadImgX(char* imgname, struct nifti_1_header *hdr, struct TDICOMdata dcm, bool iVaries) {
    //provided with a filename (imgname) and DICOM header (dcm), creates NIfTI header (hdr) and img

    if (headerDcm2Nii(dcm, hdr) == EXIT_FAILURE) return NULL;
    unsigned char * img = nii_loadImgCore(imgname, *hdr);
    if (img == NULL) return img;
#ifdef __BIG_ENDIAN__
    if ((dcm.isLittleEndian) && (hdr->bitpix > 8))
        img = nii_byteswap(img, hdr);
#else
    if ((!dcm.isLittleEndian) && (hdr->bitpix > 8))
        img = nii_byteswap(img, hdr);
#endif
    if (hdr->datatype ==DT_RGB24) img = nii_rgb2Planar(img, hdr, dcm.isPlanarRGB);//do this BEFORE Y-Flip, or RGB order can be flipped
    if (dcm.CSA.mosaicSlices > 1) {
        img = nii_demosaic(img, hdr, dcm.CSA.mosaicSlices, dcm.CSA.protocolSliceNumber1);
        /* we will do this in nii_dicom_batch #ifdef obsolete_mosaic_flip
         img = nii_flipImgY(img, hdr);
         #endif*/
    }
    if (iVaries) img = nii_iVaries(img, hdr);
    int nAcq = dcm.locationsInAcquisition;
    if ((nAcq > 1) && (hdr->dim[0] < 4) && ((hdr->dim[3]%nAcq)==0) && (hdr->dim[3]>nAcq) ) {
        hdr->dim[4] = hdr->dim[3]/nAcq;
        hdr->dim[3] = nAcq;
        hdr->dim[0] = 4;
    }
    if ((hdr->dim[0] > 3) && (dcm.patientPositionSequentialRepeats > 1)) //swizzle 3rd and 4th dimension (Philips stores time as 3rd dimension)
        img = nii_XYTZ_XYZT(img, hdr,dcm.patientPositionSequentialRepeats );
    headerDcm2NiiSForm(dcm,dcm, hdr);
    return img;
} //nii_loadImgX()

struct TDICOMdata readDICOMv(char * fname, bool isVerbose) {
    struct TDICOMdata d = clear_dicom_data();
    strcpy(d.protocolName, ""); //fill dummy with empty space so we can detect kProtocolNameGE
    FILE *file = fopen(fname, "rb");
	if (!file) {
#ifdef myUseCOut
     	std::cout<<"Unable to open file "<< fname <<std::endl;
#else
        printf( "Unable to open file %s\n", fname);
#endif
		
		return d;
	}
	fseek(file, 0, SEEK_END);
	long long fileLen=ftell(file); //Get file length
    if (fileLen < 256) {
#ifdef myUseCOut
     	std::cout<<"File too small to be a DICOM image "<< fname <<std::endl;
#else
        printf( "File too small to be a DICOM image %s\n", fname);
#endif
		return d;
	}
	fseek(file, 0, SEEK_SET);
	//Allocate memory
	unsigned char *buffer=(unsigned char *)malloc(fileLen+1);
	if (!buffer) {
		printf( "Memory error!");
        fclose(file);
		return d;
	}
	//Read file contents into buffer
	fread(buffer, fileLen, 1, file);
	fclose(file);
    if ((buffer[128] != 'D') || (buffer[129] != 'I')  || (buffer[130] != 'C') || (buffer[131] != 'M')) {
        free (buffer);
        return d;
    }
    //DEFINE DICOM TAGS
#define  kUnused 0x0001+(0x0001 << 16 )
#define  kStart 0x0002+(0x0000 << 16 )
#define  kTransferSyntax 0x0002+(0x0010 << 16)
#define  kStudyDate 0x0008+(0x0020 << 16 )
#define  kStudyTime 0x0008+(0x0030 << 16 )
#define  kAcquisitionTime 0x0008+(0x0032 << 16 )
#define  kManufacturer 0x0008+(0x0070 << 16 )
#define  kProtocolNameGE 0x0008+(0x103E << 16 )
#define  kComplexImageComponent 0x0008+(0x9208 << 16 )//'0008' '9208' 'CS' 'ComplexImageComponent'
#define  kPatientName 0x0010+(0x0010 << 16 )
#define  kPatientID 0x0010+(0x0020 << 16 )
#define  kMRAcquisitionType 0x0018+(0x0023 << 16)
#define  kZThick  0x0018+(0x0050 << 16 )
#define  kTR  0x0018+(0x0080 << 16 )
#define  kTE  0x0018+(0x0081 << 16 )
#define  kEchoNum  0x0018+(0x0086 << 16 ) //IS
#define  kZSpacing  0x0018+(0x0088 << 16 ) //'DS' 'SpacingBetweenSlices'
#define  kProtocolName  0x0018+(0x1030<< 16 )
#define  kInPlanePhaseEncodingDirection  0x0018+(0x1312<< 16 ) //CS
#define  kPatientOrient  0x0018+(0x5100<< 16 )    //0018,5100. patient orientation - 'HFS'
    //#define  kDiffusionBFactorSiemens  0x0019+(0x100C<< 16 ) //   0019;000C;SIEMENS MR HEADER  ;B_value                         ;1;IS;1
#define  kLastScanLoc  0x0019+(0x101B<< 16 )
#define  kDiffusionDirectionGEX  0x0019+(0x10BB<< 16 ) //DS
#define  kDiffusionDirectionGEY  0x0019+(0x10BC<< 16 ) //DS
#define  kDiffusionDirectionGEZ  0x0019+(0x10BD<< 16 ) //DS
#define  kPatientPosition 0x0020+(0x0032 << 16 )
#define  kSeriesNum 0x0020+(0x0011 << 16 )
#define  kAcquNum 0x0020+(0x0012 << 16 )
#define  kImageNum 0x0020+(0x0013 << 16 )
#define  kOrientation 0x0020+(0x0037 << 16 )
#define  kImageComments 0x0020+(0x4000<< 16 )// '0020' '4000' 'LT' 'ImageComments'
#define  kLocationsInAcquisitionGE 0x0021+(0x104F<< 16 )// 'SS' 'LocationsInAcquisitionGE'
#define  kSamplesPerPixel 0x0028+(0x0002 << 16 )
#define  kPlanarRGB 0x0028+(0x0006 << 16 )
#define  kDim3 0x0028+(0x0008 << 16 ) //number of frames - for Philips this is Dim3*Dim4
#define  kDim2 0x0028+(0x0010 << 16 )
#define  kDim1 0x0028+(0x0011 << 16 )
#define  kXYSpacing  0x0028+(0x0030 << 16 ) //'0028' '0030' 'DS' 'PixelSpacing'
#define  kBitsAllocated 0x0028+(0x0100 << 16 )
#define  kBitsStored 0x0028+(0x0101 << 16 )//'0028' '0101' 'US' 'BitsStored'
#define  kIsSigned 0x0028+(0x0103 << 16 )
#define  kIntercept 0x0028+(0x1052 << 16 )
#define  kSlope 0x0028+(0x1053 << 16 )
#define  kGeiisFlag 0x0029+(0x0010 << 16 ) //warn user if dreaded GEIIS was used to process image
#define  kCSAImageHeaderInfo  0x0029+(0x1010 << 16 )
    //#define  kObjectGraphics  0x0029+(0x1210 << 16 )    //0029,1210 syngoPlatformOOGInfo Object Oriented Graphics
#define  kDiffusionBFactorGE  0x0043+(0x1039 << 16 ) //IS dicm2nii's SlopInt_6_9
#define  kCoilSiemens  0x0051+(0x100F << 16 )
#define  kLocationsInAcquisition  0x0054+(0x0081 << 16 )
#define  kIconImageSequence 0x0088+(0x0200 << 16 )
#define  kDiffusionBFactor  0x2001+(0x1003 << 16 )// FL
#define  kSliceOrient  0x2001+(0x100B << 16 )//2001,100B Philips slice orientation (TRANSVERSAL, AXIAL, SAGITTAL)
#define  kLocationsInAcquisitionPhilips  0x2001+(0x1018 << 16 )
#define  kNumberOfDynamicScans  0x2001+(0x1081 << 16 )//'2001' '1081' 'IS' 'NumberOfDynamicScans'
#define  kMRAcquisitionTypePhilips 0x2005+(0x106F << 16)
#define  kAngulationAP 0x2005+(0x1071 << 16)//'2005' '1071' 'FL' 'MRStackAngulationAP'
#define  kAngulationFH 0x2005+(0x1072 << 16)//'2005' '1072' 'FL' 'MRStackAngulationFH'
#define  kAngulationRL 0x2005+(0x1073 << 16)//'2005' '1073' 'FL' 'MRStackAngulationRL'
#define  kMRStackOffcentreAP 0x2005+(0x1078 << 16)
#define  kMRStackOffcentreFH 0x2005+(0x1079 << 16)
#define  kMRStackOffcentreRL 0x2005+(0x107A << 16)
#define  kDiffusionDirectionRL 0x2005+(0x10B0 << 16)
#define  kDiffusionDirectionAP 0x2005+(0x10B1 << 16)
#define  kDiffusionDirectionFH 0x2005+(0x10B2 << 16)
#define  k2005140F 0x2005+(0x140F << 16)
#define  kWaveformSq 0x5400+(0x0100 << 16)
#define  kImageStart 0x7FE0+(0x0010 << 16 )
#define  kNest 0xFFFE +(0xE000 << 16 ) //Item follows SQ
#define  kUnnest 0xFFFE +(0xE00D << 16 ) //ItemDelimitationItem [length defined] http://www.dabsoft.ch/dicom/5/7.5/
#define  kUnnest2 0xFFFE +(0xE0DD << 16 )//SequenceDelimitationItem [length undefined]
    double zSpacing = -1.0l; //includes slice thickness plus gap
    int locationsInAcquisitionGE = 0; int locationsInAcquisitionPhilips = 0;
    long lPos = 128+4; //4-byte signature starts at 128
    uint32_t lLength;
    uint32_t groupElement = buffer[lPos] | (buffer[lPos+1] << 8) | (buffer[lPos+2] << 16) | (buffer[lPos+3] << 24);
    if (groupElement != kStart)
#ifdef myUseCOut
     	std::cout<<"DICOM appears corrupt: first group:element should be 0x0002:0x0000" <<std::endl;
#else
    printf("DICOM appears corrupt: first group:element should be 0x0002:0x0000\n");
#endif
    char vr[2];
    bool isIconImageSequence = false;
    bool isSwitchToImplicitVR = false;
    bool isSwitchToBigEndian = false;
    //bool geiisBug = false; //for buggy GEIIS http://forum.dcmtk.org/viewtopic.php?p=7162&sid=3b516cc751aae51fbb5e73184abe37c2
    bool is2005140FSQ = false; //for buggy Philips
    bool is2005140FSQwarned = false; //for buggy Philips
    bool isAtFirstPatientPosition = false; //for 3d and 4d files: flag is true for slices at same position as first slice
    int patientPositionCount = 0;
    long coilNum = 0; //Siemens can save one image per coil (H12,H13,etc) or one combined image for array (HEA;HEP)
    long echoNum = 0;
    while ((d.imageStart == 0) && ((lPos+8) <  fileLen)) {
        if (d.isLittleEndian)
            groupElement = buffer[lPos] | (buffer[lPos+1] << 8) | (buffer[lPos+2] << 16) | (buffer[lPos+3] << 24);
        else
            groupElement = buffer[lPos+1] | (buffer[lPos] << 8) | (buffer[lPos+3] << 16) | (buffer[lPos+2] << 24);
        if ((isSwitchToBigEndian) && ((groupElement & 0xFFFF) != 2)) {
            isSwitchToBigEndian = false;
            d.isLittleEndian = false;
            groupElement = buffer[lPos+1] | (buffer[lPos] << 8) | (buffer[lPos+3] << 16) | (buffer[lPos+2] << 24);
        }//transfer syntax requests switching endian after group 0002
        if ((isSwitchToImplicitVR) && ((groupElement & 0xFFFF) != 2)) {
            isSwitchToImplicitVR = false;
            d.isExplicitVR = false;
        } //transfer syntax requests switching VR after group 0001
        //uint32_t group = (groupElement & 0xFFFF);
        lPos += 4;
        if ((groupElement == kNest) || (groupElement == kUnnest) || (groupElement == kUnnest2)) {
            vr[0] = 'N';
            vr[1] = 'A';
            
            //if (groupElement == kUnnest) geiisBug = false; //don't exit if there is a proprietary thumbnail
            //printf("xxx");
            lLength = 4;
        } else if (d.isExplicitVR) {
            vr[0] = buffer[lPos]; vr[1] = buffer[lPos+1];
            if (buffer[lPos+1] < 'A') {//implicit vr with 32-bit length
                if (d.isLittleEndian)
                    lLength = buffer[lPos] | (buffer[lPos+1] << 8) | (buffer[lPos+2] << 16) | (buffer[lPos+3] << 24);
                else
                    lLength = buffer[lPos+3] | (buffer[lPos+2] << 8) | (buffer[lPos+1] << 16) | (buffer[lPos] << 24);
                lPos += 4;
            } else if ( ((buffer[lPos] == 'U') && (buffer[lPos+1] == 'N'))
                       || ((buffer[lPos] == 'O') && (buffer[lPos+1] == 'B'))
                       || ((buffer[lPos] == 'O') && (buffer[lPos+1] == 'W'))
                       ) { //VR= UN, OB, OW, SQ  || ((buffer[lPos] == 'S') && (buffer[lPos+1] == 'Q'))
                lPos = lPos + 4;  //skip 2 byte VR string and 2 reserved bytes = 4 bytes
                if (d.isLittleEndian)
                    lLength = buffer[lPos] | (buffer[lPos+1] << 8) | (buffer[lPos+2] << 16) | (buffer[lPos+3] << 24);
                else
                    lLength = buffer[lPos+3] | (buffer[lPos+2] << 8) | (buffer[lPos+1] << 16) | (buffer[lPos] << 24);
                lPos = lPos + 4;  //skip 4 byte length
            } else if   ((buffer[lPos] == 'S') && (buffer[lPos+1] == 'Q')) {
                lLength = 8; //Sequence Tag
                is2005140FSQ = (groupElement == k2005140F);
            } else { //explicit VR with 16-bit length
                if ((d.isLittleEndian)  )
                    lLength = buffer[lPos+2] | (buffer[lPos+3] << 8);
                else
                    lLength = buffer[lPos+3] | (buffer[lPos+2] << 8);
                lPos += 4;  //skip 2 byte VR string and 2 length bytes = 4 bytes
            }
        } else { //implicit VR
            if (d.isLittleEndian)
                lLength = buffer[lPos] | (buffer[lPos+1] << 8) | (buffer[lPos+2] << 16) | (buffer[lPos+3] << 24);
            else
                lLength = buffer[lPos+3] | (buffer[lPos+2] << 8) | (buffer[lPos+1] << 16) | (buffer[lPos] << 24);
            lPos += 4;  //we have loaded the 32-bit length
        } //if explicit else implicit VR
        if (lLength == 0xFFFFFFFF) lLength = 8; //SQ (Sequences) use 0xFFFFFFFF [4294967295] to denote unknown length
        //next: look for required tags
        if ((isIconImageSequence) && ((groupElement & 0x0028) == 0x0028 )) groupElement = kUnused; //ignore icon dimensions
        switch ( groupElement ) {
            case 	kTransferSyntax: {
                char transferSyntax[kDICOMStr];
                dcmStr (lLength, &buffer[lPos], transferSyntax);
                //printf("transfer syntax '%s'\n", transferSyntax);
                if (strcmp(transferSyntax, "1.2.840.10008.1.2.1") == 0)
                    ; //default isExplicitVR=true; //d.isLittleEndian=true
                else if (strcmp(transferSyntax, "1.2.840.10008.1.2.2") == 0)
                    isSwitchToBigEndian = true; //isExplicitVR=true;
                else if (strcmp(transferSyntax, "1.2.840.10008.1.2") == 0)
                    isSwitchToImplicitVR = true; //d.isLittleEndian=true
                else {
#ifdef myUseCOut
                    std::cout<<"Unsupported transfer syntax "<< transferSyntax<<std::endl;
#else
                    printf("Unsupported transfer syntax '%s'\n",transferSyntax);
#endif
                    
                    d.imageStart = 1;//abort as invalid (imageStart MUST be >128)
                }
                break;} //{} provide scope for variable 'transferSyntax
            case 	kStudyDate:
                dcmStr (lLength, &buffer[lPos], d.studyDate);
                break;
            case 	kManufacturer:
                d.manufacturer = dcmStrManufacturer (lLength, &buffer[lPos]);
                break;
            case 	kComplexImageComponent:
                d.isHasPhase = (buffer[lPos]=='P') && (toupper(buffer[lPos+1]) == 'H');
                d.isHasMagnitude = (buffer[lPos]=='M') && (toupper(buffer[lPos+1]) == 'A');
                break;
            case 	kAcquisitionTime : {
                char acquisitionTimeTxt[kDICOMStr];
                dcmStr (lLength, &buffer[lPos], acquisitionTimeTxt);
                d.acquisitionTime = atof(acquisitionTimeTxt);
                
                break; }
            case 	kStudyTime :
                dcmStr (lLength, &buffer[lPos], d.studyTime);
                break;
            case 	kPatientName :
                dcmStr (lLength, &buffer[lPos], d.patientName);
                break;
            case 	kPatientID :
                dcmStr (lLength, &buffer[lPos], d.patientID);
                break;
            case 	kProtocolNameGE: {
                if (d.manufacturer == kMANUFACTURER_GE)
                    dcmStr (lLength, &buffer[lPos], d.protocolName);
                break; }
                
            case 	kProtocolName : {
                if (strlen(d.protocolName) < 1) //GE uses a generic session name here: do not overwrite kProtocolNameGE
                    dcmStr (lLength, &buffer[lPos], d.protocolName);
                break; }
            case 	kPatientOrient :
                dcmStr (lLength, &buffer[lPos], d.patientOrient);
                break;
            case 	kLastScanLoc :
                d.lastScanLoc = dcmStrFloat(lLength, &buffer[lPos]);
                break;
                /*case kDiffusionBFactorSiemens :
                 if (d.manufacturer == kMANUFACTURER_SIEMENS)
                 printf(">>>>%f\n,",dcmStrFloat(lLength, &buffer[lPos]));
                 
                 break;*/
            case kDiffusionDirectionGEX :
                if (d.manufacturer == kMANUFACTURER_GE)  d.CSA.dtiV[0][1] =  dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case kDiffusionDirectionGEY :
                if (d.manufacturer == kMANUFACTURER_GE)  d.CSA.dtiV[0][2] =  dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case kDiffusionDirectionGEZ :
                if (d.manufacturer == kMANUFACTURER_GE) {
                    d.CSA.dtiV[0][3] =  dcmStrFloat(lLength, &buffer[lPos]);
                    d.CSA.numDti = 1;
                }
                break;
            case 	kPatientPosition :
                if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (is2005140FSQ)) {
#ifdef myUseCOut
                    if (!is2005140FSQwarned)
                        std::cout<<"Warning: Philips R3.2.2 can report different positions for the same slice. Attempting patch." <<std::endl;
#else
                    if (!is2005140FSQwarned)
                        printf("Warning: Philips R3.2.2 can report different positions for the same slice. Attempting patch.\n");
#endif
                    
                    is2005140FSQwarned = true;
                } else {
                    patientPositionCount++;
                    isAtFirstPatientPosition = true;
                    if (isnan(d.patientPosition[1]))
                        dcmMultiFloat(lLength, (char*)&buffer[lPos], 3, &d.patientPosition[0]); //slice position
                    else {
                        dcmMultiFloat(lLength, (char*)&buffer[lPos], 3, &d.patientPositionLast[0]); //slice direction for 4D
                        if ((isFloatDiff(d.patientPositionLast[1],d.patientPosition[1]))  ||
                            (isFloatDiff(d.patientPositionLast[2],d.patientPosition[2]))  ||
                            (isFloatDiff(d.patientPositionLast[3],d.patientPosition[3])) ) {
                            isAtFirstPatientPosition = false; //this slice is not at position of 1st slice
                            if (d.patientPositionSequentialRepeats == 0) //this is the first slice with different position
                                d.patientPositionSequentialRepeats = patientPositionCount-1;
                        } //if different position from 1st slice in file
                    } //if not first slice in file
                } //not after 2005,140F
                break;
            case 	kInPlanePhaseEncodingDirection:
                d.phaseEncodingRC = toupper(buffer[lPos]); //first character is either 'R'ow or 'C'ol
                break;
            case 	kSeriesNum:
                d.seriesNum =  dcmStrInt(lLength, &buffer[lPos]);
                break;
            case 	kAcquNum:
                d.acquNum = dcmStrInt(lLength, &buffer[lPos]);
                break;
            case 	kImageNum:
                d.imageNum = dcmStrInt(lLength, &buffer[lPos]);
                break;
            case 	kPlanarRGB:
                d.isPlanarRGB = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kDim3:
                d.xyzDim[3] = dcmStrInt(lLength, &buffer[lPos]);
                break;
            case 	kSamplesPerPixel:
                d.samplesPerPixel = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kDim2:
                d.xyzDim[2] = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kDim1:
                d.xyzDim[1] = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kXYSpacing:
                dcmMultiFloat(lLength, (char*)&buffer[lPos], 2, d.xyzMM);
                
                break;
            case 	kImageComments:
                dcmStr (lLength, &buffer[lPos], d.imageComments);
                break;
            case 	kLocationsInAcquisitionGE:
                locationsInAcquisitionGE = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kBitsAllocated :
                d.bitsAllocated = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kBitsStored :
                d.bitsStored = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kIsSigned : //http://dicomiseasy.blogspot.com/2012/08/chapter-12-pixel-data.html
                d.isSigned = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kTR :
                d.TR = dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case 	kTE :
                d.TE = dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case kEchoNum :
                echoNum =  dcmStrInt(lLength, &buffer[lPos]);
                break;
            case 	kZSpacing :
                zSpacing = dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case 	kSlope :
                d.intenScale = dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case 	kIntercept :
                d.intenIntercept = dcmStrFloat(lLength, &buffer[lPos]);
                break;
            case 	kZThick :
                d.xyzMM[3] = dcmStrFloat(lLength, &buffer[lPos]);
                break;
                
            case 	kCoilSiemens : {
                if (d.manufacturer == kMANUFACTURER_SIEMENS) {
                    //see if image from single coil "H12" or an array "HEA;HEP"
                    char coilStr[kDICOMStr];
                    dcmStr (lLength, &buffer[lPos], coilStr);
                    //long coilNum = 0;
                    char *ptr;
                    dcmStrDigitsOnly(coilStr);
                    coilNum = strtol(coilStr, &ptr, 10);
                    if (*ptr != '\0')
                        coilNum = 0;
                }
                break; }
            case 	kLocationsInAcquisition :
                d.locationsInAcquisition = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case 	kLocationsInAcquisitionPhilips:
                locationsInAcquisitionPhilips = dcmInt(lLength,&buffer[lPos],d.isLittleEndian);
                break;
            case kIconImageSequence:
                isIconImageSequence = true;
                break;
            case 	kNumberOfDynamicScans:
                d.numberOfDynamicScans =  dcmStrInt(lLength, &buffer[lPos]);
                break;
            case	kMRAcquisitionType:
                if (lLength > 1) d.is3DAcq = (buffer[lPos]=='3') && (toupper(buffer[lPos+1]) == 'D');
                break;
            case	kMRAcquisitionTypePhilips: //kMRAcquisitionType
                if (lLength > 1) d.is3DAcq = (buffer[lPos]=='3') && (toupper(buffer[lPos+1]) == 'D');
                break;
            case	kAngulationRL:
                d.angulation[1] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kAngulationAP:
                d.angulation[2] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kAngulationFH:
                d.angulation[3] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kMRStackOffcentreRL:
                d.stackOffcentre[1] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kMRStackOffcentreAP:
                d.stackOffcentre[2] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kMRStackOffcentreFH:
                d.stackOffcentre[3] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kSliceOrient: {
                char orientStr[kDICOMStr];
                dcmStr (lLength, &buffer[lPos], orientStr);
                if (toupper(orientStr[0])== 'S')
                    d.sliceOrient = kSliceOrientSag; //sagittal
                else if (toupper(orientStr[0])== 'C')
                    d.sliceOrient = kSliceOrientCor; //coronal
                else
                    d.sliceOrient = kSliceOrientTra; //transverse (axial)
                break; }
            case	kDiffusionBFactor:
                if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (isAtFirstPatientPosition)) {
                    d.CSA.numDti++; //increment with BFactor: on Philips slices with B=0 have B-factor but no diffusion directions
                    if ((d.CSA.numDti > 0) && (d.CSA.numDti <= kMaxDTIv))
                        d.CSA.dtiV[d.CSA.numDti-1][0] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                }
                break;
            case    kDiffusionDirectionRL:
                if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (isAtFirstPatientPosition) && (d.CSA.numDti > 0) && (d.CSA.numDti <= kMaxDTIv))
                    d.CSA.dtiV[d.CSA.numDti-1][1] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case kDiffusionDirectionAP:
                if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (isAtFirstPatientPosition) && (d.CSA.numDti > 0) && (d.CSA.numDti <= kMaxDTIv))
                    d.CSA.dtiV[d.CSA.numDti-1][2] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                break;
            case	kDiffusionDirectionFH:
                if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (isAtFirstPatientPosition) && (d.CSA.numDti > 0) && (d.CSA.numDti <= kMaxDTIv))
                    d.CSA.dtiV[d.CSA.numDti-1][3] = dcmFloat(lLength, &buffer[lPos],d.isLittleEndian);
                //http://www.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:DICOM_for_DWI_and_DTI
                break;
            case 	kWaveformSq:
                d.imageStart = 1; //abort!!!
#ifdef myUseCOut
                std::cout<<"Warning: Unable to extract sound wave forms" <<std::endl;
#else
                printf("Warning: Unable to extract sound wave forms\n");
#endif
                break;
            case 	kCSAImageHeaderInfo:
                readCSAImageHeader(&buffer[lPos], lLength, &d.CSA, isVerbose);
                break;
                //case kObjectGraphics:
                //    printf("---->%d,",lLength);
                //    break;
            case kDiffusionBFactorGE :
                if (d.manufacturer == kMANUFACTURER_GE) d.CSA.dtiV[0][0] =  dcmStrInt(lLength, &buffer[lPos]);
                break;
            case kGeiisFlag:
                if ((lLength > 4) && (buffer[lPos]=='G') && (buffer[lPos+1]=='E') && (buffer[lPos+2]=='I')  && (buffer[lPos+3]=='I')) {
                    //read a few digits, as bug is specific to GEIIS, while GEMS are fine
#ifdef myUseCOut
                    std::cout<<"Warning: GEIIS violates the DICOM standard. Inspect results and admonish your vendor." <<std::endl;
#else
                    printf("Warning: GEIIS violates the DICOM standard. Inspect results and admonish your vendor.\n");
#endif
                    isIconImageSequence = true;
                    //geiisBug = true; //compressed thumbnails do not follow transfer syntax! GE should not re-use pulbic tags for these proprietary images http://sonca.kasshin.net/gdcm/Doc/GE_ImageThumbnails
                    
                }
                break;
            case 	kOrientation :
                dcmMultiFloat(lLength, (char*)&buffer[lPos], 6, d.orient);
                break;
            case 	kImageStart:
                //if ((!geiisBug) && (!isIconImageSequence)) //do not exit for proprietary thumbnails
                if (!isIconImageSequence) //do not exit for proprietary thumbnails
                    d.imageStart = (int)lPos;
                //geiisBug = false;
                isIconImageSequence = false;
                break;
        } //switch/case for groupElement
        //printf("VR=%c%c tag=%04x,%04x length=%lu, pos=%ld x=%d\n",vr[0],vr[1],groupElement & 65535,groupElement>>16, lLength, lPos, d.xyzDim[1]);
        lPos = lPos + (lLength);
    }
    free (buffer);
    d.dateTime = (atof(d.studyDate)* 1000000) + atof(d.studyTime);
    //printf("slices in Acq %d %d\n",d.locationsInAcquisition,locationsInAcquisitionPhilips);
    if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (d.locationsInAcquisition == 0))
        d.locationsInAcquisition = locationsInAcquisitionPhilips;
        if ((d.manufacturer == kMANUFACTURER_GE) && (d.locationsInAcquisition == 0))
            d.locationsInAcquisition = locationsInAcquisitionGE;
            if (zSpacing > 0) d.xyzMM[3] = zSpacing; //use zSpacing if provided: depending on vendor, kZThick may or may not include a slice gap
                //printf("patientPositions = %d XYZT = %d slicePerVol = %d numberOfDynamicScans %d\n",patientPositionCount,d.xyzDim[3], d.locationsInAcquisition, d.numberOfDynamicScans);
                if ((d.manufacturer == kMANUFACTURER_PHILIPS) && (patientPositionCount > d.xyzDim[3]))
                    printf("Please check slice thicknesses: Philips R3.2.2 bug can disrupt estimation (%d positions reported for %d slices)\n",patientPositionCount, d.xyzDim[3]); //Philips reported different positions for each slice!
                    if ((d.imageStart > 144) && (d.xyzDim[1] > 1) && (d.xyzDim[2] > 1))
                        d.isValid = true;
                        if ((d.xyzMM[1] > FLT_EPSILON) && (d.xyzMM[2] < FLT_EPSILON)) {
                            printf("Please check voxel size\n");
                            d.xyzMM[2] = d.xyzMM[1];
                        }
    if ((d.xyzMM[2] > FLT_EPSILON) && (d.xyzMM[1] < FLT_EPSILON)) {
        printf("Please check voxel size\n");
        d.xyzMM[1] = d.xyzMM[2];
    }
    
    if ((d.xyzMM[3] < FLT_EPSILON)) {
        printf("Unable to determine slice thickness: please check voxel size\n");
        d.xyzMM[3] = 1.0;
    }
    if (coilNum > 0) //segment images with multiple coils
        d.seriesNum = d.seriesNum + (100*coilNum);
        if (echoNum > 2) //segment images with multiple echoes
            d.seriesNum = d.seriesNum + (100*echoNum);
            if (isVerbose) {
                printf("Patient Position %f %f %f\n",d.patientPosition[1],d.patientPosition[2],d.patientPosition[3]);
                printf("DICOM acq %d img %d ser %ld dim %dx%dx%d mm %gx%gx%g offset %d dyn %d loc %d valid %d ph %d mag %d posReps %d nDTI %d 3d %d\n",d.acquNum,d.imageNum,d.seriesNum,d.xyzDim[1],d.xyzDim[2],d.xyzDim[3],d.xyzMM[1],d.xyzMM[2],d.xyzMM[3],d.imageStart, d.numberOfDynamicScans, d.locationsInAcquisition, d.isValid, d.isHasPhase, d.isHasMagnitude,d.patientPositionSequentialRepeats, d.CSA.numDti, d.is3DAcq);
            }
    return d;
} // readDICOM()

struct TDICOMdata readDICOM(char * fname) {
    return readDICOMv(fname, false);
}


