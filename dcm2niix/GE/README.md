## About

dcm2niix attempts to convert GE DICOM format images to NIfTI. The current generation DICOM files generated be GE equipment is quite impoverished relative to other vendors. Therefore, the amount of information dcm2niix is able to extract is relatively limited. Hopefully, in the future GE will provide more details that are critical for brain scientists.

## Diffusion Tensor Notes

The [NA-MIC Wiki](https://www.na-mic.org/wiki/NAMIC_Wiki:DTI:DICOM_for_DWI_and_DTI#Private_vendor:_GE) provides a nice description of the GE diffusion tags. In brief, the B-value is stored as the first element in the array of 0043,1039. The DICOM elements 0019,10bb, 0019,10bc and 0019,10bd provide the gradient direction relative to the frequency, phase and slice. As noted by Jaemin Shin (GE), the GE convention for reported diffusion gradient direction has always been in “MR physics” logical coordinate, i.e Freq (X), Phase (Y), Slice (Z). Note that this is neither “with reference to the scanner bore” (like Siemens or Philips) nor “with reference to the imaging plane” (as expected by FSL tools). This is the main source of confusion. This explains why the dcm2niix function geCorrectBvecs() checks whether the DICOM tag In-plane Phase Encoding Direction (0018,1312) is 'ROW' or 'COL'. In addition, it will generate the warning 'reorienting for ROW phase-encoding untested' if you acquire DTI data with the phase encoding in the ROW direction. If you do test this feature, please report your findings as a Github issue. Assuming you have COL phase encoding, dcm2niix should provide [FSL format](http://justinblaber.org/brief-introduction-to-dwmri/) [bvec files](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/FAQ#What_conventions_do_the_bvecs_use.3F).

## Slice Timing

Knowing the relative timing of the acquisition for each 2D slice in a 3D volume is useful for [slice time correction](https://www.mccauslandcenter.sc.edu/crnl/tools/stc) of both fMRI and DTI data. Unfortunately, current GE software does not provide a consistent way to record this.

[Some sequences](https://afni.nimh.nih.gov/afni/community/board/read.php?1,154006) encode the RTIA Timer (0021,105E) element. For example, [this dataset DV24](https://github.com/nikadon/cc-dcm2bids-wrapper/tree/master/dicom-qa-examples/ge-mr750-slice-timing) includes timing data, while [this DV26 dataset does not](https://github.com/neurolabusc/dcm_qa_nih). Even with the sequences that do encode the RTIA Timer, there is some debate regarding the accuracy of this element. In the example listed, the slice times are clearly wrong in the first volume. Therefore, dcm2niix always estimates slice times based on the 2nd volume in a time series.

In general, fMRI acquired using GE product sequence (PSD) “epi” with the multiphase option will store slice timing in the Trigger Time (DICOM 0018,1060) element. The current version of dcm2niix ignores this field, as no examples are available. In contrast, the popular PSD “epiRT” (BrainWave RT, fMRI/DTI package provided by Medical Numerics) does not save this tag (though in some cases it saves the RTIA Timer). Examples are [available](https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage#Slice_timing_correction) for both the “epiRT” and “epi” sequences.

## User Define Data GE  (0043,102A)

This private element of the DICOM header is used to determine the phase encoding polarity. Specifically, we need to know the "Ky traversal direction" (top-down, or bottom up) and the phase encoding polarity. Unfortunately, this data is stored in a complicated, proprietary structure, that has changed with different releases of GE software. [Click here to see the definition for this structure](https://github.com/ScottHaileRobertson/GE-MRI-Tools/blob/master/GePackage/%2BGE/%2BPfile/%2BHeader/%2BRDB15/rdbm.h).

## Total Readout Time

One often wants to determine [echo spacing, bandwidth, ](https://support.brainvoyager.com/brainvoyager/functional-analysis-preparation/29-pre-processing/78-epi-distortion-correction-echo-spacing-and-bandwidth) and total read-out time for EPI data so they can be undistorted. Total readout time is influence by parallel acceleration factor, bandwidth, number of EPI lines, and partial Fourier. Not all of these parameters are available from the GE DICOM images, so a user needs to check the scanner console.

## GE Protocol Data Block

In addition to the public DICOM tags, previous versions of dcm2niix attempted to decode the proprietary GE Protocol Data Block (0025,101B). This is essentially a [GZip format](http://www.onicos.com/staff/iz/formats/gzip.html) file embedded inside the DICOM header. Unfortunately, this data seems to be [unreliable](https://github.com/rordenlab/dcm2niix/issues/163) and therefore this strategy is not used anymore. The notes below regarding the usage of this data block are provided for historical purposes.

 - The VIEWORDER tag is used to set the polarity of the BIDS tag PhaseEncodingDirection, with VIEWORDER of 1 suggesting bottom up phase encoding. Unfortunately, users can separately reverse the phase encoding direction making this tag unreliable.
 - The SLICEORDER tag could be used to set the SliceTiming for the BIDS tag PhaseEncodingDirection, with a SLICEORDER of 1 suggesting interleaved acquisition.
 - There are reports that newer versions of GE equipement (e.g. DISCOVERY MR750 / 24\MX\MR Software release:DV24.0_R01_1344.a) are now storing an [XML](https://groups.google.com/forum/#!msg/comp.protocols.dicom/mxnCkv8A-i4/W_uc6SxLwHQJ) file within the Protocolo Data Block (compressed). In theory this might also provide useful information.

## Sample Datasets

 - [A validation dataset for dcm2niix commits](https://github.com/neurolabusc/dcm_qa_nih).
 - [Examples of phase encoding polarity, slice timing and diffusion gradients](https://github.com/nikadon/cc-dcm2bids-wrapper/tree/master/dicom-qa-examples/).
 - The dcm2niix (wiki)[https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage] includes examples of diffusion data, slice timing, and other variations.