## PLUG-IN LIBRARY

The plug-in file is a dynamic library file (a so file on the Linux system, a dylib file on the Mac OSX system and a dll file in Microsoft Windows®) that implements specific functions called internally by the Engine at specific times during the pipeline. This is a major advantage of the FRIEND Engine, because, when in need of additional features, users can focus on writing just the necessary functionality for the specific needs of their own research. This allows customization of the neurofeedback tool, using encapsulated codes that run additional functionalities from external libraries, leaving the Engine code intact. This characteristic favors usability and code maintenance, because errors in a plug-in library are also encapsulated in that library and do not affect other plug-ins. This framework makes it easier to setup pilot experiments and to explore new hypotheses. We are also planning on maintaining a repository on the Internet for sharing plug-ins and front-ends to facilitate the spread of new technologies.

The plug-in library must implement all the processing needed to calculate the feedback response. The Engine defines six functions that are called at predefined time points during the pipeline execution. Not all of those six functions must be implemented in a plug-in library, just the ones necessary to properly calculate the feedback information. We recommend advanced users to code those functions in C/C++ because that avoids compatibility errors during the execution of the plug-in's functions by the engine. The main parameter of these functions is an object of the studyParameters class type, defined in the file vardb.cpp (FES/FES/PlugIn/vardb.cpp), which encloses all the information used by the processing pipeline to correctly identify the files, directories and expected number of volumes. Another important parameter passed to the plug-in is a pointer referencing the dynamically allocated memory created by inside the initialization function of the plug-in to be used for temporary calculations. The list of functions a plug-in can define are:

### 1. TRAIN

Executed when the front-end issues the TRAIN command. Normally, this function is used to analyze a complete data set (e.g. a training run). This function is generally quite time demanding, and is therefore not a "real-time" operation. It is normally called at the end of the acquisition of the training run. This function builds a model based on the training data that can be used on a subsequent classification run.

### 2. TEST

Executed when the Engine needs to calculate the feedback value for a given volume. There are two predefined values that the Engine returns to the front-end: information about the condition of the processed image volume indicating a specific type of image volume classification (useful for multi voxel patter analysis [MVPA]), and the feedback information value. The interpretation of these two variables is left to the front-end, which must implement how this information will be passed to the participant.

### 3. INITIALIZATION

Called after the engine reads information from the study configuration file (see Sato et al., 2013). All the memory data structure used by the plug-in must be initialized here. A pointer reference for this data structures has to be returned in the function argument. This reference will be used in further plug-in function calls issued by the FRIEND Engine.

### 4. FINALIZATION

Called right before ending the processing of a session. All the memory data structure allocated in the Initialization function must be destroyed here.

### 5. VOLUME

Called before the pre-processing of each volume.

### 6. POST PREPROCESSING

Called after the pre-processing of each volume.

### 7. PluginHandler class

The class PluginHandler, located in FES/FES/FRIEND\_Engine/PlugInHandler.h is responsible for handling the interface between engine and the plug-ins. The function definitions are also located in this header file and listed below :

```c

// plug in function pointer definitions
typedef int(*InitializationFunction)(studyParams&vardb, void *&userData);
typedef int(*TrainFunction)(studyParams&vardb, void *&userData);
typedef int(*TestFunction)(studyParams&vardb, int index, float &classnum, float &projection, void *&userData);
typedef int(*FinalizationFunction)(studyParams&vardb, void *&userData);
typedef int(*VolumeFunction)(studyParams&vardb, int index, char *fileName, void *&userData);
typedef int(*AfterPreprocessingFunction)(studyParams&vardb, int index, char *fileName, void *&userData);

```


<a name="available"></a>
## AVAILABLE PLUG-INS

The FRIEND Engine distribution comes with four plug-ins: one for the SVM pipeline (libBrainDecoding), using the libSVM library (Chang &Lin, 2011); one for the ROI pipeline (libROI), used in the Matlab and the Medieval frontend examples; one for the functional connectivity between two ROIs (libConnectivity) and one (libMotor) that extracts ROI information from two ROIs located in the motor cortex area (left and right), used in the character finger tapping virtual scenario. The source code of these plug-ins are located in the directory FES/FES/PlugIn.

In table 2 we list the functions implemented in the libROI plug-in. We recommend the users that want or need to write their plug-in to start studying the files **ROIPlugIn.cpp** and **ROIPlugIn.h** , in FES/FES/PlugIn/ROIbasePlugIn directory


> Table 2. Functions implemented in the libROI plug-in.

| Plug-in | LibROI |
| --- | --- |
| initializeROIProcessing | Initialization function that creates the data structures necessary for processing ROIs and transforms the ROI MNI mask to native space. |
| processROI | Feedback function that calculates the percent signal change of a ROI in the current volume compared to the mean activation of a baseline block. |
| finalizeROIProcessing | Finalization function that destroys all data structures created in the initialization function. |

#### Code snippets from the LibROI plugin

Here we present code snippets from the libROI plug-in used in the Matlab frontend and in the Medieval virtual scenario frontends of the distribution. The code snippet 1 contains the function processROI, which calls a class method that actually performs the calculations. It also demonstrates the use of some methods of the studyParameters class and volume manipulation with functions used within the FSL source code. Other functions and objects, like the one that calculates the mean of a ROI are also shown.

The targetValue variable used at the end of the method processVolume is read within the initializeROIProcessing function, presented in code snippet 2. The studyParams variable grants the plug-in library access to the setup configuration information of the study contained in the studyparams.txt file. Users can create their own configuration information in the configuration file and read it in the plug-in initialization function.

**Code snippet 1**. Main function of the libROI plug-in that calculates the feedback value.

```c

// libROI Plug-in function that calculates the feedback value that the engine will pass to the frontend
int roiProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &projection)
{
   char processedFile[200];
   int idxInterval = vdb.interval.returnInterval(index);

   volume <float> v ;

   // gets the motion corrected and Gaussian smoothed file
   vdb.getMCGVolumeName(processedFile , index);
   read_volume(v , string(processedFile));

   // if in baseline condition, adds current image to the previous sum image
   classnum = vdb.getClass(index);
   projection = 0;

    if ( vdb. interval. isBaselineCondition ( index ))
    {

      if ( vdb. interval. intervals [idxInterval]. start == index )
          meanbaseline = v;
      else meanbaseline += v;

      // when finishing adding images on baseline condition block,
      // divide the sum image by the size of the block to create the
      // mean volume
      if ( vdb.interval.intervals[idxInterval].end == index )
      {
         meanbaseline = ( vdb.interval.intervals[idxInterval].end 
                       - vdb.interval.intervals[idxInterval].start + 1);

         // obtain ROI mean value. 
         // The MeanCalculation object variable was initialized previously with the 
         // ROI mask in Initialization function.
         //
         // calculateMeans method obtains the mean for each ROI in the current volume
         meanCalculation.calculateMeans(meanbaseline);

         // Calculates the mean of the ROI for the current mean volume.
         lastBaselineValue = meanCalculation.means[0];
      }
    }
    // task condition. Taking the mean of the volume and calculating the PSC
    else {
      meanCalculation.calculateMeans(v);

      // Percent signal change calculation
      projection = PSC(meanCalculation.means[0], lastBaselineValue);

      // Divide feedback value by user-defined target value. The
      // front-end will use the feedback value for neurofeedback
      // display
      projection = projection / targetValue ;
      fprintf(stderr ,"Projection value = %f\n" , projection);
    }
    return 0;
}

// calculate the percent signal change value
float roiProcessing::PSC(float value, float base)
{
    return ( base ) ? (( value - base ) / base): 0;
}

// plugin function for calculationg feedback value
DLLExport processROI(studyParams &vdb, int index, float &classnum, float &projection, void *&userData)
{
   roiProcessing *roiVar = ( roiProcessing *) userData;
   roiVar->processVolume(vdb , index , classnum , projection);
   return 0;
}

```


**Code snippet 2**. Initialization function of the libROI plug-in.
```c

// initializes the object variables. This function brings the mni mask to subject space
int roiProcessing::initialization(studyParams &vdb)
{
    // verifies if the configuration file was read
    if(vdb.readedIni.IsEmpty())
   {
       fprintf(stderr ,"the studyparams.txt file was not read.\n");
       return 1;
   }

    // reads the specific plugin information. readedIni object contains all the information from the
   // configuration file
   strcpy(vdb.mniMask , vdb.readedIni.GetValue("FRIEND" ,"ActivationLevelMask"));
   strcpy(vdb.mniTemplate , vdb.readedIni.GetValue("FRIEND" ,"ActivationLevelMaskReference"));
   targetValue = vdb.readedIni.GetDoubleValue("FRIEND" ,"ActivationLevel");
   int masktype = vdb.readedIni.GetLongValue("FRIEND" ,"ActivationLevelMaskType");
   
   fprintf(stderr ,"mnimask = %s\n" , vdb.mniMask);
   fprintf(stderr ,"mnitemp = %s\n" , vdb.mniTemplate);
   fprintf(stderr ,"masktype = %d\n" , masktype);
   
    if (( fileExists ( vdb.mniMask )) && ( fileExists ( vdb.mniTemplate )) && ( masktype == 2 ))
    {
      char outputFile[500], prefix[500]="_RFI2" , name[500];
      
      extractFileName(vdb.mniMask , name);
      for (int t = 0; t < strlen( name ); t++)
      if (name[t] == '.') name[t] = '_';
      
      sprintf(outputFile ,"%s%s%s.nii" , vdb.inputDir , name , vdb.trainFeatureSuffix);
      
      fprintf(stderr ,"Calculating the native template %s\n" , outputFile);
      // bringing mni mask to subject space
      MniToSubject(vdb.rfiFile , vdb.mniMask , vdb.mniTemplate , outputFile , prefix);

      // loads the reference mask
      meanCalculation.loadReference(outputFile);
    }
    else if((fileExists(vdb.mniMask))&&(masktype == 1))
    {
          // loads the reference mask
          fprintf(stderr ,"Loading native space mask %s\n" , vdb.mniMask);
          meanCalculation.loadReference(vdb.mniMask);
    }
   lastBaselineValue = 0 ;
   return 0 ;
}

// plugin function for initializating the roi processing object
DLLExport initializeROIProcessing(studyParams &vdb, void *&userData)
{
   roiProcessing *roiVar = new roiProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

```

The following code snippets are from the Matlab frontend showing how information must be exchanged with the engine via a non-blocked TCP/IP communication in the ROI processing pipeline. The first command issued is the "NEWSESSION" command that creates a new session in memory and returns the session id that uniquely identifies the newly created session. As already explained in the main paper, a session is an independent location in the memory of the computer running the engine, capable of storing all the information needed to be sent back to the frontend, such as neurofeedback information and motion corrected volume parameters. We use the Matlab functions fprintf and fgetl to send and receive, respectively, messages via TCP/IP communication.

The variables mainThread and responseThread, presented in the following code snippets, are two Matlab TCP/IP objects. The connection established through the mainThread variable is the first made and lasts until the end of the acquisition run processing. The principal commands are sent through that connection.  The connection established through the responseThread is temporary. Each time an information is needed, such the feedback information for a specific volume, the connection is opened and after the information is acquired, the connection is closed.

**Code snippet 3**. Creating a new session in the engine

```matlab

fprintf(mainThread ,'NEWSESSION');

% reading session id
sessionID = fgetl(mainThread);

% reading acknowledge
response = fgetl(mainThread);

```

The next command issued by the frontend is the "PLUGIN", which sends the plug-in library filename and the names of the functions that the engine must call to accomplish the objective of the neurofeedback study. This list of function names must be sent in a pre-specified order (train, test, initialization, finalization, volume and post-processing as indicated by the Matlab % commentary directives to the right of each programming line). Not implemented functions must be specified as "no" in the place where the name of the function should be present.


**Code snippet 4**. Configuring the plugin to be used in the experiment

```matlab

% sending the PLUG-IN command and parameters
fprintf(mainThread, 'PLUGIN');
fprintf(mainThread, 'libROI');
fprintf(mainThread, 'no'); % train function
fprintf(mainThread, 'processROI'); % test function
fprintf(mainThread, 'initializeROIProcessing'); % initialization function
fprintf(mainThread, 'finalizeProcessing'); %finalization function
fprintf(mainThread, 'no'); % volume function
fprintf(mainThread, 'no'); % post preprocessing function

% getting the acknowledge
response = fgetl(mainThread);

```


The next command is "NBPREPROC", which initiates the preprocessing steps of the FRIEND pipeline in asynchronous (non-blocked) mode. Note that after this command, the front-end has to regularly query the Engine for the termination of this step.

**Code snippet 5.** Starting the PREPROC asynchronously

```matlab

% sending PREPROC non-blocked command
fprintf(mainThread, 'NBPREPROC');

% getting the acknowledge
response = fgetl(mainThread);

```

**Code snippet 6**. Getting a feedback value

```matlab

% open a communication channel with the engine
fopen(responseThread);

% sending the session command to create a new workspace
fprintf(responseThread, 'SESSION');

% sending the session id
fprintf(responseThread, '%s', sessionID);

response = fgetl(responseThread);

% sending the TEST sub-command
fprintf(responseThread, 'TEST');

% sending the volume index of the feedback value
fprintf(responseThread, '%d', actualVolume);

% getting feedback information
class = str2double(fgetl(responseThread));
percentage = str2double(fgetl(responseThread));

% getting acknowledge
response = fgetl(responseThread);

% closing the connection
fclose(responseThread);

```


Next, the "NBFEEDBACK" command is sent, initiating the processing of the volumes of the actual run. The front-end client issues various SESSION/GRAPHPARS commands through the responseThread variable to query for translation and rotation parameters of each motion corrected volume to be presented in graphs placed in the interface, and SESSION/TEST commands to control the thermometer to be presented to the participant, as showed in code snippet 5. If all the volume files in the run were processed, the SESSION/GRAPHPARS command returns an "END" token, indicating the termination of the feedback process.


**Code snippet 7**. Example of how to use the SET command

```matlab

% changing the mask type to native space
fprintf(mainThread ,'SET');
fprintf(mainThread ,'ActivationLevelMaskType');
fprintf(mainThread ,'1');
response = fgetl(mainThread);

% changing the mask to the generated by the funcional localizer run
fprintf(mainThread ,'SET');
fprintf(mainThread ,'ActivationLevelMask');
fprintf(mainThread ,'glmdirtstats_features_RUN01\_bin');
response = fgetl(mainThread);

% changing the mask to the generated by the funcional localizer run
fprintf(mainThread ,'SET');
fprintf(mainThread ,'Prefix');
fprintf(mainThread ,'outputdirRUN02/DRIN-');
response = fgetl(mainThread);

```

> **Code snippet 7** show how to use the SET command, one of the commands to change pieces of the read in memory study\_params.txt. In this snippet, the Matlab frontend changes the ActivationLevelMaskType to 1, indicating that the mask informed in the variable ActivationLevelMask is in native space. It also changes this mask, to the binary version of the result of the FEATURESELECTION command. The last SET command changes the Prefix variable, indicating that the volumes of the actual RUN should be find in RUN02 directory. The _glmdir_ and _outputdir_ are replaced in the engine for the actual values. _glmdir_ is replaced by the glm output directory, inside the subject directory  (in the engine distribution is FES/FES/TestData/SUBJ002/glm). The _outputdir_ is replaced by the subject directory(in the engine distribution is FES/FES/TestData/SUBJ002).