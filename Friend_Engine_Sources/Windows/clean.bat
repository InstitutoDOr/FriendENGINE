echo Removing intermediate files and directories of FriendEngine
rmdir FriendEngine\release /s /q
rmdir FriendEngine\debug /s /q
rmdir FriendEngine\x64 /s /q
rmdir FriendEngine\FriendEngine\Release /s /q
rmdir FriendEngine\FriendEngine\debug /s /q
rmdir FriendEngine\FriendEngine\x64 /s /q
del   FriendEngine\FriendEngine.sdf
del   FriendEngine\FriendEngine.v12.suo /AH

echo Removing intermediate files and directories of libROI
rmdir libroi\release /s /q
rmdir libroi\debug /s /q
rmdir libroi\x64 /s /q
rmdir libROI\libROI\Release /s /q
rmdir libROI\libROI\debug /s /q
rmdir libROI\libROI\x64 /s /q
del   libroi\libroi.sdf
del   libroi\libroi.suo
del   libroi\libroi.v12.suo /AH

echo Removing intermediate files and directories of libMotor
rmdir libmotor\release /s /q
rmdir libmotor\debug /s /q
rmdir libmotor\x64 /s /q
rmdir libmotor\libmotor\Release /s /q
rmdir libmotor\libmotor\debug /s /q
rmdir libmotor\libmotor\x64 /s /q
del   libmotor\libmotor.sdf
del   libmotor\libmotor.suo
del   libmotor\libmotor.v12.suo /AH

echo Removing intermediate files and directories of libEmotionRoi
rmdir libEmotionROI\release /s /q
rmdir libEmotionROI\debug /s /q
rmdir libEmotionROI\x64 /s /q
rmdir libEmotionROI\libEmotionROI\Release /s /q
rmdir libEmotionROI\libEmotionROI\debug /s /q
rmdir libEmotionROI\libEmotionROI\x64 /s /q
del   libEmotionROI\libEmotionROI.sdf
del   libEmotionROI\libEmotionROI.v12.suo /AH

echo Removing intermediate files and directories of libConnectivity
rmdir libConnectivity\release /s /q
rmdir libConnectivity\debug /s /q
rmdir libConnectivity\x64 /s /q
rmdir libConnectivity\libConnectivity\Release /s /q
rmdir libConnectivity\libConnectivity\debug /s /q
rmdir libConnectivity\libConnectivity\x64 /s /q
del   libConnectivity\libConnectivity.sdf
del   libConnectivity\libConnectivity.v12.suo /AH

echo Removing intermediate files and directories of libBrainDecoding
rmdir libBrainDecoding\release /s /q
rmdir libBrainDecoding\debug /s /q
rmdir libBrainDecoding\x64 /s /q
rmdir libBrainDecoding\libBrainDecoding\Release /s /q
rmdir libBrainDecoding\libBrainDecoding\debug /s /q
rmdir libBrainDecoding\libBrainDecoding\x64 /s /q
del   libBrainDecoding\libBrainDecoding.sdf
del   libBrainDecoding\libBrainDecoding.v12.suo /AH

echo Removing intermediate files and directories of libMDD
rmdir libMDD\release /s /q
rmdir libMDD\debug /s /q
rmdir libMDD\x64 /s /q
rmdir libMDD\libMDD\Release /s /q
rmdir libMDD\libMDD\debug /s /q
rmdir libMDD\libMDD\x64 /s /q
del   libMDD\libMDD.sdf
del   libMDD\libMDD.v12.suo /AH

echo Removing intermediate files and directories of libAttention
rmdir libAttention\release /s /q
rmdir libAttention\debug /s /q
rmdir libAttention\x64 /s /q
rmdir libAttention\libAttention\Release /s /q
rmdir libAttention\libAttention\debug /s /q
rmdir libAttention\libAttention\x64 /s /q
del   libAttention\libAttention.sdf
del   libAttention\libAttention.v12.suo /AH

echo Removing intermediate files and directories of buildAll
rmdir buildAll\release /s /q
rmdir buildAll\debug /s /q
rmdir buildAll\x64 /s /q
del   buildAll\buildAll.sdf
del   buildAll\buildAll.v12.suo /AH

echo Removing executable files
del Application\*.exe
del Application\*.ilk
del Application\*.dll
del Application\*.iobj
del Application\*.ipdb

echo Deleting the testdata intermediate files and directories
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\input /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\glm /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\log /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\svm /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\preproc /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\preproc_RUN01 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\preproc_RUN02 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\preproc_RUN03 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\preproc_RUN04 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN01 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN02 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN03 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN04 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\copied_RUN01 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\copied_RUN02 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\copied_RUN03 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\copied_RUN04 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN01 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN02 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN03 /s /q
rmdir ..\Friend_Engine_Sources\TestData\SUBJ002\activations_RUN04 /s /q

del ..\Friend_Engine_Sources\TestData\SUBJ002\MNI2RFI_RFI2.mat
del ..\Friend_Engine_Sources\TestData\SUBJ002\RFI2MNI_RFI2.mat
del ..\Friend_Engine_Sources\TestData\SUBJ002\RFI2MNI_RFI2.nii
del ..\Friend_Engine_Sources\TestData\SUBJ002\study_params_RUN01.txt
del ..\Friend_Engine_Sources\TestData\SUBJ002\study_params_RUN02.txt
del ..\Friend_Engine_Sources\TestData\SUBJ002\study_params_RUN03.txt
del ..\Friend_Engine_Sources\TestData\SUBJ002\study_params_RUN04.txt

del ..\Friend_Engine_Sources\TestData\SUBJ002\ROIsMap_RUN01.nii
del ..\Friend_Engine_Sources\TestData\SUBJ002\ROIsMap_RUN02.nii
del ..\Friend_Engine_Sources\TestData\SUBJ002\ROIsMap_RUN03.nii
del ..\Friend_Engine_Sources\TestData\SUBJ002\ROIsMap_RUN04.nii
del ..\\Friend_Engine_Sources\TestData\hmat_spm_final_std.nii


echo Removing fsl.dll files
rmdir ..\..\fsl\release /s /q
rmdir ..\..\fsl\debug /s /q
del   ..\..\fsl\fsl.sdf
del   ..\..\fsl\fsl.opensdf

echo Removing intermediate files and directories of gdlib
rmdir ..\..\fsl\gdlib\release /s /q
rmdir ..\..\fsl\gdlib\debug /s /q
rmdir ..\..\fsl\gdlib\x64 /s /q
rmdir ..\..\fsl\gdlib\gdlib\Release /s /q
rmdir ..\..\fsl\gdlib\gdlib\debug /s /q
rmdir ..\..\fsl\gdlib\gdlib\x64 /s /q
del   ..\..\fsl\gdlib\gdlib.sdf
del   ..\..\fsl\gdlib\gdlib.v12.suo /AH
