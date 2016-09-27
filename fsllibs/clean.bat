echo Removing intermediate files and directories of basisfield
rmdir basisfield\release /s /q
rmdir basisfield\debug /s /q
rmdir basisfield\basisfield\Release /s /q
rmdir basisfield\basisfield\debug /s /q
del   basisfield\basisfield.sdf
del   basisfield\basisfield.v12.suo

echo Removing intermediate files and directories of commandsLib
rmdir commandsLib\release /s /q
rmdir commandsLib\debug /s /q
rmdir commandsLib\commandsLib\Release /s /q
rmdir commandsLib\commandsLib\debug /s /q
del   commandsLib\commandsLib.sdf
del   commandsLib\commandsLib.v12.suo

echo Removing intermediate files and directories of buildlibs
rmdir buildlibs\release /s /q
rmdir buildlibs\debug /s /q
rmdir buildlibs\buildlibs\Release /s /q
rmdir buildlibs\buildlibs\debug /s /q
del   buildlibs\buildlibs.sdf
del   buildlibs\buildlibs.v12.suo

echo Removing intermediate files and directories of cprob
rmdir cprob\release /s /q
rmdir cprob\debug /s /q
rmdir cprob\cprob\Release /s /q
rmdir cprob\cprob\debug /s /q
del   cprob\cprob.sdf
del   cprob\cprob.v12.suo

echo Removing intermediate files and directories of fslio
rmdir fslio\release /s /q
rmdir fslio\debug /s /q
rmdir fslio\fslio\Release /s /q
rmdir fslio\fslio\debug /s /q
del   fslio\fslio.sdf
del   fslio\fslio.v12.suo

echo Removing intermediate files and directories of meshclass
rmdir meshclass\release /s /q
rmdir meshclass\debug /s /q
rmdir meshclass\meshclass\Release /s /q
rmdir meshclass\meshclass\debug /s /q
del   meshclass\meshclass.sdf
del   meshclass\meshclass.v12.suo

echo Removing intermediate files and directories of miscmaths
rmdir miscmaths\release /s /q
rmdir miscmaths\debug /s /q
rmdir miscmaths\miscmaths\Release /s /q
rmdir miscmaths\miscmaths\debug /s /q
del   miscmaths\miscmaths.sdf
del   miscmaths\miscmaths.v12.suo

echo Removing intermediate files and directories of newimage
rmdir newimage\release /s /q
rmdir newimage\debug /s /q
rmdir newimage\newimage\Release /s /q
rmdir newimage\newimage\debug /s /q
del   newimage\newimage.sdf
del   newimage\newimage.v12.suo

echo Removing intermediate files and directories of newmat
rmdir newmat\release /s /q
rmdir newmat\debug /s /q
rmdir newmat\newmat\Release /s /q
rmdir newmat\newmat\debug /s /q
del   newmat\newmat.sdf
del   newmat\newmat.v12.suo

echo Removing intermediate files and directories of niftiio
rmdir niftiio\release /s /q
rmdir niftiio\debug /s /q
rmdir niftiio\niftiio\Release /s /q
rmdir niftiio\niftiio\debug /s /q
del   niftiio\niftiio.sdf
del   niftiio\niftiio.v12.suo

echo Removing intermediate files and directories of utils
rmdir utils\release /s /q
rmdir utils\debug /s /q
rmdir utils\utils\Release /s /q
rmdir utils\utils\debug /s /q
del   utils\utils.sdf
del   utils\utils.v12.suo

echo Removing intermediate files and directories of warpfns
rmdir warpfns\release /s /q
rmdir warpfns\debug /s /q
rmdir warpfns\warpfns\Release /s /q
rmdir warpfns\warpfns\debug /s /q
del   warpfns\warpfns.sdf
del   warpfns\warpfns.v12.suo

echo Removing intermediate files and directories of zlib
rmdir zlib\release /s /q
rmdir zlib\debug /s /q
rmdir zlib\zlib\Release /s /q
rmdir zlib\zlib\debug /s /q
del   zlib\zlib.sdf
del   zlib\zlib.v12.suo

echo Removing intermediate files and directories of znzlib
rmdir znzlib\release /s /q
rmdir znzlib\debug /s /q
rmdir znzlib\znzlib\Release /s /q
rmdir znzlib\znzlib\debug /s /q
del   znzlib\znzlib.sdf
del   znzlib\znzlib.v12.suo

del Libs\*.* /s /q