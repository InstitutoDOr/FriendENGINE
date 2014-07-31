set graphpic [ image create photo -file [ file rootname $argv ].ppm ]
$graphpic write  [ file rootname $argv ].gif -format gif
exit
   