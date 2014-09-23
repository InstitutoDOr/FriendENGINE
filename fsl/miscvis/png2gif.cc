#include <stdio.h>
#include "libvis/miscplot.h"

int main()
{
  gdImagePtr im;

  FILE  *inpng  = fopen("hello.png", "rb" );
  FILE  *outgif = fopen("hello.gif", "wb" );
  im = gdImageCreateFromPng(inpng);
  gdImageGif(im, outgif); 
}
