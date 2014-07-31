#include <iostream>
//#include "fmribmain.h"
#include "newimageall.h"
#include <time.h>

using namespace NEWIMAGE;
//using namespace GENERALIO;


//  template <class T>
//  int fmrib_main(int argc, char* argv[]) {
//    cout << "Inside fmrib_main" << endl;
//    cout << "Data type  : sizeof() = " << sizeof(T) << endl;
//    T val;
//    val = (T) -4.6;
//    cout << "Result of val = -4.6 is " << val << endl;
//    return 0;
//  }

//  const short& star(short *it) {


void timewaste(long n) {
  cerr << "Timewaste answer = ";
  double sum=0, sigma=6.3;
  for (double t=0; t<(double) n; t++) {
    sum += exp(-(t*t)/sigma);
  }
  cerr << sum << endl;
}

int dumbextrap(const volume<int>& vol, int x, int y, int z)
{
  return -202 + x + y + z;
}

float dumbinterp(const volume<int>& vol, float x, float y, float z)
{
  return (float) -vol((int)x,(int)y,(int)z);
}


int main(int argc, char *argv[]) 
{

  try {

    cout << "Test for min and max" << endl;
    volume<float> v1, v2, v3;
    read_volume(v1,argv[1]);
    read_volume(v2,argv[2]);
    print_volume_info(v1,argv[1]);
    cout << "Coordinates of min voxel are: ("<<v1.mincoordx()<<","<<v1.mincoordy()<<","<<v1.mincoordz()<<")"<<endl;
    cout << "Coordinates of max voxel are: ("<<v1.maxcoordx()<<","<<v1.maxcoordy()<<","<<v1.maxcoordz()<<")"<<endl;
    print_volume_info(v2,argv[2]);
    cout << "Coordinates of min voxel are: ("<<v2.mincoordx()<<","<<v2.mincoordy()<<","<<v2.mincoordz()<<")"<<endl;
    cout << "Coordinates of max voxel are: ("<<v2.maxcoordx()<<","<<v2.maxcoordy()<<","<<v2.maxcoordz()<<")"<<endl;
    v3 = min(v1,v2);
    print_volume_info(v3,"Min of v1 and v2");
    cout << "Coordinates of min voxel are: ("<<v3.mincoordx()<<","<<v3.mincoordy()<<","<<v3.mincoordz()<<")"<<endl;
    cout << "Coordinates of max voxel are: ("<<v3.maxcoordx()<<","<<v3.maxcoordy()<<","<<v3.maxcoordz()<<")"<<endl;
    v3 = max(v1,v2);
    print_volume_info(v3,"Max of v1 and v2");
    cout << "Coordinates of min voxel are: ("<<v3.mincoordx()<<","<<v3.mincoordy()<<","<<v3.mincoordz()<<")"<<endl;
    cout << "Coordinates of max voxel are: ("<<v3.maxcoordx()<<","<<v3.maxcoordy()<<","<<v3.maxcoordz()<<")"<<endl;
    return 0;



    cout << "Test for robust min and max" << endl;
    volume<float> vol;
    read_volume(vol,argv[1]);
    print_volume_info(vol,argv[1]);
    cout << "Robust range = " << vol.robustmin() << " " << vol.robustmax() << endl;
    return 0;




  // cout << "Inside main - about to call_fmrib_main" << endl;
  // return call_fmrib_main(DT_COMPLEX, argc, argv);

    int *test = new int[64*64*20*3];
    cerr << "First test " << test[0] << " " << test[1] << endl;
    cerr << "Here endth the test" << endl;

    volume4D<int> vol1(64,64,20,3);
    cerr << "TEST " << vol1(1,0,0,0) << endl;
    for (int t=0; t<3; t++) {
      for (int z=0; z<20; z++) {
	for (int y=0; y<64; y++) {
	  for (int x=0; x<64; x++) {
	    cerr << vol1(x,y,z,t) << "  ";
	  }
	}
      }
    }
    cerr << "Before save_volume4D" << endl;
//      vol1 = 1;
//      vol1 += vol1;
    save_volume4D(vol1, "vol1");

//   if (argc<6) { 
//     cerr << "Usage: " << argv[0] << " <input-file> <output-file> dim1 dim2 dim3" << endl; exit(1); 
//   }

//   volume<int> v1;
//   int dim1 = atoi(argv[3]), dim2 = atoi(argv[4]), dim3 = atoi(argv[5]);
//   v1.swapdimensions(dim1,dim2,dim3);
//   v1.defineuserinterpolation(dumbinterp);
//   v1.setinterpolationmethod(userinterpolation);
//   volume<int> vtr(2*v1.xsize(),2*v1.ysize(),2*v1.zsize());
//   vtr.setdims(v1.xdim(),v1.ydim(),v1.zdim());
//   affine_transform(v1,vtr,IdentityMatrix(4));

//   exit(0);

// //    read_volume(v1,argv[1]);
// //    cout << "Read Volume : (X,Y,Z) = (" << v1.xsize() << "," << v1.ysize()
// //         << "," << v1.zsize() << ") and (x,y,z) = (" << v1.getx() << ","
// //         << v1.gety() << "," << v1.getz() << ")" << endl;
// //    //print_volume_info(v1,argv[1]);
// //    cout << "COG is " << v1.cog().t() << endl;
// //    cout << "Background is " << v1.backgroundval() << endl;

// //    save_volume(v1,"testvol");

//   cout << endl << endl;
//   volume<short> v2(3,2,1);
//   v2.setpadvalue(1);
//   v2 = 0;
//   v2(0,0,0) = -2;
//   v2(1,0,0) = -3;
//   v2(2,0,0) = -4;
//   v2(0,1,0) = 2;
//   v2(1,1,0) = 3;
//   v2(2,1,0) = 4;


// //    save_volume(v2,"testvol");

//   volume<double> vd;
//   read_volume(vd,"/usr/people/mark/reg/stroke/ab");
//   save_volume(vd,"dbltestvol");

// //    volume4D<double> td;
// //    read_volume4D(td,"/usr/people/mark/motionc/data/clean");
// //    for (unsigned int t=0; t<td.tsize(); t++) {
// //      if ((t % 2)==0)  td[t] *= -1;
// //    }

//   timewaste(10000000);

//   volume4D<int> tnew;
//   volume<int> vsh, vsh1;
//   copyconvert(vd,vsh);
//   vsh.setROIlimits(10,50,14,210,200,17);  // also activates ROI
//   //vsh.activateROI();
//   vsh1.copyROIonly(vsh);
//   tnew.addvolume(vsh1);
//   tnew.addvolume(vsh1 * -2.0);
//   tnew.addvolume(vsh1);
//   tnew.deletevolume(2);

//   timewaste(10000000);

//   {
//     volume4D<double> dnew(100,100,30,100);
//     dnew = (double) vsh1(0,0,0);
//     timewaste(10000000);
//   }

//   timewaste(10000000);

//   volume<int> vmask(vsh1);
//   vmask = 0;
//   vmask(2,2,2) = 1;
//   vmask(2,2,3) = 1;
//   {
//     Matrix tmat = tnew.matrix();
//     cout << "MATRIX SIZE IS : " << endl << tmat.Ncols() << " by " 
// 	 << tmat.Nrows() << endl;
//     tmat(2,2)=57000;
//     tnew.setmatrix(tmat);
//     timewaste(10000000);
//   }

//   timewaste(10000000);

//   save_volume4D(tnew,"tstest");

//   return 0;


// //    cout << endl << "EXTRASLICE" << endl;

// //    v2.setextrapolationmethod(extraslice);
// //    for (int k=0; k<1; k++) {
// //      for (int j=-4; j<8; j++) {
// //        for (int i=-4; i<7; i++) {
// //  	cout << v2(i,j,k) << " , ";
// //        }
// //        cout << endl;
// //      }
// //      cout << endl << endl;
// //    }
  
// //    cout << endl << "PERIODIC" << endl;

// //    v2.setextrapolationmethod(periodic);
// //    for (int k=0; k<1; k++) {
// //      for (int j=-4; j<8; j++) {
// //        for (int i=-4; i<7; i++) {
// //  	cout << v2(i,j,k) << " , ";
// //        }
// //        cout << endl;
// //      }
// //      cout << endl << endl;
// //    }

// //    cout << endl << "MIRROR" << endl;

// //    v2.setextrapolationmethod(mirror);
// //    for (int k=0; k<1; k++) {
// //      for (int j=-4; j<8; j++) {
// //        for (int i=-4; i<7; i++) {
// //  	cout << v2(i,j,k) << " , ";
// //        }
// //        cout << endl;
// //      }
// //      cout << endl << endl;
// //    }

// //    cout << endl << "BOUNDSASSERT" << endl;

// //    v2.setextrapolationmethod(boundsassert);
// //    for (int k=0; k<1; k++) {
// //      for (int j=0; j<8; j++) {
// //        for (int i=0; i<7; i++) {
// //  	cout << v2(i,j,k) << " , ";
// //        }
// //        cout << endl;
// //      }
// //      cout << endl << endl;
// //    }

//   time_t t0, t1, t2, t3, t4, t5, t6, t7, t8;
//   float sum1=0.0, sum2=0.0, sum3=0.0, sum4=0.0, sum5=0.0, sum6=0.0,
//     sum7=0.0, sum8=0.0;

// //    MJIMAGE::volume v4(200,250,200);
// //    v4 = 2;

// //    for (int n=0; n<20; n++) 
// //    for (int k=0; k<v4.zsize(); k++) {
// //      for (int j=0; j<v4.ysize(); j++) {
// //        for (int i=0; i<v4.xsize(); i++) {
// //  	sum3+= v4(i,j,k);
// //        }
// //      }
// //    }

// //    float *v5 = new float[200*250*200];
// //    float *pend = v5 + 200*250*200;
// //    for (float *it=v5; it!=pend; ++it) {  *it=2.0; } 
// //    for (int n=0; n<20; n++) {
// //      for (float *it=v5; it!=pend; ++it) {  sum1+= *it; } 
// //    }

//   time(&t0);

// //    volume<float> v3(200,250,200);
// //    const volume<float>& v3c(v3);
// //    v3 = 2;

// //    for (int n=0; n<20; n++) 
// //    sum2 += v3.sum();

// //    for (int n=0; n<20; n++) 
// //    for (volume<float>::fast_const_iterator it=v3c.fbegin(), 
// //  	 pend=v3c.fend();   it!=pend; ++it) 
// //      {  sum1+= *it; }

// //    cout << endl << "ABOUT TO SET ROI" << endl;

// //    int rx=0, ry=0, rz=0;
//   //v3.setROIlimits(1,1,1,100,100,100);
//   //v3.activateROI();
// //    for (int n=0; n<20; n++) 
// //    for (volume<float>::iterator it=v3.begin(), 
// //  	 pend=v3.end();   it!=pend; ++it) 
// //      {  
// //        sum1+= *it; 
// //        //it.getposition(rx,ry,rz);
// //        // if (n==1) cout << "Position is ("<<rx<<","<<ry<<","<<rz<<")"<<endl;
// //      }

// //    volume<float>::fast_const_iterator pend = v3.fend();
// //    for (int n=0; n<20; n++) 
// //    for (volume<float>::fast_const_iterator it=v3.fbegin(); it.iter!=pend.iter; ++it) 
// //      {  sum1+= *it; }

//   time(&t1);
// //    sum2=v3.sum();

// //    float *pend = &v3(199,249,199);
// //    for (int n=0; n<20; n++) {
// //      for (float *it=&v3(0,0,0); it!=pend; ++it) {  sum1+= *it; } 
// //    }

//   time(&t2);
// //    for (int n=0; n<20; n++) 
// //    for (int k=v3.limits(2); k<=v3.limits(5); k++) {
// //      for (int j=v3.limits(1); j<=v3.limits(4); j++) {
// //        for (int i=v3.limits(0); i<=v3.limits(3); i++) {
// //  	sum3+= v3.value(i,j,k);
// //        }
// //      }
// //    }
// //    for (int n=0; n<20; n++) 
// //      for (int i=v3.ROIlimits()[0]; i<=v3.ROIlimits()[3]; i++) {
// //        for (int j=v3.ROIlimits()[1]; j<=v3.ROIlimits()[4]; j++) {
// //  	for (int k=v3.ROIlimits()[2]; k<=v3.ROIlimits()[5]; k++) {
// //  	  sum3+= v3.value(i,j,k);
// //  	}
// //        }
// //      }

//   time(&t3);
// //    for (int n=0; n<20; n++) 
// //    for (int k=v3.ROIlimits()[2]; k<=v3.ROIlimits()[5]; k++) {
// //      for (int j=v3.ROIlimits()[1]; j<=v3.ROIlimits()[4]; j++) {
// //        for (int i=v3.ROIlimits()[0]; i<=v3.ROIlimits()[3]; i++) {
// //  	sum4+= v3(i,j,k);
// //        }
// //      }
// //    }

//   time(&t4);
// //    for (volume<float>::fast_const_iterator it=v3.fbegin(); it<v3.fend();
// //    ++it) {  sum5+= *it; }

//   time(&t5);
// //    sum6=v3.sum();

//   time(&t6);
// //    for (int k=v3.ROIlimits()[2]; k<=v3.ROIlimits()[5]; k++) {
// //      for (int j=v3.ROIlimits()[1]; j<=v3.ROIlimits()[4]; j++) {
// //        for (int i=v3.ROIlimits()[0]; i<=v3.ROIlimits()[3]; i++) {
// //  	sum7+= v3.value(i,j,k);
// //        }
// //      }
// //    }

//   time(&t7);
// //    for (int k=v3.ROIlimits()[2]; k<=v3.ROIlimits()[5]; k++) {
// //      for (int j=v3.ROIlimits()[1]; j<=v3.ROIlimits()[4]; j++) {
// //        for (int i=v3.ROIlimits()[0]; i<=v3.ROIlimits()[3]; i++) {
// //  	sum8+= v3(i,j,k);
// //        }
// //      }
// //    }

//   time(&t8);

//   cout << "Timings are : " << t1-t0 << " " << t2-t1 << " " << t3-t2 
//        << " " << t4-t3 << " " << t5-t4 << " " << t6-t5 << " " << t7-t6 
//        << " " << t8-t7 << endl;
//   cout << "Sums are : " << sum1 << " " << sum2 << " " << sum3 << " "
//        << sum4 << " " << sum5 << " " << sum6 << " " << sum7 << " " << sum8
//        << endl;

//   cout << "Limits = ";
//   for (int n=0; n<6; n++) { cout << v2.ROIlimits()[n] << " , "; }
//   cout << endl;

  } 
  catch (string msg)
        { cerr << "Exception message::" << msg << endl; throw; }
    
}




