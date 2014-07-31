#if !defined (MISCPLOT_H)
#define MISCPLOT_H

#include <string>
#include "newmatap.h"
#include "miscmaths/histogram.h" 
#include "gd.h"

#ifndef FALSE
#define FALSE false
#endif
#ifndef TRUE
#define TRUE true
#endif


namespace MISCPLOT{

	extern unsigned long sc_init[64];
  class miscplot
  {
    public:
      
   		//constructor
      miscplot(){
				req_xsize=0; req_ysize=0;explabel=string("");bp_colctr=0;
				null_shift=0.0;minmaxscale=0.0;spacing=0;ymin=0.0;ymax=0.0;
				bp_whiskerlength = 1.5; bp_notched = TRUE; 
				histogram_bins = 0;scat_ctr=0;gridswapdefault=false;
				Ylabel_fmt = string("");
				for(int ctr=0; ctr<64; ++ctr)
					sc[ctr]=MISCPLOT::sc_init[ctr];
      };

      //destructor
      ~miscplot() {
				GDCglobals_reset();
      }    
			
			static void Timeseries(const NEWMAT::Matrix& mat, string filename, 
				string title, float tr = 0.0, int ysize = 150, int width = 4, 
				int prec = 2, bool sci = false)
				{	miscplot plot;
					plot.timeseries(mat, filename, title, tr, ysize, width, prec, sci);
				}

      void timeseries(const NEWMAT::Matrix& mat, string filename, 
				string title, float tr = 0.0, int ysize = 150, int width = 4, 
				int prec = 2, bool sci = false);

      void histogram(const NEWMAT::Matrix& mat, string filename, string title);

      void boxplot(string filename, string title);
      void boxplot(const NEWMAT::Matrix& mat, string filename, string title);
      void boxplot(const NEWMAT::ColumnVector& vec, string filename, 
				string title);

      void gmmfit(const NEWMAT::Matrix& mat, Matrix& mu, Matrix& sig, 
				Matrix& pi, string filename, string title, bool mtype = false, 
        float offset = 0.0, float detailfactor = 0.0); 
    
      inline void ggmfit(const NEWMAT::Matrix& mat, Matrix& mu, 
				Matrix& sig, Matrix& pi, string filename, string title, 
				float offset = 0.0, float detailfactor = 0.0){
					this->gmmfit(mat, mu, sig, pi, filename, title, true, 
					offset, detailfactor);
				}
    
      inline void add_label(string txt){ 
				labels.push_back(txt);}
      inline void remove_labels(int i) {
				for (int j=1;j<=i;j++) labels.pop_back();}
			inline void clear_labels(){
				labels.clear();}
      inline void add_xlabel(string txt){ 
				xlabels.push_back(txt);}
      inline void remove_xlabel(){ 
				xlabels.pop_back();}
			inline void clear_xlabel(){
					xlabels.clear();}
      inline void add_ylabel(string txt){ 
				ylabels.push_back(txt);}
			inline void clear_ylabel(){
				ylabels.clear();}
      void setscatter(Matrix &data, int width=15);
      void deletescatter();
			void GDCglobals_reset();
      void add_bpdata(const NEWMAT::Matrix& mat);
      void add_bpdata(const NEWMAT::ColumnVector& vec);
      inline void set_xysize(int xsize, int ysize){ 
				req_xsize=xsize; req_ysize=ysize;}
      inline void set_nullshift(double val){
				null_shift=val;};
      inline void set_minmaxscale(float val){
				minmaxscale = val;};
      inline void set_yrange(float val1, float val2){
				ymin = val1; ymax=val2;};
      inline void set_spacing(int val){
				spacing = val;};
      inline void set_bpnotches(bool val){
				bp_notched  = val;};
      inline void set_bpwhiskerlength(float val){
				bp_whiskerlength = val;};
      inline void set_histogram_bins(int val){
				histogram_bins = val;};
			inline void grid_swapdefault(){
				gridswapdefault = true;};
			inline void set_Ylabel_fmt(string what){
				Ylabel_fmt = what;}
			inline void col_replace(int which, unsigned long what){
				if(which>=0 && which<64)
					sc[which] = what;
			}

    private:
			unsigned long sc[64];
      vector<string> labels;
      vector<string> xlabels;
      vector<string> ylabels;
      vector<float> bp_min;
      vector<float> bp_max;
      vector<float> bp_median;
      vector<float> bp_medhi;
      vector<float> bp_medlo;
      vector<float> bp_wishi;
      vector<float> bp_wislo;
      vector<float> bp_iqr;
      vector<float> bp_q1;
      vector<float> bp_q3;
      vector<int> bp_outlierindex;
      vector<float> bp_outliervalue;
      int scat_ctr;

      string explabel;
      int req_xsize,req_ysize;
      double null_shift;
      double minmaxscale;
      float ymin,ymax;
			string Ylabel_fmt;
			
      int spacing;
      int bp_colctr;
      int histogram_bins;

      bool bp_notched;
			bool gridswapdefault;
      float bp_whiskerlength, bp_maxall, bp_minall;

      gdImagePtr outim;
      void add_legend(void* ptr, unsigned long cmap[], bool inside=false);
  };
}
#endif
