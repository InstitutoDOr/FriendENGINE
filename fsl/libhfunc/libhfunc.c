#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "libhfunc.h"

int avw_read(char* filename,struct dsr *header)
{
  FILE *fp;
  char fname[1024];
		
  xtrt_filepath(filename,fname);
  strcat(fname,".hdr");
  if((fp=fopen(fname,"rb"))==NULL)
    return(-1);

  if(!fread(header,sizeof(*header),1,fp))
    return(-1);

  fclose(fp);
  return(0);
}
int avw_write(char* filename,struct dsr *header)
{
  FILE *fp;
  char fname[1024];

  xtrt_filepath(filename,fname);
  strcat(fname,".hdr");
  if((fp=fopen(fname,"w+b"))==NULL)
    return(-1);	

  if(fwrite(header,sizeof(*header),1,fp)!=1)
    return(-1);

  fclose(fp);
  return(0);
}

int avw_xdim(struct dsr *header){return(header->dime.dim[1]);}
int avw_ydim(struct dsr *header){return(header->dime.dim[2]);}
int avw_zdim(struct dsr *header){return(header->dime.dim[3]);}
int avw_vdim(struct dsr *header){return(header->dime.dim[4]);}
int avw_dt(struct dsr *header){return(header->dime.datatype);}
float avw_xvox(struct dsr *header){return(header->dime.pixdim[1]);}
float avw_yvox(struct dsr *header){return(header->dime.pixdim[2]);}
float avw_zvox(struct dsr *header){return(header->dime.pixdim[3]);}
float avw_tr(struct dsr *header){return(header->dime.pixdim[4]);}
char avw_orient(struct dsr *header){return(header->hist.orient);}

void avw_set_dim(struct dsr *header,int x,int y,int z,int v)
{
  header->dime.dim[1]=x;
  header->dime.dim[2]=y;
  header->dime.dim[3]=z;
  header->dime.dim[4]=v;
}
void avw_get_dim(struct dsr *header,int *x,int *y,int *z,int *v)
{
  *x=header->dime.dim[1];
  *y=header->dime.dim[2];
  *z=header->dime.dim[3];
  *v=header->dime.dim[4];
}
void avw_set_dt(struct dsr *header,int dt)
{
  switch(dt){	
  case(DT_UNSIGNED_CHAR):
    header->dime.datatype=DT_UNSIGNED_CHAR;
    header->dime.bitpix=8;
    break;
  case(DT_SIGNED_SHORT):
    header->dime.datatype=DT_SIGNED_SHORT;
    header->dime.bitpix=16;
    break;
  case(DT_FLOAT):
    header->dime.datatype=DT_FLOAT;
    header->dime.bitpix=32;
    break;
  case(DT_COMPLEX):
    header->dime.datatype=DT_COMPLEX;
    header->dime.bitpix=64;
    break;
  default:
    header->dime.datatype=dt;
  }
}
void avw_get_dt(struct dsr *header,int *dt)
{
  *dt=header->dime.datatype;
}
void avw_set_maxmin(struct dsr *header,int max,int min)
{
  header->dime.glmax=max;
  header->dime.glmin=min;
}
void avw_get_maxmin(struct dsr *header,int *max,int *min)
{
  *max=header->dime.glmax;
  *min=header->dime.glmin;
}
void avw_set_vox(struct dsr *header,float vx,float vy,float vz)
{
  header->dime.pixdim[1]=vx;
  header->dime.pixdim[2]=vy;
  header->dime.pixdim[3]=vz;
}
void avw_get_vox(struct dsr *header,float *vx,float *vy,float *vz)
{
  *vx=header->dime.pixdim[1];
  *vy=header->dime.pixdim[2];
  *vz=header->dime.pixdim[3];
}
void avw_set_orient(struct dsr *header,char orient)
{
  header->hist.orient=orient;
}
void avw_get_orient(struct dsr *header,char *orient)
{
  *orient=header->hist.orient;
}
void avw_set_tr(struct dsr *header,float tr)
{
  header->dime.pixdim[4]=tr;
}
void avw_get_tr(struct dsr *header,float *tr)
{
  *tr=header->dime.pixdim[4];
}
void avw_set_study(struct dsr *header,char *string)
{
  strncpy(header->hk.db_name,string,17);
  header->hk.db_name[17]='\0';
}
void avw_get_study(struct dsr *header,char *string)
{
  strcpy(string,header->hk.db_name);
}
void avw_set_scannum(struct dsr *header,char *string)
{
  strncpy(header->hist.scannum,string,9);
  header->hist.scannum[9]='\0';
}
void avw_get_scannum(struct dsr *header,char *string)
{
  strcpy(string,header->hist.scannum);
}
void avw_set_date(struct dsr *header,char *string)
{
  strncpy(header->hist.exp_date,string,9);
  header->hist.exp_date[9]='\0';
}
void avw_get_date(struct dsr *header,char *string)
{
  strcpy(string,header->hist.exp_date);
}
void avw_set_patient(struct dsr *header,char *string)
{
  strncpy(header->hist.patient_id,string,9);
  header->hist.patient_id[9]='\0';
}
void avw_get_patient(struct dsr *header,char *string)
{
  strcpy(string,header->hist.patient_id);
}
void avw_get_descrip(struct dsr *header,char *string)
{
  strcpy(string,header->hist.descrip);
}
void avw_set_descrip(struct dsr *header,char *string)
{
  strncpy(header->hist.descrip,string,79);
  header->hist.descrip[79]='\0';
}
int avw_compare_dim(struct dsr *header1,struct dsr *header2)
{
  int i;
  for(i=1;i<=4;i++){
    if(header1->dime.dim[i]!=header2->dime.dim[i])return(1);
  }
  return(0);
}
void avw_swap_xy(struct dsr *header)
{
  int x;
  float vx;
  x=header->dime.dim[1];
  header->dime.dim[1]=header->dime.dim[2];
  header->dime.dim[2]=x;
  vx=header->dime.pixdim[1];
  header->dime.pixdim[1]=header->dime.pixdim[2];
  header->dime.pixdim[2]=vx;
}
void avw_init(struct dsr *header,int x,int y,int z,int v,int t)
{
  short *sbuf;

  header->hk.sizeof_hdr=348;
  strcpy(header->hk.data_type,"");		
  strcpy(header->hk.db_name,"");
  header->hk.extents=16384;
  header->hk.session_error=0;
  header->hk.regular='r';
  header->hk.hkey_un0=' ';

  header->dime.dim[0]=4;
  header->dime.dim[1]=x;
  header->dime.dim[2]=y;
  header->dime.dim[3]=z;
  header->dime.dim[4]=v;
  header->dime.dim[5]=0;
  header->dime.dim[6]=0;
  header->dime.dim[7]=0;
  strcpy(header->dime.vox_units,"mm");
  strcpy(header->dime.cal_units,"");
  header->dime.unused1=0;
  header->dime.dim_un0=0;
  header->dime.pixdim[0]=0.0;
  header->dime.pixdim[1]=1.0;
  header->dime.pixdim[2]=1.0;
  header->dime.pixdim[3]=1.0;
  header->dime.pixdim[4]=0.0;
  header->dime.pixdim[5]=0.0;
  header->dime.pixdim[6]=0.0;
  header->dime.pixdim[7]=0.0;
  header->dime.vox_offset=0.0;
  header->dime.funused1=1.0;
  header->dime.funused2=0.0;
  header->dime.funused3=0.0;
  header->dime.cal_max=0.0;
  header->dime.cal_min=0.0;
  header->dime.compressed=0;
  header->dime.verified=0;
  header->dime.datatype=t;
  switch(t){	
  case(DT_UNSIGNED_CHAR):
    header->dime.datatype=DT_UNSIGNED_CHAR;
    header->dime.bitpix=8;
    header->dime.glmax=255;
    header->dime.glmin=0;
    break;
  case(DT_SIGNED_SHORT):
    header->dime.datatype=DT_SIGNED_SHORT;
    header->dime.bitpix=16;
    header->dime.glmax=32767;
    header->dime.glmin=-32767;
    break;
  case(DT_FLOAT):
    header->dime.datatype=DT_FLOAT;
    header->dime.bitpix=32;
    header->dime.glmax=1;
    header->dime.glmin=-1;
    break;
  case(DT_COMPLEX):
    header->dime.datatype=DT_COMPLEX;
    header->dime.bitpix=64;
    header->dime.glmax=1;
    header->dime.glmin=-1;
    break;
  }

  strcpy(header->hist.descrip,"");
  strcpy(header->hist.aux_file,"");
  header->hist.orient=0;
  sbuf=(short *)calloc(5,sizeof(short));
  memcpy((void *)header->hist.originator,(void *)sbuf,5*sizeof(short));
  free(sbuf);
  strcpy(header->hist.generated,"");
  strcpy(header->hist.scannum,"");
  strcpy(header->hist.patient_id,"");
  strcpy(header->hist.exp_date,"");
  strcpy(header->hist.exp_time,"");
  strcpy(header->hist.hist_un0,"");
  header->hist.views=0;
  header->hist.vols_added=0;
  header->hist.start_field=0;
  header->hist.field_skip=0;
  header->hist.omax=0;
  header->hist.omin=0;
  header->hist.smax=0;
  header->hist.smin=0;  
}
void xtrt_filepath(char *in,char *out)
{
  int len;	
  len=strlen(in);
  if(!strcmp(in+len-4,".img"))len-=4;
  if(!strcmp(in+len-4,".hdr"))len-=4;
  strcpy(out,in);
  out[len]='\0';
}
