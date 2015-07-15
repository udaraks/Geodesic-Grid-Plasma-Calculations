#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "silo.h"
#include "stdlib.h"

#define PI acos(-1)

using namespace std;
/*******************************************************************************
 * GLOBAL CONSTANTS:  ******************************************************************************/
const int iprf = 6;
const int iprm1 = 32;
const int iprm2 = 122;
const int iprm3 = 482;
const int iprm4 = 1922;
const int iprm5 = 7682;
const int iprm6 = 30722;
const int iprn2 = 240;
const int iprn3 = 960;
const int iprn4 = 3840;
const int iprn5 = 15360;
const int iprn6 = 61440;
const int ipro = 5;
const int iprv = 7682; //=iprm5
const int ipru = 15360; //=iprn5
const int iprz = 340; //341
const int ipri = 8;
const int iprc = 4;
const int npri = 5;

const int gridlevel_max = 5;  // = ipro
int lmax, gridlevel, nvertx, nfacex;
float *vlat, *vlong ;
float **xv, **yv, **zv,**zzv, **cocx, **cocy, **cocz, **v0, **s[7],**sb ,****t[7],****tb  ;
int **vofv, **voff, **fnxtf , **vofv_new,**voff_new,**voff_new2,*order , *ans , **fofv , **fnxte,**fnxtf_new, **vofe,**eoff,**eofv;
float** buav,**fw,**gw,**hw, *xbuf,*ybuf,*zbuf;
int nv1,i,j,i1,i2,i3,nz1,nz = iprz;
float f1,f2,f3,g1,g2,g3,h1,h2,h3,v1,v2,v3;
float** bu,**bub, **xclc,**xcld,**buc,**bud,*bulc,*buld,**xub[4],**tbc[4],*sbc,**sc,****tc,*vc,*vx,*vy,*vz;
float **plasma_prim[9], **plasma_prim2[9],**plasma_vector[7],**plasma_vect_out[7],**plasma_vect_out2[7];

const char file1_name[] = "pentd1";
const char file2_name[] = "pentd2";
const char file3_name[] = "pentd3";
const char file4_name[] = "dataf1.050910";
const char plasma_silo[] = "vectors.silo";
const char mesh_name[] = "helio_mesh";
const char var1_name[] = "Divu";
const char var2_name[] = "Velocity";
const char var3_name[] = "CurlB_B2";
const char var4_name[] = "Magnetic_field";
char var1_sub1[] = "Divu";
char var2_sub1[] = "U_x";
char var2_sub2[] = "U_y";
char var2_sub3[] = "U_z";
char var3_sub1[] = "CurlB_x";
char var3_sub2[] = "CurlB_y";
char var3_sub3[] = "CurlB_z";
char var4_sub1[] = "B_x";
char var4_sub2[] = "B_y";
char var4_sub3[] = "B_z";

int WriteSiloScalar(DBfile *plamsafile, float **var1, const char *var_name,   char *var_sub1);
int WriteSiloVector(DBfile *plasmafile, float **var1, float **var2, float **var3, const char *var_name, char *var_sub1, char *var_sub2, char *var_sub3);
void vecop0(int imode);
void vform2(int imode ,float** buav ,float** fw ,float** gw ,float** hw) ;
void venocvb1(int imode,int nz1,float** xclc,float** xcld,float* sbc,float** tbc[4],float**buc,float**bud,float* vx,float* vy,float* vz) ;
void venocv1(int imode,int nz1,float**xclc,float**sc,float**** tc,float**buc,float*vx,float*vy,float*vz);

float **Create2DFloat(int n, int m)
{
   float **array;
   int i;
   array = new float *[n + 1];
   array[0] = new float[n * m + 1];
   array[1] = array[0];
   for(i = 2; i <= n; i++) array[i] = array[i - 1] + m;
   return array;
};

void Delete2DFloat(float **array)
{
   delete[] array[0];
   delete[] array;
};

float ****Create4DFloat(int d1, int d2, int d3, int d4)
{
   float ****arr;

   int i1,i2,i3;

   arr = new float ***[d1 + 1];
   arr[0] = new float **[d1 * d2 + 1];
   arr[0][0] = new float *[d1 * d2 * d3 + 1];
   arr[0][0][0] = new float [d1 * d2 * d3 * d4 + 1];

   arr[1] = arr[0];
   arr[1][1] = arr[0][0];
   arr[1][1][1] = arr[0][0][0];

   for( i3 = 2; i3 <= d3; i3++ ) {
	arr[1][1][i3] = arr[1][1][i3-1] + d4;
   }

   for( i2 = 2; i2 <= d2; i2++ ) {
	arr[1][i2] = arr[1][i2-1] + d3;
	for( i3 = 1; i3 <= d3; i3++ ) {
		arr[1][i2][i3] = arr[1][i2][i3-1] + d4;
	}
   }

   for(i1 = 2; i1 <= d1; i1++) {
		arr[i1] = arr[i1 - 1] + d2;

	   for(i2 = 1; i2 <= d2; i2++) {
			arr[i1][i2] = arr[i1][i2 - 1] + d3;

		   for(i3 = 1; i3 <= d3; i3++) { 
				arr[i1][i2][i3] = arr[i1][i2][i3-1] + d4;

		}
	}
   }		
	
// float ****arr = new float [d1 + 1][d2 + 1][d3 + 1][d4 + 1];
   return arr;
}

void Delete4DFloat(float ****arr)
{
   delete[] arr[0][0][0];
   delete[] arr[0][0];
   delete[] arr[0];
   delete[] arr;
}

int **Create2DInt(int n, int m)
{
   int **array;
   int i;
   array = new int *[n + 1];
   array[0] = new int[n * m + 1];
   array[1] = array[0];
   for(i = 2; i <= n; i++) array[i] = array[i - 1] + m;
   return array;
};

void Delete2DInt(int **array)
{
   delete[] array[0];
   delete[] array;
};

int Pow2(int exponent)
{
   int i;
   double power = 1;
   for(i = 1; i <= exponent; i++) power *= 2;
   return power;
};

bool Match3(int i1, int i2, int i3, int j1, int j2, int j3)
{
   if((i1 != j1) && (i1 != j2) && (i1 != j3)) return false;
   if((i2 != j1) && (i2 != j2) && (i2 != j3)) return false;
   if((i3 != j1) && (i3 != j2) && (i3 != j3)) return false;
   return true;
};

int main(void)
{
   bool found;
   int var, nvertx_max, nfacex_max, vert, vert1,vert2, vert3,vert4, face, nfacex_cur, nedgex, iv, ic, ivmax, l;
   nvertx_max = 30 * Pow2(2 * (gridlevel_max - 1)) + 2; // 30 * 2^(2*4) +2  =  7682  =  iprm5  = iprv (no of vertices)
   nfacex_max = (nvertx_max - 2) * 2; // 15360 = iprn5 = ipru (no of triangle)
//------------------------------------------------------------------------------
   ifstream parmfile;
   parmfile.open(file1_name, ifstream::in | ifstream::binary);

// skip 4 bytes plus nfc
   parmfile.seekg(4);
   parmfile.seekg(sizeof(int), ifstream::cur); // skipping nfc

// grid size
   parmfile.read((char *)&gridlevel, sizeof(int));     // norder
   parmfile.read((char *)&nvertx, sizeof(int));	       // nv
   parmfile.read((char *)&nfacex, sizeof(int));        // nu 
   parmfile.read((char *)&lmax, sizeof(int)); lmax++;  // nz

//------------------------------------------------------------------------------
   vlat = new float[nvertx + 1];
   vlong = new float[nvertx + 1];
   vofv = Create2DInt(6, nvertx);
   voff = Create2DInt(3, nfacex);
   fnxtf = Create2DInt(3, nfacex);
   ifstream vertfile;
   vertfile.open(file2_name, ifstream::in | ifstream::binary);

// skip 4 bytes, read vertex latitudes
   vertfile.seekg(4);
   vertfile.read((char *)&vlat[1], nvertx * sizeof(float));                        //te

// read vertex longitudes
   vertfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
   vertfile.read((char *)&vlong[1], nvertx * sizeof(float));                       //fi

// skip to the start of the vofv record, position at vofv[iv][1],
// or we would overwrite the last record, which is the same as 0th
   vertfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
   for(iv = 1; iv <= 6; iv++) {
      vertfile.read((char *)&vofv[iv][1], nvertx * sizeof(int));                  // icon
      vertfile.seekg((nvertx_max - nvertx) * sizeof(int), ifstream::cur);

   };

// skip ifc (???) for now
   vertfile.seekg(24 * nvertx_max * sizeof(int), ifstream::cur);

// skip ra (???) for now
   vertfile.seekg(24 * nvertx_max * sizeof(float), ifstream::cur);

// read voff
   for(iv = 1; iv <= 3; iv++) {
      vertfile.read((char *)&voff[iv][1], nfacex * sizeof(int));                  //itri(k,3,2) ????
      vertfile.seekg((nfacex_max - nfacex) * sizeof(int), ifstream::cur);
   };
// read fnxtf
   for(ic = 1; ic <= 3; ic++) {
      vertfile.read((char *)&fnxtf[ic][1], nfacex * sizeof(int));                 //itri(k,3,2) ????
      vertfile.seekg((nfacex_max - nfacex) * sizeof(int), ifstream::cur);
   };
   vertfile.close();

//------------------------------------------------------------------------------
   xv = Create2DFloat(lmax, nvertx);
   yv = Create2DFloat(lmax, nvertx);
   zv = Create2DFloat(lmax, nvertx);
   zzv = Create2DFloat(lmax, nvertx);
   v0 = Create2DFloat(lmax, nvertx);

   for(var = 1; var <= 6; var++) 
   {
	s[var] = Create2DFloat(lmax, nvertx);  

   }

   sb = Create2DFloat(lmax, nvertx);

   for(var = 1; var <= 6; var++) 
   {
	t[var] = Create4DFloat(lmax, 3, 3, nvertx);  
   }

   tb = Create4DFloat(lmax, 3, 3, nvertx);  

   ifstream coorfile;
   coorfile.open(file3_name, ifstream::in | ifstream::binary);

   for(l = 1; l <= lmax; l++) {
      coorfile.seekg(4, ifstream::cur);

      coorfile.read((char *)&xv[l][1], nvertx * sizeof(float));   //xcl1(1) - x
      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
      coorfile.read((char *)&yv[l][1], nvertx * sizeof(float));   //xcl1(2) - y
      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
      coorfile.read((char *)&zv[l][1], nvertx * sizeof(float));   //xcl1(3) -z
      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
      coorfile.read((char *)&zzv[l][1], nvertx * sizeof(float));   //xcl1(4) -???
      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
      coorfile.read((char *)&v0[l][1], nvertx * sizeof(float));   // volume 
      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);

      for(var = 1; var <= 6; var++) {   // side areas
   	coorfile.read((char *)&s[var][l][1], nvertx * sizeof(float));   
      	coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
      }

      coorfile.read((char *)&sb[l][1], nvertx * sizeof(float));   // bottom area
      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);

      // side normal/tangential vectors
      for(var = 1; var <= 6; var++) {
   			
     	for(ic=1 ; ic<=3 ; ic++)  {

   	  	for(iv=1 ; iv<=3 ; iv++)  {
			      coorfile.read((char *)&t[var][l][ic][iv][1], nvertx * sizeof(float));   
			      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
		}	
	 }
    }

	// bottom normal/tangential vectors
	for(ic=1 ; ic<=3 ; ic++)
	{
		for(iv=1 ; iv<=3 ; iv++)
		{   	
		      coorfile.read((char *)&tb[l][ic][iv][1], nvertx * sizeof(float));   
		      coorfile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);	
		}
	}
      coorfile.seekg(4, ifstream::cur);

   };
   coorfile.close();


//------------------------------------------------------------------------------
   for(var = 1; var <= 8; var++) 
   {
	plasma_prim[var] = Create2DFloat(lmax, nvertx);
        plasma_prim2[var] = Create2DFloat(lmax, nvertx);

   }

   ifstream datafile;
   datafile.open(file4_name, ifstream::in | ifstream::binary);
   datafile.seekg(103 * sizeof(float) + sizeof(int) + 8, ifstream::cur);
   for(l = 1; l <= lmax; l++) {
      datafile.seekg(4, ifstream::cur);
      for(var = 1; var <= 8; var++) {
         datafile.read((char *)&plasma_prim[var][l][1], nvertx * sizeof(float));
         datafile.seekg((nvertx_max - nvertx) * sizeof(float), ifstream::cur);
      };
      datafile.seekg(4, ifstream::cur);
   };
   datafile.close();

lmax = lmax-1;

// In data file it gives velocity * density

   for(l = 1; l <= lmax; l++) {
      for(vert = 1; vert <= (nvertx); vert++) {

  		plasma_prim[2][l][vert] = plasma_prim[2][l][vert] / plasma_prim[1][l][vert];
  		plasma_prim[3][l][vert] = plasma_prim[3][l][vert] / plasma_prim[1][l][vert];
  		plasma_prim[4][l][vert] = plasma_prim[4][l][vert] / plasma_prim[1][l][vert];

      }
   }

/*
    ofstream file;
    file.open ("pt.dat");
    file.precision(6);
    file.setf(ios::fixed | ios::showpoint);
    cout.precision(10);
    cout.setf(ios::fixed | ios::showpoint);
*/

   for(var = 1; var <= 6; var++) 
   {
	plasma_vector[var] = Create2DFloat(lmax, nvertx);
        plasma_vect_out[var] = Create2DFloat(lmax, nvertx);
        plasma_vect_out2[var] = Create2DFloat(lmax, nvertx);
   }


//**********
//calculate Div u

   for(l = 1; l <= lmax; l++) {
      for(vert = 1; vert <= (nvertx); vert++) {

  		plasma_vector[1][l][vert] = plasma_prim[2][l][vert] ;
  		plasma_vector[2][l][vert] = plasma_prim[3][l][vert] ;
  		plasma_vector[3][l][vert] = plasma_prim[4][l][vert] ;

      }
   }

   vecop0(1);  // calculation

   for(l = 1; l <= lmax; l++) {
      for(vert = 1; vert <= (nvertx); vert++) {

  		plasma_vect_out[1][l][vert] = plasma_vector[4][l][vert] ;
  		plasma_vect_out[2][l][vert] = plasma_vector[5][l][vert] ;
  		plasma_vect_out[3][l][vert] = plasma_vector[6][l][vert] ;

      }
   }

//**********
//calculate Curl B/B^2

   for(l = 1; l <= lmax; l++) {
      for(vert = 1; vert <= (nvertx); vert++) {
  		plasma_vector[1][l][vert] = plasma_prim[5][l][vert] ;
  		plasma_vector[2][l][vert] = plasma_prim[6][l][vert] ;
  		plasma_vector[3][l][vert] = plasma_prim[7][l][vert] ;
      }
   }
   vecop0(2); // calculation of curl


   for(l = 1; l <= lmax; l++) {
      for(vert = 1; vert <= (nvertx); vert++) {
 		plasma_vect_out[4][l][vert] = plasma_vector[4][l][vert] /  
(plasma_prim[5][l][vert]*plasma_prim[5][l][vert] + plasma_prim[6][l][vert]*plasma_prim[6][l][vert] + plasma_prim[7][l][vert]*plasma_prim[7][l][vert])  ;

  		plasma_vect_out[5][l][vert] = plasma_vector[5][l][vert] /  
(plasma_prim[5][l][vert]*plasma_prim[5][l][vert] + plasma_prim[6][l][vert]*plasma_prim[6][l][vert] + plasma_prim[7][l][vert]*plasma_prim[7][l][vert])  ;

  		plasma_vect_out[6][l][vert] = plasma_vector[6][l][vert] / 
 (plasma_prim[5][l][vert]*plasma_prim[5][l][vert] + plasma_prim[6][l][vert]*plasma_prim[6][l][vert] + plasma_prim[7][l][vert]*plasma_prim[7][l][vert])  ;

      }
   }

//**********
   float x1,x2,x3,y1,y2,y3,z1,z2,z3,a1,b1,c1,c,d1,e1,coc1,x4,y4,z4,x5,y5,z5,x6,y6,z6,x12,y12,z12,x13,y13,z13,r5;
   cocx = Create2DFloat(lmax, nvertx*2);
   cocy = Create2DFloat(lmax, nvertx*2);
   cocz = Create2DFloat(lmax, nvertx*2);
   voff_new = Create2DInt(6, nvertx); // vertices of hex/pent faces
//lmax = 10 ;

// Finding centers of triangles

int ind2 , ind3, i,count,temp ,j;

	for(l = 1; l <= lmax-1; l++) {

		ind2 = 0,ind3 =1, i,count=0,temp=0 ,j=1;

		for(vert = 1; vert <= nvertx; vert++) {

			//points of a triangle   // points between top and bottom layers
			x1 = (xv[l][vert] + xv[l+1][vert])*0.5;
			y1 = (yv[l][vert] + yv[l+1][vert])*0.5;
			z1 = (zv[l][vert] + zv[l+1][vert])*0.5;
			
			count = 0.0;

			 
			 // find points around the mid point

			ivmax = (vofv[6][vert] == 0 ? 5 : 6); 

    

			for(iv = 1; iv <= ivmax; iv++) {

				ind2++;// vertex index

				voff_new[iv][ind3] = ind2; // initial assign


	
			// old adjacent points
       			 vert2 = vofv[iv][vert];
        		 vert3 = (iv == ivmax ? vofv[1][vert] : vofv[iv + 1][vert]);

			//vert2 = vofv[1][vert]; // next vertex of vert
			x2 = (xv[l][vert2] + xv[l+1][vert2])*0.5;
			y2 = (yv[l][vert2] + yv[l+1][vert2])*0.5;
			z2 = (zv[l][vert2] + zv[l+1][vert2])*0.5;

			//vert3 = vofv[2][vert]; // 2nd vertex connected 2 vert
			x3 = (xv[l][vert3] + xv[l+1][vert3])*0.5;
			y3 = (yv[l][vert3] + yv[l+1][vert3])*0.5;
			z3 = (zv[l][vert3] + zv[l+1][vert3])*0.5;

		//   side lengths 2-3, 3-1, 1-2 of the triangle

		      a1 = pow( (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3), 0.5);
		      b1 = pow( (x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+(z3-z1)*(z3-z1), 0.5);
		      c1 = pow( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2), 0.5);

		//    distance from side 3-1 to the center of outer circle 

		      c = acos( -(c1*c1-a1*a1-b1*b1) / (2.0*a1*b1) );
		      d1 = b1/2.0 - a1/2.0 * cos(c);
		      e1 = a1/2.0 * sin(c);
		      coc1 = e1 - d1 * tan(PI/2.0-c);

		//   cross product - to get the vector normal to 1-3 surface
		      x12 = (x2-x1);
		      y12 = (y2-y1);
		      z12 = (z2-z1);
		      x13 = (x3-x1);
		      y13 = (y3-y1);
		      z13 = (z3-z1);
		      x4 = y13*z12-z13*y12;
		      y4 = z13*x12-x13*z12;
		      z4 = x13*y12-y13*x12;
		      x5 = y4*z13-z4*y13;
		      y5 = z4*x13-x4*z13;
		      z5 = x4*y13-y4*x13;
		      r5 = pow( (x5*x5+y5*y5+z5*z5), 0.5);

		//    center of side 3-1
		      x6=(x3+x1)*0.5;
		      y6=(y3+y1)*0.5;
		      z6=(z3+z1)*0.5;


         
		//    Centers of triangles 
		      cocx[l][ind2] = 50.0 * (x6 + coc1* x5/r5);
		      cocy[l][ind2] = 50.0 * (y6 + coc1* y5/r5);
		      cocz[l][ind2] = 50.0 * (z6 + coc1* z5/r5);

			for(i = 1; i < ind2; i++){	
				if( (fabs(cocz[l][ind2] - cocz[l][i]) <=0.01) && (fabs(cocy[l][ind2] - cocy[l][i]) <=0.01)  && (fabs(cocx[l][ind2] -cocx[l][i]) <=0.01) ) { 
					voff_new[iv][ind3] = i; 
					count++;
					temp = -1;
					break;
				}  		
			}

			ind2 = ind2 + temp;  // stop increment if its the same vertex
			
			temp = 0;	

			j++;	
		} //iv
		 ind3++;
		}

	}
	cout<<"Total Hex/Pent faces: "<<nvertx<<"	Total vertices: "<<ind2<<endl;

	nfacex = nvertx;  // New number of faces - center of face is nvertx
	nvertx = ind2;  // new no. of vertices (vertices of each pent/hex)
	nedgex = 3 * nvertx / 2;  // no of edges

	voff_new2 = Create2DInt(6, nfacex); // vertices of hex/pent faces

	int k = 12;


// put the 12 pent faces at the beginning of voff_new2 array

// get the hex elements of voff
 for(vert = 1; vert <= nfacex; vert++) {

	if(voff_new[6][vert] != 0)
	{
			k++;
		for(iv=1;iv<=6; iv++)
		{
			voff_new2[iv][k] = voff_new[iv][vert];
		}
	}
}

k=0;
// get the pentagons to 1st 12 elements of voff

 for(vert = 1; vert <= nfacex; vert++) {

	if(voff_new[6][vert] == 0)
	{
		//cout<<vert<<endl;
			k++;
		for(iv=1;iv<=6; iv++)
		{
			voff_new2[iv][k] = voff_new[iv][vert];
		}
	}
}

// put the 12 pent faces data values at the beginning of data arrays

 for(l = 1; l <= lmax; l++) {
	k=12;
      for(vert = 1; vert <= nfacex; vert++) {
	if(voff_new[6][vert] != 0)
	  {
			k++;
     	 	for(var = 1; var <= 8; var++) {
         		plasma_prim2[var][l][k] = plasma_prim[var][l][vert]; 
	        }

     	 	for(var = 1; var <= 6; var++) {
			plasma_vect_out2[var][l][k] = plasma_vect_out[var][l][vert];	 
	       }
	  }

	};

   };

 for(l = 1; l <= lmax; l++) {
	k=0;
      for(vert = 1; vert <= nfacex; vert++) {
	if(voff_new[6][vert] == 0)
	  {
			k++;
     	 	for(var = 1; var <= 8; var++) {
         		plasma_prim2[var][l][k] = plasma_prim[var][l][vert];         
	      }
     	 	for(var = 1; var <= 6; var++) {
			plasma_vect_out2[var][l][k] = plasma_vect_out[var][l][vert];	
	       }
	  }

	};

   };


int edge, edge1, edge2, face1, face2,face3, face4, iv1, iv2, tempi1;
////////////////////////////////

// connectivity arrays

   fnxtf_new = Create2DInt(6, nfacex);
   eoff = Create2DInt(6, nfacex);
   fnxte = Create2DInt(2, nedgex);
   vofe = Create2DInt(2, nedgex);
   fofv = Create2DInt(3, nvertx);
   eofv = Create2DInt(3, nvertx);

// calculate fofv

   memset(fofv[0], 0, (nfacex) * sizeof(int));

   for(face = 1; face <= nfacex; face++) {
      for(iv = 1; iv <= (face > 12 ? 6 : 5); iv++) {
         vert = voff_new2[iv][face];
         if(!fofv[1][vert]) fofv[1][vert] = face;
         else if(!fofv[2][vert]) fofv[2][vert] = face;
         else fofv[3][vert] = face;
      };
   };

// calculate fnxtf_new, fnxte, vofe
   edge = 1;
   for(face = 1; face <= nfacex; face++) {
      for(iv1 = 1; iv1 <= (face > 12 ? 6 : 5); iv1++) {
         iv2 = iv1 + 1;
         if((iv1 == 5 && face <= 12) || (iv1 == 6 && face > 12)) iv2 = 1;
         vert1 = voff_new2[iv1][face];
         vert2 = voff_new2[iv2][face];
         face1 = fofv[1][vert1]; // two neighbors of this face sharing vert1
         if(face1 == face) {
            face1 = fofv[2][vert1];
            face2 = fofv[3][vert1];
         }
         else {
            face2 = fofv[2][vert1];
            if(face2 == face) face2 = fofv[3][vert1];
         };
         face3 = fofv[1][vert2]; // two neighbors of this face sharing vert2
         if(face3 == face) {
            face3 = fofv[2][vert2];
            face4 = fofv[3][vert2];
         }
         else {
            face4 = fofv[2][vert2];
            if(face4 == face) face4 = fofv[3][vert2];
         };
         if(face1 == face3 || face1 == face4) fnxtf_new[iv1][face] = face1;
         else fnxtf_new[iv1][face] = face2;
         if(face < fnxtf_new[iv1][face]) { // insert an edge
            fnxte[1][edge] = face;
            fnxte[2][edge] = fnxtf_new[iv1][face]; // fnxte[1] < fnxte[2]
            vofe[1][edge] = vert1;
            vofe[2][edge] = vert2;
            edge++;
         };
      };
   };

// calculate eoff
   for(edge = 1; edge <= nedgex; edge++) {
      face1 = fnxte[1][edge];
      face2 = fnxte[2][edge];
      ic = 1;
      while(fnxtf_new[ic][face1] != face2) ic++;
      eoff[ic][face1] = edge;
      ic = 1;
      while(fnxtf_new[ic][face2] != face1) ic++;
      eoff[ic][face2] = edge;

   };

// calculate eofv
   memset(eofv[0], 0, (3 * nvertx + 1) * sizeof(int));
   for(edge = 1; edge <= nedgex; edge++) {
      for(iv = 1; iv <= 2; iv++) {
         vert = vofe[iv][edge];
         if(!eofv[1][vert]) eofv[1][vert] = edge;
         else if(!eofv[2][vert]) eofv[2][vert] = edge;
         else eofv[3][vert] = edge;
      };
   };

//////////////////////////////////

//------------------------------------------------------------------------------
   int nnodes, nzones, lznodelist, err, ind;
   int *znodelist;
   int zshapecnt[2], zshapesize[2], zshapetype[2];
   float *coords[3];
   DBfile *plasmafile;
   DBfacelist *fl;
   plasmafile = DBCreate(plasma_silo, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);

// vertices
   nnodes = (nvertx) * (lmax - 3);  
   coords[0] = new float[nnodes]; // x
   coords[1] = new float[nnodes]; // y
   coords[2] = new float[nnodes]; // z
   ind = 0;
   for(l = 3; l <= lmax-1; l++) {
      for(vert = 1; vert <= (nvertx); vert++) {
         coords[0][ind] =  (cocx[l][vert]);
         coords[1][ind] =  (cocy[l][vert]);
         coords[2][ind] =  (cocz[l][vert]);
         ind++;
      };
   };

// volume cells
   nzones = 2 * nfacex * (lmax - 4);
   lznodelist = (12 * 6 + 12 * 8 + 2 * (nfacex - 12) * 8) * (lmax - 4);
   znodelist = new int[lznodelist];
   zshapecnt[0] = 12 * (lmax - 4);
   zshapesize[0] = 6;
   zshapetype[0] = DB_ZONETYPE_PRISM;
   zshapecnt[1] = nzones - zshapecnt[0];
   zshapesize[1] = 8;
   zshapetype[1] = DB_ZONETYPE_HEX;
   ind = 0;
// prism elements of pentagonal cells
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 1; face <= 12; face++) {
         edge = eoff[1][face];
         vert1 = vofe[1][edge]; // not shared with the next edge
         vert2 = vofe[2][edge];
         if(vert1 == vofe[1][eoff[2][face]] || vert1 == vofe[2][eoff[2][face]]) {
            vert1 = vofe[2][edge];
            vert2 = vofe[1][edge];
         };
         edge = eoff[2][face];
         vert3 = vofe[1][edge];
         if(vert3 == vert2) vert3 = vofe[2][edge];
         znodelist[ind++] = vert1 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert1 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 2) * nvertx - 1;
      };
   };
// hex elements of pentagonal cells
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 1; face <= 12; face++) {
         edge = eoff[3][face];
         vert1 = vofe[1][edge]; // not shared with the next edge
         vert2 = vofe[2][edge];
         if(vert1 == vofe[1][eoff[4][face]] || vert1 == vofe[2][eoff[4][face]]) {
            vert1 = vofe[2][edge];
            vert2 = vofe[1][edge];
         };
         edge = eoff[4][face];
         vert3 = vofe[1][edge];
         if(vert3 == vert2) vert3 = vofe[2][edge];
         edge = eoff[5][face];
         vert4 = vofe[1][edge];
         if(vert4 == vert3) vert4 = vofe[2][edge];
         znodelist[ind++] = vert1 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert4 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert1 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert4 + (l - 2) * nvertx - 1;
      };
   };
// hex elements of hexagonal cells
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 13; face <= nfacex; face++) {
         edge = eoff[1][face];
         vert1 = vofe[1][edge]; // not shared with the next edge
         vert2 = vofe[2][edge];
         if(vert1 == vofe[1][eoff[2][face]] || vert1 == vofe[2][eoff[2][face]]) {
            vert1 = vofe[2][edge];
            vert2 = vofe[1][edge];
         };
         edge = eoff[2][face];
         vert3 = vofe[1][edge];
         if(vert3 == vert2) vert3 = vofe[2][edge];
         edge = eoff[3][face];
         vert4 = vofe[1][edge];
         if(vert4 == vert3) vert4 = vofe[2][edge];
         znodelist[ind++] = vert1 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert4 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert1 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert4 + (l - 2) * nvertx - 1;
         vert1 = vert4;
         edge = eoff[4][face];
         vert2 = vofe[1][edge];
         if(vert2 == vert1) vert2 = vofe[2][edge];
         edge = eoff[5][face];
         vert3 = vofe[1][edge];
         if(vert3 == vert2) vert3 = vofe[2][edge];
         edge = eoff[6][face];
         vert4 = vofe[1][edge];
         if(vert4 == vert3) vert4 = vofe[2][edge];
         znodelist[ind++] = vert1 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert4 + (l - 3) * nvertx - 1;
         znodelist[ind++] = vert1 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert2 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert3 + (l - 2) * nvertx - 1;
         znodelist[ind++] = vert4 + (l - 2) * nvertx - 1;
      };
   };
// external faces
   fl = DBCalcExternalFacelist2(znodelist, nnodes, 0, 0, 0, zshapetype,
      zshapesize, zshapecnt, 2, NULL, 0);
   err = DBPutFacelist(plasmafile, "fl1", fl->nfaces, 3, fl->nodelist,
      fl->lnodelist, 0, fl->zoneno, fl->shapesize, fl->shapecnt, fl->nshapes,
      NULL, NULL, 0);
   err = DBPutZonelist2(plasmafile, "zl1", nzones, 3, znodelist, lznodelist,
      0, 0, 0, zshapetype, zshapesize, zshapecnt, fl->nshapes, NULL);
   err = DBPutUcdmesh(plasmafile, mesh_name, 3, NULL, coords, nnodes, nzones,
     "zl1", "fl1", DB_FLOAT, NULL);

   DBFreeFacelist(fl);
   delete[] znodelist;
   delete[] coords[0];
   delete[] coords[1];
   delete[] coords[2];
//----------------------------------------------------------------------
// write data to silo

// Divu
WriteSiloScalar(plasmafile, plasma_vect_out2[1], var1_name,   var1_sub1);

//velocity
WriteSiloVector(plasmafile, plasma_prim2[2], plasma_prim2[3], plasma_prim2[4],  var2_name, var2_sub1, var2_sub2, var2_sub3);

//Curl B /B^2
WriteSiloVector(plasmafile, plasma_vect_out2[4], plasma_vect_out2[5], plasma_vect_out2[6],  var3_name, var3_sub1, var3_sub2, var3_sub3);

//magnetic field
WriteSiloVector(plasmafile, plasma_prim2[5], plasma_prim2[6], plasma_prim2[7],  var4_name, var4_sub1, var4_sub2, var4_sub3);

//-----------------------------------------------------------------------
   DBClose(plasmafile);
   delete[] vlat; delete[] vlong;
   Delete2DInt(vofv); Delete2DInt(voff); Delete2DInt(fnxtf); Delete2DInt(voff_new);
   Delete2DFloat(xv); Delete2DFloat(yv); Delete2DFloat(zv); Delete2DFloat(zzv); 
   Delete2DFloat(cocx); Delete2DFloat(cocy); Delete2DFloat(cocz);Delete2DFloat(v0);Delete2DFloat(sb);Delete4DFloat(tb);
};

// write one scalar variable to a SILO file
int WriteSiloScalar(DBfile *plamsafile, float **var1, const char *var_name,   char *var_sub1)
{
   int err, l, ind, face;
   int nzones;
   float *vars[1];
   char *varnames[1];
   nzones = 2 * nfacex * (lmax - 4);
   vars[0] = new float[nzones];
   varnames[0] = var_sub1;
   ind = 0;
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 1; face <= 12; face++) vars[0][ind++] = var1[l][face];
   };
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 1; face <= 12; face++) vars[0][ind++] = var1[l][face];
   };
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 13; face <= nfacex; face++) {
         vars[0][ind++] = var1[l][face];
         vars[0][ind++] = var1[l][face];
      };
   };
   err = DBPutUcdvar(plamsafile, var_name, mesh_name, 1, varnames, vars,
      nzones, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
   delete[] vars[0];
   return 0;
};

// write one 3-component vector variable to a SILO file
int WriteSiloVector(DBfile *plasmafile, float **var1, float **var2, float **var3,
   const char *var_name, char *var_sub1, char *var_sub2, char *var_sub3)
{
   int err, l, ind, face;
   int nzones;
   float *vars[3];
   char *varnames[3];
   nzones = 2 * nfacex * (lmax - 4);
   vars[0] = new float[nzones];
   vars[1] = new float[nzones];
   vars[2] = new float[nzones];
   varnames[0] = var_sub1;
   varnames[1] = var_sub2;
   varnames[2] = var_sub3;
   ind = 0;
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 1; face <= 12; face++) {
         vars[0][ind] = var1[l][face];
         vars[1][ind] = var2[l][face];
         vars[2][ind] = var3[l][face];
         ind++;
      };
   };
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 1; face <= 12; face++) {
         vars[0][ind] = var1[l][face];
         vars[1][ind] = var2[l][face];
         vars[2][ind] = var3[l][face];
         ind++;
      };
   };
   for(l = 3; l <= lmax - 2; l++) {
      for(face = 13; face <= nfacex; face++) {
         vars[0][ind] = var1[l][face];
         vars[1][ind] = var2[l][face];
         vars[2][ind] = var3[l][face];
         ind++;
         vars[0][ind] = var1[l][face];
         vars[1][ind] = var2[l][face];
         vars[2][ind] = var3[l][face];
         ind++;
      };
   };
   err = DBPutUcdvar(plasmafile, var_name, mesh_name, 3, varnames, vars,
      nzones, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
   delete[] vars[0]; delete[] vars[1]; delete[] vars[2];
   return 0;
};

//--------------------------------------------------------------------------------------------------
// functions for vector operations
//imode = 1 => divergence
//imode = 2  => curl
//imode = 3 => gradient

void vecop0(int imode)
{
//   input  bu
//     output cu (xub)
	xclc = Create2DFloat(3, iprv);
	xcld = Create2DFloat(3, iprv);
	buc = Create2DFloat(3, iprv);
	bud = Create2DFloat(3, iprv);
	buav = Create2DFloat(3, iprv);	

	for(int var = 1; var <= 3; var++) 
  	{
		xub[var] = Create2DFloat(iprz+1, iprv);
       	 	tbc[var] = Create2DFloat(3, iprv);
  	 }

	tc = Create4DFloat(6, 3 , 3, iprv);
	sc = Create2DFloat(6, iprv);
	sbc = new float[iprv+1];
	vc = new float[iprv+1];
	vx = new float[iprv+1];
	vy = new float[iprv+1];
	vz = new float[iprv+1];
	xbuf = new float[iprv+1];
	ybuf = new float[iprv+1];
	zbuf = new float[iprv+1];
	fw = Create2DFloat(3, iprv);
	gw = Create2DFloat(3, iprv);
	hw = Create2DFloat(3, iprv);

//      load to level 1..................................function.......
//      buc/side, bu1/bottom

for(nz1 =1;nz1<=iprz ; nz1++) {
      for(i1=1 ;i1<=6;i1++) {
      		for(i2=1 ;i2<=3;i2++) {
			for(nv1 =1 ; nv1 <= iprv ; nv1++ )
			{
				 tc[i1][1][i2][nv1]=t[i1][nz1][1][i2][nv1];
				 tc[i1][2][i2][nv1]=t[i1][nz1][2][i2][nv1];
				 tc[i1][3][i2][nv1]=t[i1][nz1][3][i2][nv1];
		   	}
		}

		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
        		 sc[i1][nv1]=s[i1][nz1][nv1];
	   	}
	} // i1

//     (1) vc/volume, sbc/bottom surface, tbc/bottom normal
          
	for(nv1 =1 ; nv1 <= iprv ; nv1++ )
	{
		 vc[nv1]= v0[nz1][nv1];
		 sbc[nv1]=sb[nz1][nv1];
   	}

       for(i2=1 ;i2<=3;i2++) {
	       for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{	
			 tbc[1][i2][nv1]= tb[nz1][1][i2][nv1];
			 tbc[2][i2][nv1]= tb[nz1][2][i2][nv1];
			 tbc[3][i2][nv1]= tb[nz1][3][i2][nv1];
	   	}
   	}

		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		  xclc[1][nv1]=xv[nz1][nv1];
		  xclc[2][nv1]=yv[nz1][nv1];
		  xclc[3][nv1]=zv[nz1][nv1];

		  if(nz1!=1) xcld[1][nv1]=xv[nz1-1][nv1];
		  if(nz1!=1) xcld[2][nv1]=yv[nz1-1][nv1];
		  if(nz1!=1) xcld[3][nv1]=zv[nz1-1][nv1];

		  if(nz1==1) xcld[1][nv1]=xv[nz1][nv1] - 0.1;
		  if(nz1==1) xcld[2][nv1]=yv[nz1][nv1] - 0.1;
		  if(nz1==1) xcld[3][nv1]=zv[nz1][nv1] - 0.1;

		  buc[1][nv1] = plasma_vector[1][nz1][nv1];
		  buc[2][nv1] = plasma_vector[2][nz1][nv1];
		  buc[3][nv1] = plasma_vector[3][nz1][nv1];

		//  bulc[nv1]=bub[nz1][nv1];

		  if(nz1!=1) bud[1][nv1] = plasma_vector[1][nz1-1][nv1];
		  if(nz1!=1) bud[2][nv1] = plasma_vector[2][nz1-1][nv1];
		  if(nz1!=1) bud[3][nv1] = plasma_vector[3][nz1-1][nv1];

		  if(nz1==1) bud[1][nv1] = plasma_vector[1][nz1][nv1];
		  if(nz1==1) bud[2][nv1] = plasma_vector[2][nz1][nv1];
		  if(nz1==1) bud[3][nv1] = plasma_vector[3][nz1][nv1];
	   	}
	
	//     integration.....................................................
	//     side surface
	      venocv1(imode,nz1,xclc,sc,tc,buc,vx,vy,vz);

	//     store for nz1
		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		      plasma_vector[4][nz1][nv1] = vx[nv1];
		      plasma_vector[5][nz1][nv1] = vy[nv1];
		      plasma_vector[6][nz1][nv1] = vz[nv1];
		}

	//     bottom surface
	      venocvb1(imode,nz1,xclc,xcld,sbc,tbc,buc,bud,vx,vy,vz);

	//     store for nz1
		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		      xub[1][nz1][nv1]=vx[nv1];
		      xub[2][nz1][nv1]=vy[nv1];
		      xub[3][nz1][nv1]=vz[nv1];
		}
	}

	//    z loop end

	//  outer surface........................................nz+1
		nz = iprz;
	        nz1= nz+1;

	      for(nv1 =1 ; nv1 <= iprv ; nv1++ )	sbc[nv1]=sb[nz1][nv1];

	      for(i2=1;i2<=3;i2++)
		{
		       for(nv1 =1 ; nv1 <= iprv ; nv1++ )
			{
				 tbc[1][i2][nv1]=tb[nz1][1][i2][nv1];
				 tbc[2][i2][nv1]=tb[nz1][2][i2][nv1];
				 tbc[3][i2][nv1]=tb[nz1][3][i2][nv1];
		  	}
		}

	//     bu1/bottom
	       for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		  xclc[1][nv1]=xv[nz][nv1] + 0.1;
		  xclc[2][nv1]=yv[nz][nv1] + 0.1;
		  xclc[3][nv1]=zv[nz][nv1] + 0.1;

		  xcld[1][nv1]=xv[nz][nv1];
		  xcld[2][nv1]=yv[nz][nv1];
		  xcld[3][nv1]=zv[nz][nv1];

		  buc[1][nv1] = plasma_vector[1][nz][nv1];
		  buc[2][nv1] = plasma_vector[2][nz][nv1];
		  buc[3][nv1] = plasma_vector[3][nz][nv1];

		  bud[1][nv1] = plasma_vector[1][nz][nv1];
		  bud[2][nv1] = plasma_vector[2][nz][nv1];
		  bud[3][nv1] = plasma_vector[3][nz][nv1];
	      }

	      venocvb1(imode,nz1,xclc,xcld,sbc,tbc,buc,bud,vx,vy,vz);

	//     store for nz1
	       for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		      xub[1][nz1][nv1]=vx[nv1];
		      xub[2][nz1][nv1]=vy[nv1];
		      xub[3][nz1][nv1]=vz[nv1];
	       }

//	     data transfer for xub
//	     final summation...................................side+top+bottom

	for(nz1 =1;nz1<=iprz ; nz1++) {
	//	     add bottom and top
		     for(nv1 =1 ; nv1 <= iprv ; nv1++ )
			{
			      plasma_vector[4][nz1][nv1] = plasma_vector[4][nz1][nv1]+xub[1][nz1][nv1]-xub[1][nz1+1][nv1];
			      plasma_vector[5][nz1][nv1] = plasma_vector[5][nz1][nv1]+xub[2][nz1][nv1]-xub[2][nz1+1][nv1];
			      plasma_vector[6][nz1][nv1] = plasma_vector[6][nz1][nv1]+xub[3][nz1][nv1]-xub[3][nz1+1][nv1];
		       }

	//	    divide by volume
		     for(nv1 =1 ; nv1 <= iprv ; nv1++ )
			{
			      plasma_vector[4][nz1][nv1] = plasma_vector[4][nz1][nv1]/v0[nz1][nv1];
			      plasma_vector[5][nz1][nv1] = plasma_vector[5][nz1][nv1]/v0[nz1][nv1];
			      plasma_vector[6][nz1][nv1] = plasma_vector[6][nz1][nv1]/v0[nz1][nv1];
		        }
		//    z loop end
	}
}

void venocv1(int imode,int nz1,float**xclc,float**sc,float****tc,float**buc,float*vx,float*vy,float*vz)
{
//    input     buc
//     output    vx,vy,vz
//     function in flux formula at each grid......................
      vform2(imode,buc,fw,gw,hw);

//     clear vx,vy,vz
	for(nv1 =1 ; nv1 <= iprv ; nv1++ )
	{
	      vx[nv1]=0.0;
	      vy[nv1]=0.0;
	      vz[nv1]=0.0;
	}

//    sum nfc side surfaces......................................i
	 for(i=1;i<=6;i++) 
	{
		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{    
		      j = vofv[i][nv1];
		        if(j != 0) 
			{
		//     function value on surface (f, g, h)
		//     f((a+b)/2)=(f(a)+f(g))/2

			      f1=(fw[1][nv1]+fw[1][j])/2.0;
			      f2=(fw[2][nv1]+fw[2][j])/2.0;
			      f3=(fw[3][nv1]+fw[3][j])/2.0;
			      g1=(gw[1][nv1]+gw[1][j])/2.0;
			      g2=(gw[2][nv1]+gw[2][j])/2.0;
			      g3=(gw[3][nv1]+gw[3][j])/2.0;
			      h1=(hw[1][nv1]+hw[1][j])/2.0;
			      h2=(hw[2][nv1]+hw[2][j])/2.0;
			      h3=(hw[3][nv1]+hw[3][j])/2.0;

			//     normal vector inner product
			      v1=f1*tc[i][1][1][nv1]+g1*tc[i][2][1][nv1]+h1*tc[i][3][1][nv1];
			      v2=f2*tc[i][1][1][nv1]+g2*tc[i][2][1][nv1]+h2*tc[i][3][1][nv1];
			      v3=f3*tc[i][1][1][nv1]+g3*tc[i][2][1][nv1]+h3*tc[i][3][1][nv1];

			      xbuf[nv1]=v1*sc[i][nv1];
			      ybuf[nv1]=v2*sc[i][nv1];
			      zbuf[nv1]=v3*sc[i][nv1];
			}
		}

	//    sum for present cell i
	      for(nv1 =1 ; nv1 <= iprv ; nv1++ ) {
		      j=vofv[i][nv1];

///*******************
		      if(j != 0) 
			{
			      vx[nv1]=vx[nv1]+xbuf[nv1];
			      vy[nv1]=vy[nv1]+ybuf[nv1];
			      vz[nv1]=vz[nv1]+zbuf[nv1];
			}
		}
	}

}

void venocvb1(int imode,int nz1,float** xclc,float** xcld,float* sbc,float** tbc[4],float**buc,float**bud,float* vx,float* vy,float* vz) 
{

//     input     buc
//     output    vx,vy,vz
//     bottom surface................................................
//     average of buc, bud
    	for(nv1 =1 ; nv1 <= iprv ; nv1++ )
	{
		buav[1][nv1] = (buc[1][nv1] + bud[1][nv1]) * 0.5;
		buav[2][nv1] = (buc[2][nv1] + bud[2][nv1]) * 0.5;
		buav[3][nv1] = (buc[3][nv1] + bud[3][nv1]) * 0.5;
	}

//     function in flux formula on surface
	vform2(imode,buav,fw,gw,hw);

	for(nv1 =1 ; nv1 <= iprv ; nv1++ )
	{
	//     function value on surface (f, g, h)
	      f1 = fw[1][nv1];
	      f2 = fw[2][nv1];
	      f3 = fw[3][nv1];
	      g1 = gw[1][nv1];
	      g2 = gw[2][nv1];
	      g3 = gw[3][nv1];
	      h1 = hw[1][nv1];
	      h2 = hw[2][nv1];
	      h3 = hw[3][nv1];
	//     normal vector inner product
	      v1 = f1*tbc[1][1][nv1] + g1*tbc[2][1][nv1] + h1*tbc[3][1][nv1];
	      v2 = f2*tbc[1][1][nv1] + g2*tbc[2][1][nv1] + h2*tbc[3][1][nv1];
	      v3 = f3*tbc[1][1][nv1] + g3*tbc[2][1][nv1] + h3*tbc[3][1][nv1];

	//     surface area
	      vx[nv1] = v1*sbc[nv1];
	      vy[nv1] = v2*sbc[nv1];
	      vz[nv1] = v3*sbc[nv1];
	}
}

void vform2(int imode ,float** buav ,float** fw ,float** gw ,float** hw) 
{
//     input       buav
//     output      fw,gw,hw

	if(imode == 1) // Divegence
	{
		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		      fw[1][nv1] = buav[1][nv1] ;
		      fw[2][nv1] = buav[1][nv1] ;
		      fw[3][nv1] = buav[1][nv1] ;
		      gw[1][nv1] = buav[2][nv1] ;
		      gw[2][nv1] = buav[2][nv1] ; 
		      gw[3][nv1] = buav[2][nv1] ;
		      hw[1][nv1] = buav[3][nv1] ;
		      hw[2][nv1] = buav[3][nv1] ; 
		      hw[3][nv1] = buav[3][nv1] ;
		}
	}

	if(imode == 2) // Curl
	{
		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		      fw[1][nv1] = 0 ;
		      fw[2][nv1] = -buav[3][nv1] ;
		      fw[3][nv1] = buav[2][nv1] ;
		      gw[1][nv1] = buav[3][nv1] ;
		      gw[2][nv1] = 0; 
		      gw[3][nv1] = -buav[1][nv1] ;
		      hw[1][nv1] = -buav[2][nv1] ;
		      hw[2][nv1] = buav[1][nv1] ; 
		      hw[3][nv1] = 0 ;
		}
	}

	if(imode == 3) // Gradient
	{
		for(nv1 =1 ; nv1 <= iprv ; nv1++ )
		{
		      fw[1][nv1] = buav[1][nv1] ;
		      fw[2][nv1] = 0 ;
		      fw[3][nv1] = 0 ;
		      gw[1][nv1] = 0 ;
		      gw[2][nv1] = buav[2][nv1] ; 
		      gw[3][nv1] = 0 ;
		      hw[1][nv1] = 0 ;
		      hw[2][nv1] = 0 ; 
		      hw[3][nv1] = buav[3][nv1] ;
		}
	}
}