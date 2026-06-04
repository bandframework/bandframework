#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>

using namespace std;

int main(){
	unsigned int NPars,NObs;
	unsigned int ipar,iy;
	double ALPHA,sensitivity=1.0;
	FILE *fptr;
	string parname, obsname;
   //printf("Enter the number of model parameters: ");
   //scanf("%d",&NPars);
   //printf("Enter the number of observables ");
   //scanf("%d",&NObs);
   NPars=6;
   NObs=10;
   ALPHA=0.0;
   
	fptr=fopen("smooth_data/Info/prior_info.txt","w");
	fprintf(fptr,"# par_name  dist_type  xmin/centroid  xmax/width  SensitivityScale\n");
	for(ipar=0;ipar<NPars;ipar++){
		parname="par"+to_string(ipar);
		if(ipar%2==0)
			fprintf(fptr,"%7s gaussian 0 100  %6.5f\n",parname.c_str(),sensitivity);
		else
			fprintf(fptr,"%7s uniform -100 100  %6.5f\n",parname.c_str(),sensitivity);
	}
	fclose(fptr);
   
	fptr=fopen("smooth_data/Info/observable_info.txt","w");
   fprintf(fptr,"# ObservableName   ALPHA\n");
	for(iy=0;iy<NObs;iy++){
		obsname="obs"+to_string(iy);
		fprintf(fptr,"%7s %g\n",obsname.c_str(),ALPHA);
	}
	fclose(fptr);
	
  return 0;
}


