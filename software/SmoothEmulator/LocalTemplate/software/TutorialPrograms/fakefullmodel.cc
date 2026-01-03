#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include "msu_smoothutils/parametermap.h"
#include <filesystem>

using namespace std;
using namespace NMSUUtils;
int main(){
	double GSCALE=sqrt(3.0);
	double exp_unc;
	char dummy[400];
	unsigned int NObs,NPars;
	unsigned int NTrain,itrain,ic,ipar,maxrank=5,iy;
	double LAMBDA,xminread,xmaxread,sensitivityread,alpharead;
	vector<vector<double>> A;
	vector<vector<double>> xtrain, thetatrain,Ytrain;
	vector<double> theta,exptheta,SigmaY;
	vector<string> priortype;
	string expfilename,filename,command;
	vector<string> obsname;
	vector<string> parname;
   vector<double> alpha;
	vector<double> xmin,xmax,x0,Rgauss,thetatrue,ytrue,sensitivity;
	char parname_read[200],obsname_read[200],typeread[20];
	FILE *fptr;
	//NMSUUtils::Crandy randy(-time(NULL));
	NMSUUtils::Crandy randy(-123);

   // Read in prior info
	fptr=fopen("smooth_data/Info/prior_info.txt","r");
   NPars=0;
   fgets(dummy,400,fptr);
   do{
      fscanf(fptr,"%s",parname_read);
      if(!feof(fptr)){
         parname.push_back(parname_read);
         fscanf(fptr,"%s",typeread);
         priortype.push_back(typeread);
         fscanf(fptr,"%lf",&xminread);
         fscanf(fptr,"%lf",&xmaxread);
         fscanf(fptr,"%lf",&sensitivityread);
         xmin.push_back(xminread);
         xmax.push_back(xmaxread);
         sensitivity.push_back(sensitivityread);
         fgets(dummy,400,fptr);
         Rgauss.push_back(xmax[NPars]);
         x0.push_back(xmin[NPars]);
         NPars+=1;
      }
   }while(!feof(fptr));
   fclose(fptr);
   
   // Read in Observable Info
   fptr=fopen("smooth_data/Info/observable_info.txt","r");
   NObs=0;
   alpha.resize(0);
   obsname.resize(0);
   fscanf(fptr,"%s",obsname_read);
   do{
      while(obsname_read[0]=='#'){
         fgets(dummy,400,fptr);
         fscanf(fptr,"%s",obsname_read);
      }
      if(!feof(fptr)){
         obsname.push_back(obsname_read);
         NObs+=1;
         fscanf(fptr,"%lf",&alpharead);
         alpha.push_back(alpharead);
         fgets(dummy,400,fptr);
         fscanf(fptr,"%s",obsname_read);
      }
   }while(!feof(fptr));
   fclose(fptr);
   NObs=obsname.size();
   if(NObs!=obsname.size()){
      printf("alpha and obsname have different sizes, %zu != %zu\n",alpha.size(),obsname.size());
      exit(1);
   }
   
	printf("NPars=%d, NObs=%d\n",NPars,NObs);
	//printf("Enter Lambda for fake model: ");
	///scanf("%lf",&LAMBDA);
   LAMBDA=3.0;

	thetatrue.resize(NPars);
	ytrue.resize(NObs);
	SigmaY.resize(NObs);
	xmin.resize(NPars);
	xmax.resize(NPars);
	x0.resize(NPars);
	Rgauss.resize(NPars);
	obsname.resize(NObs);
	Ytrain.resize(NObs);

	// Find NTrain
	NTrain=0;
	bool existence;
	do{
		string filename="smooth_data/FullModelRuns/run"+to_string(NTrain);
		filesystem::path f{filename};
		existence=filesystem::exists(f);
		if(existence){
			NTrain+=1;
		}
	}while(existence);
	xtrain.resize(NTrain);
	thetatrain.resize(NTrain);
	for(itrain=0;itrain<NTrain;itrain++){
		xtrain[itrain].resize(NPars);
		thetatrain[itrain].resize(NPars);
	}
	CLog::Info("Creating Fake Model Data for "+to_string(NTrain)+" Training Pts\n");

	// Create smooth object
	NBandSmooth::CSmooth smooth(NPars,maxrank);
	// Create A for smooth object
	A.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		A[iy].resize(smooth.NCoefficients);
		for(ic=0;ic<smooth.NCoefficients;ic++){
			A[iy][ic]=100.0*randy.ran_gauss();
		}
	}
	// Create some arrays
   Ytrain.resize(NTrain);
   for(itrain=0;itrain<NTrain;itrain++){
      Ytrain[itrain].resize(NObs);
   }
	for(iy=0;iy<NObs;iy++){
		obsname[iy]="obs"+to_string(iy);
	}
	
   // Observable uncertainties
	// Set experimental value to theta=0.2
	fptr=fopen("smooth_data/Info/experimental_info.txt","w");
	for(iy=0;iy<NObs;iy++){
		ytrue[iy]=smooth.CalcY(A[iy],LAMBDA,thetatrue);
      exp_unc=2.0; // 2% of the strength 100
		fprintf(fptr,"%s\t%g\t%g 0.0\n",
		obsname[iy].c_str(),ytrue[iy],exp_unc);
	}
	fclose(fptr);
   
	// Read Training Points, Calc Y at points, and write observables for those points
	for(itrain=0;itrain<NTrain;itrain++){
		filename="smooth_data/FullModelRuns/run"+to_string(itrain)+"/model_parameters.txt";
		fptr=fopen(filename.c_str(),"r");
		for(ipar=0;ipar<NPars;ipar++){
			fscanf(fptr,"%s %lf",dummy,&xtrain[itrain][ipar]);
			if(priortype[ipar]=="uniform"){
            thetatrain[itrain][ipar]=-1.0+2.0*(xtrain[itrain][ipar]-xmin[ipar])/(xmax[ipar]-xmin[ipar]);
			}
			else{
            thetatrain[itrain][ipar]=(xtrain[itrain][ipar]-x0[ipar])/(Rgauss[ipar]*GSCALE);
			}
		}
		fclose(fptr);
		filename="smooth_data/FullModelRuns/run"+to_string(itrain)+"/obs.txt";
		fptr=fopen(filename.c_str(),"w");
		for(iy=0;iy<NObs;iy++){
			Ytrain[itrain][iy]=smooth.CalcY(A[iy],LAMBDA,thetatrain[itrain]);
         Ytrain[itrain][iy]+=alpha[iy]*SigmaY[iy]*randy.ran_gauss();
			fprintf(fptr,"%s %lf\n",obsname[iy].c_str(),Ytrain[itrain][iy]);
		}
		fclose(fptr);
	}
	//
	
	
	// Now make some data for later testing, not used for tuning
	vector<vector<double>> thetatest,xtest;
	unsigned int itest,Ntest=100;
	double y;
	
	thetatest.resize(Ntest);
	xtest.resize(Ntest);
	for(itest=0;itest<Ntest;itest++){
		thetatest[itest].resize(NPars);
		xtest[itest].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			if(priortype[ipar]=="uniform"){
				thetatest[itest][ipar]=-1.0+2.0*randy.ran();
            xtest[itest][ipar]=xmin[ipar]
            +0.5*(1.0+thetatest[itest][ipar])*(xmax[ipar]-xmin[ipar]);
			}
			else{
				thetatest[itest][ipar]=randy.ran_gauss()/GSCALE;
            if(itest<NTrain){
               thetatest[itest][ipar]=thetatrain[itest][ipar];
            }
            xtest[itest][ipar]=x0[ipar]+Rgauss[ipar]*thetatest[itest][ipar]*GSCALE;
			}
		}
	}
	
	
   for(itest=0;itest<Ntest;itest++){
      string testdirname="smooth_data/FullModelTestingRuns/run"+to_string(itest);
      command="mkdir -p "+testdirname;
      system(command.c_str());
      
      filename=testdirname+"/model_parameters.txt";
      fptr=fopen(filename.c_str(),"w");
      for(ipar=0;ipar<NPars;ipar++){
         fprintf(fptr,"%s %17.10e\n",parname[ipar].c_str(),xtest[itest][ipar]);
      }
      fclose(fptr);
      
      filename=testdirname+"/obs.txt";
      fptr=fopen(filename.c_str(),"w");
      for(iy=0;iy<NObs;iy++){
         y=smooth.CalcY(A[iy],LAMBDA,thetatest[itest]);
         fprintf(fptr,"%s %17.10e\n",obsname[iy].c_str(),y);
      }
      fclose(fptr);
		
	}

	// write experimental data (for use with mcmc)
	FILE *expfptr;
	expfilename="smooth_data/Info/experimental_info.txt";
	expfptr=fopen(expfilename.c_str(),"w");
	exptheta.resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
		exptheta[ipar]=0.2;
	}
	for(iy=0;iy<NObs;iy++){
		y=smooth.CalcY(A[iy],LAMBDA,exptheta);
		fprintf(expfptr,"%s  %lf  10.0 0.0\n",obsname[iy].c_str(),y);
	}
	
	fclose(expfptr);
	
	
	return 0;
}
