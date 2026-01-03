#include "msu_smooth/smooth.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CSmooth::CSmooth(unsigned int NPars_Set,unsigned int maxrank_set){
	NPars=NPars_Set;
	MaxRank=maxrank_set;
	InitArrays();
}

CSmooth::CSmooth(){
	CparameterMap *parmap;
   parmap=new CparameterMap();
	parmap->ReadParsFromFile("smooth_data/Options/emulator_options.txt");
	NPars=parmap->getI("SmoothEmulator_NPars",0);
	MaxRank=parmap->getI("Smooth_MAXRANK",5);
	if(MaxRank>5){
		CLog::Info("Inside CSmooth::InitArrays(), MaxRank="+to_string(MaxRank)+" is too big, being reset to 5\n");
		MaxRank=5;
	}
	InitArrays();
}

void CSmooth::InitArrays(){
	unsigned int ic,j,isame,ir;
	vector<unsigned int> countsame;
	vector<unsigned int> dummy;
	vector<unsigned int> i;
	factorial.resize(MaxRank+1);
	factorial[0]=factorial[1]=1;
	for(j=2;j<=MaxRank;j++)
		factorial[j]=j*factorial[j-1];

	i.resize(MaxRank+1);
	countsame.resize(MaxRank);

	ic=0;
	rank.resize(ic+1);
	rank[ic]=0;
	dupfactor.push_back(1.0);
	IPar.resize(ic+1);
	IPar[0].resize(rank[ic]);

	ic+=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		dupfactor.push_back(0);
		IPar.resize(ic+1);
		rank.resize(ic+1);
		rank[ic]=1;
		IPar[ic].resize(rank[ic]);
		for(ir=0;ir<rank[ic];ir++){
			IPar[ic][ir]=i[ir];
		}
		for(j=0;j<rank[ic];j++)
			countsame[j]=0;
		isame=0;
		countsame[isame]=1;
		for(j=1;j<rank[ic];j++){
			if(i[j]!=i[j-1])
				isame+=1;
			countsame[isame]+=1;
		}
		dupfactor[ic]=factorial[rank[ic]];

		for(j=0;j<rank[ic];j++)
			dupfactor[ic]/=factorial[countsame[j]];
		ic+=1;
	}

	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			dupfactor.push_back(0);
			IPar.resize(ic+1);
			rank.resize(ic+1);
			rank[ic]=2;
			IPar[ic].resize(rank[ic]);
			for(ir=0;ir<rank[ic];ir++){
				IPar[ic][ir]=i[ir];
			}
			for(j=0;j<rank[ic];j++)
				countsame[j]=0;
			isame=0;
			countsame[isame]=1;
			for(j=1;j<rank[ic];j++){
				if(i[j]!=i[j-1])
					isame+=1;
				countsame[isame]+=1;
			}
			dupfactor[ic]=factorial[rank[ic]];
			for(j=0;j<rank[ic];j++)
				dupfactor[ic]/=factorial[countsame[j]];
			ic+=1;
		}
	}

	if(MaxRank>=3){
		for(i[0]=0;i[0]<NPars;i[0]++){
			for(i[1]=0;i[1]<=i[0];i[1]++){
				for(i[2]=0;i[2]<=i[1];i[2]++){
					dupfactor.push_back(0);
					IPar.resize(ic+1);
					rank.resize(ic+1);
					rank[ic]=3;
					IPar[ic].resize(rank[ic]);
					for(ir=0;ir<rank[ic];ir++){
						IPar[ic][ir]=i[ir];
					}
					for(j=0;j<rank[ic];j++)
						countsame[j]=0;
					isame=0;
					countsame[isame]=1;
					for(j=1;j<rank[ic];j++){
						if(i[j]!=i[j-1])
							isame+=1;
						countsame[isame]+=1;
					}
					dupfactor[ic]=factorial[rank[ic]];
					for(j=0;j<rank[ic];j++)
						dupfactor[ic]/=factorial[countsame[j]];
					ic+=1;
				}
			}
		}
	}

	if(MaxRank>=4){
		for(i[0]=0;i[0]<NPars;i[0]++){
			for(i[1]=0;i[1]<=i[0];i[1]++){
				for(i[2]=0;i[2]<=i[1];i[2]++){
					for(i[3]=0;i[3]<=i[2];i[3]++){
						dupfactor.push_back(0);
						IPar.resize(ic+1);
						rank.resize(ic+1);
						rank[ic]=4;
						IPar[ic].resize(rank[ic]);
						for(ir=0;ir<rank[ic];ir++){
							IPar[ic][ir]=i[ir];
						}
						for(j=0;j<rank[ic];j++)
							countsame[j]=0;
						isame=0;
						countsame[isame]=1;
						for(j=1;j<rank[ic];j++){
							if(i[j]!=i[j-1])
								isame+=1;
							countsame[isame]+=1;
						}
						dupfactor[ic]=factorial[rank[ic]];
						for(j=0;j<rank[ic];j++)
							dupfactor[ic]/=factorial[countsame[j]];
						ic+=1;
					}
				}
			}
		}
	}

	if(MaxRank>=5){
		for(i[0]=0;i[0]<NPars;i[0]++){
			for(i[1]=0;i[1]<=i[0];i[1]++){
				for(i[2]=0;i[2]<=i[1];i[2]++){
					for(i[3]=0;i[3]<=i[2];i[3]++){
						for(i[4]=0;i[4]<=i[3];i[4]++){
							dupfactor.push_back(0);
							IPar.resize(ic+1);
							rank.resize(ic+1);
							rank[ic]=5;
							IPar[ic].resize(rank[ic]);
							for(ir=0;ir<rank[ic];ir++){
								IPar[ic][ir]=i[ir];
							}
							for(j=0;j<rank[ic];j++)
								countsame[j]=0;
							isame=0;
							countsame[isame]=1;
							for(j=1;j<rank[ic];j++){
								if(i[j]!=i[j-1])
									isame+=1;
								countsame[isame]+=1;
							}
							dupfactor[ic]=factorial[rank[ic]];
							for(j=0;j<rank[ic];j++)
								dupfactor[ic]/=factorial[countsame[j]];
							ic+=1;
						}
					}
				}
			}
		}
	}
	NCoefficients=ic;
	if(NCoefficients!=IPar.size()){
		CLog::Fatal("size mismatch in CSmooth::InitArrays()\n");
	}
}

double CSmooth::GetRFactor(double LAMBDA,vector<double> &theta){
	unsigned int ipar,NPars=theta.size();
	double r2=0.0,answer;
	for(ipar=0;ipar<NPars;ipar++)
		r2+=theta[ipar]*theta[ipar];
	answer=exp(-0.5*r2/(LAMBDA*LAMBDA));
	return answer;
}

double CSmooth::GetT(unsigned int ic,double LAMBDA,vector<double> &theta){
	unsigned int ir;
	double answer,rfactor=1.0;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=rfactor*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
	for(ir=0;ir<rank[ic];ir++){
		answer*=theta[IPar[ic][ir]]/LAMBDA;
	}
	return answer;
}

void CSmooth::Copy(CSmooth *smooth){
	MaxRank=smooth->MaxRank;
	NPars=smooth->NPars;
	NCoefficients=smooth->NCoefficients;
	dupfactor=smooth->dupfactor;
	IPar=smooth->IPar;
	rank=smooth->rank;
	factorial=smooth->factorial;
}
