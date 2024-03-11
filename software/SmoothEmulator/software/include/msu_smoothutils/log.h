#ifndef __INCLUDE_LOG_H_
#define __INCLUDE_LOG_H_

#include <string>
#include <cstdlib>

using namespace std;
namespace NMSUUtils{

	class CLog{
	public:
		static bool INTERACTIVE;
		static string logfilename;
		static FILE *fptr;
		static void Init(char *logfilename_in);
		static void Init(string &logfilename_in);
		//static void Fatal(string &message);
		static void Info(string message);
		static void Fatal(string message);
		static const int CHARLENGTH=300;
		//static void Info(string &message);
		//static void Fatal(char *message);
		//static void Info(char *message);
	};
}

#endif
