//----------------------------------------------------------------
#ifndef VCFFileH
#define	VCFFileH
//----------------------------------------------------------------
#include <map>
#include <string>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
typedef struct
{
    int iRef;
    int iAlt;
    void * pUserData;

}VCF_ENTRY;
//----------------------------------------------------------------
class CVCFFile
{
    private:

    public:        
        
        int nEntries;
        map<string,map<int,VCF_ENTRY> > mEntries;
        
        bool Open(const string & sFileName,double dGMAFThreshold);
};
//----------------------------------------------------------------
#endif