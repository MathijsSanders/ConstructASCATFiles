//----------------------------------------------------------------
#ifndef OptionsH
#define OptionsH
//----------------------------------------------------------------
#include <string>
#include <vector>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class COptions
{
    private:
    public:
    
        bool            bCountDuplicates;

        int             iMinAlignmentScore;
        int             iMinBaseScore;
        int             nThreads;
        
        double          dGMAFThreshold;
        
        string          sSNPFileName;
        string          sOutputPrefix;
        
        vector<string>  vBamFileNames;
        vector<string>  vSampleNames;
        
        COptions(void)
        {
            bCountDuplicates=false;
            
            iMinAlignmentScore=40;
            iMinBaseScore=30;
            nThreads=1;
            
            dGMAFThreshold=0.01;
        }
};
//----------------------------------------------------------------
#endif
