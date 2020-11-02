//----------------------------------------------------------------
#ifndef AlgorithmH
#define AlgorithmH
//----------------------------------------------------------------
#include <deque>
#include "VCFFile.h"
#include "Options.h"
//----------------------------------------------------------------
typedef struct
{
    int     iBamFile;
    string  sTarget;

}JOB;
//----------------------------------------------------------------
typedef struct
{
    int iCov;
    int iAltCov;
    
}RESULTS_ENTRY;
//----------------------------------------------------------------
class CAlgorithm
{
    private:
        
        bool            bCountDuplicates;
        
        int             iMinAlignmentScore;
        int             iMinBaseScore;
        int             nJobsDone;
        
        RESULTS_ENTRY * pResults;

        pthread_mutex_t mutex;
        
        vector<string>  vBamFileNames;
        deque<JOB>      jobs;
        
        CVCFFile        vcfFile;
        
        static void * Worker(void * pUserData);
        
        bool InitJobs(const COptions & options);
        void RunJobs(const COptions & options);
        bool WriteResults(const COptions & options);
        
    public:
        
        CAlgorithm(void);
        ~CAlgorithm(void);
        
        bool Run(const COptions & options);
};
//----------------------------------------------------------------
#endif 
