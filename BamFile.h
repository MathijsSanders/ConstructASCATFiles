//----------------------------------------------------------------
#ifndef BamFileH
#define	BamFileH
//----------------------------------------------------------------
#include <string>
#include "bam.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
#define BAM_FSUPPLEMENTARY 2048

typedef struct
{
    int                     iTargetId;
    int                     iTargetPos;
    
    int                     nAlignments;
    const bam_pileup1_t *   pAlignments;
    
}PLP_DATA;
//----------------------------------------------------------------
class CBamFile
{
    private:
        
        bool            bInitialized;
        
        uint32_t        dwMinQual;
        
        string          sFileName;

        bamFile         pBamFile;
        bam_header_t *  pBamHeader;
        bam_index_t *   pBamIndex;
        bam_plp_t       pPileupObject;
        bam_iter_t      pBamRegion;
        
        PLP_DATA        pileupData;
        
        static int PileupCallback(void * pUserData,bam1_t * pAlignment);
        
    public:
        
        CBamFile(void);
        ~CBamFile(void);
        
        bool            Open(const string & sFileName,uint32_t dwMinQual,bool bCountDuplicates);
        void            Close(void);

        const char *    GetFileName(void);

        int             GetTargetCount(void);
        int             GetTargetLen(int iTargetId);
        const char *    GetTargetName(int iTargetId);
        PLP_DATA *      GetPileupData(void);

        bool            SetRegion(const string & sRegion);
        bool            PileupRegion(void);

        int             Read(bam1_t * pAlignment);
        int             ReadRegion(bam1_t * pAlignment);
};
//----------------------------------------------------------------
#endif
