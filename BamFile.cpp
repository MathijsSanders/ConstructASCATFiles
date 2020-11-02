//----------------------------------------------------------------
// Name        : BamFile.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include "BamFile.h"
//----------------------------------------------------------------
int CBamFile::PileupCallback(void * pUserData,bam1_t * pAlignment)
{
    return ((CBamFile*)pUserData)->ReadRegion(pAlignment);
}
//----------------------------------------------------------------
CBamFile::CBamFile(void) : bInitialized(false) , pBamFile(NULL) , pBamHeader(NULL) , pBamIndex(NULL) ,pPileupObject(NULL) , pBamRegion(NULL)
{
}
//----------------------------------------------------------------        
CBamFile::~CBamFile(void)
{
    Close();
}
//----------------------------------------------------------------        
bool CBamFile::Open(const string & sFileName,uint32_t dwMinQual,bool bCountDuplicates)
{
    int             iMask;

    bamFile         pBamFile;
    bam_header_t *  pBamHeader;
    bam_index_t *   pBamIndex;
    bam_plp_t       pPileupObject;

    //----------------------------------------------------------------
    //Open bam file
    //----------------------------------------------------------------
    
    pBamFile=bam_open(sFileName.c_str(),"r");
            
    if(pBamFile==NULL)
    {
        cerr << "Error: Could not open bam file" << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Read bam header
    //----------------------------------------------------------------
    
    pBamHeader=bam_header_read(pBamFile);
            
    if(pBamHeader==NULL)
    {
        cerr << "Error: Could not read header bam file" << endl;
        bam_close(pBamFile);
        return false;
    }
    
    //----------------------------------------------------------------
    //Load bam index
    //----------------------------------------------------------------
    
    pBamIndex=bam_index_load(sFileName.c_str());

    if(pBamIndex==NULL)
    {
        cerr << "Warning: Could not load index bam file (Creating)" << endl;

        //----------------------------------------------------------------
        //Build bam index if it does not exist
        //----------------------------------------------------------------
        
        bam_index_build(sFileName.c_str());
        
        pBamIndex=bam_index_load(sFileName.c_str());
        
        if(pBamIndex==NULL)
        {
            cerr << "Error: Could not build index bam file" << endl;
            bam_header_destroy(pBamHeader);
            bam_close(pBamFile);
            return false;
        }
    }

    pPileupObject=bam_plp_init(PileupCallback,(void*)this);
    
    if(pPileupObject==NULL)
    {
        cerr << "Error: Could not init pileup object bam file" << endl;
        bam_index_destroy(pBamIndex);
        bam_header_destroy(pBamHeader);
        bam_close(pBamFile);
        return false;
    }
    
    iMask=BAM_FUNMAP|BAM_FQCFAIL|BAM_FSUPPLEMENTARY;
    
    if(bCountDuplicates==false)
    {
        iMask|=BAM_FDUP;
    }
    
    bam_plp_set_maxcnt(pPileupObject,2147483647);
    bam_plp_set_mask(pPileupObject,iMask);
    
    //----------------------------------------------------------------
    //Close previous bam file
    //----------------------------------------------------------------
    
    Close();
    
    //----------------------------------------------------------------
    //Save results
    //----------------------------------------------------------------
    
    this->dwMinQual=dwMinQual;
    this->sFileName=sFileName;
    this->pBamFile=pBamFile;
    this->pBamHeader=pBamHeader;
    this->pBamIndex=pBamIndex;
    this->pPileupObject=pPileupObject;
    
    bInitialized=true;
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return true;
}
//----------------------------------------------------------------                
void CBamFile::Close(void)
{
    bInitialized=false;

    if(pBamRegion)
    {
        bam_iter_destroy(pBamRegion);
        pBamRegion=NULL;
    }
    
    if(pPileupObject)
    {
        bam_plp_destroy(pPileupObject);
        pPileupObject=NULL;
    }
    
    if(pBamIndex)
    {
        bam_index_destroy(pBamIndex);
        pBamIndex=NULL;
    }

    if(pBamHeader)
    {
        bam_header_destroy(pBamHeader);
        pBamHeader=NULL;
    }

    if(pBamFile)
    {
        bam_close(pBamFile);
        pBamFile=NULL;
    }
    
    sFileName.clear();
}
//----------------------------------------------------------------
const char * CBamFile::GetFileName(void)
{
    return sFileName.c_str();
}
//----------------------------------------------------------------
int CBamFile::GetTargetCount(void)
{
    if(bInitialized==false)
    {
        return 0;
    }
    
    return pBamHeader->n_targets;
}
//----------------------------------------------------------------
int CBamFile::GetTargetLen(int iTargetId)
{
    if(bInitialized==false || iTargetId<0 || iTargetId>=pBamHeader->n_targets)
    {
        return 0;
    }
    
    return pBamHeader->target_len[iTargetId];
}
//----------------------------------------------------------------        
const char * CBamFile::GetTargetName(int iTargetId)
{
    if(bInitialized==false || iTargetId<0 || iTargetId>=pBamHeader->n_targets)
    {
        return NULL;
    }

    return pBamHeader->target_name[iTargetId];
}
//----------------------------------------------------------------
PLP_DATA * CBamFile::GetPileupData(void)
{
    if(bInitialized==false)
    {
        return NULL;
    }

    return &pileupData;
}
//----------------------------------------------------------------
bool CBamFile::SetRegion(const string & sRegion)
{
    int         iTargetId;
    int         iBegin;
    int         iEnd;
    bam_iter_t  pBamRegion;
    
    if(bInitialized==false)
    {
        return false;
    }

    if(bam_parse_region(pBamHeader,sRegion.c_str(),&iTargetId,&iBegin,&iEnd)==-1)
    {
        return false;
    }
    
    pBamRegion=bam_iter_query(pBamIndex,iTargetId,iBegin,iEnd);    
    
    if(pBamRegion==NULL)
    {
        return false;
    }
    
    bam_plp_reset(pPileupObject);
    
    if(this->pBamRegion)
    {
        bam_iter_destroy(this->pBamRegion);
    }
    
    this->pBamRegion=pBamRegion;
    
    return true;
}
//----------------------------------------------------------------
bool CBamFile::PileupRegion(void)
{
    return (pileupData.pAlignments=bam_plp_auto(pPileupObject,&pileupData.iTargetId,&pileupData.iTargetPos,&pileupData.nAlignments))!=NULL;
}
//----------------------------------------------------------------
int CBamFile::Read(bam1_t * pAlignment)
{
    int iRet;

    for(iRet=bam_read1(pBamFile,pAlignment);iRet>=0 && (pAlignment->core.tid<0 || pAlignment->core.qual<dwMinQual);iRet=bam_read1(pBamFile,pAlignment));

    return iRet;
}
//----------------------------------------------------------------
int CBamFile::ReadRegion(bam1_t * pAlignment)
{
    int iRet;
    
    for(iRet=bam_iter_read(pBamFile,pBamRegion,pAlignment);iRet>=0 && (pAlignment->core.tid<0 || pAlignment->core.qual<dwMinQual);iRet=bam_iter_read(pBamFile,pBamRegion,pAlignment));

    return iRet;
}
//----------------------------------------------------------------

