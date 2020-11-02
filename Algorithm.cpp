//----------------------------------------------------------------
// Name        : Algorithm.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : 
//----------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "BamFile.h"
#include "ProgressBar.h"
#include "Algorithm.h"
//----------------------------------------------------------------
typedef struct
{
    int iBase;
    int iBaseScore;
    
}BASE;
//----------------------------------------------------------------
const bool bIsBase[]={0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0};
//----------------------------------------------------------------
void * CAlgorithm::Worker(void * pUserData)
{

    int                                     nAlignments;
    int                                     iAlignment;
    int                                     iMinBaseScore;
    int                                     iAlt;

    RESULTS_ENTRY                           result;
    
    const bam_pileup1_t *                   pAlignments;        
    PLP_DATA *                              pPileupData;
    map<int,VCF_ENTRY> *                    pTargetEntries;
    CAlgorithm *                            pAlgorithm;
    
    JOB                                     job;
    BASE                                    base;
    
    string                                  sBamFileName;
    string                                  sQueryName;
    
    unordered_map<string,BASE>              bases;

    CBamFile                                bamFile;
   
    map<int,VCF_ENTRY>::iterator            itTargetEntries;
    map<int,VCF_ENTRY>::iterator            itTargetEntriesEnd;
    unordered_map<string,BASE>::iterator    itBases;
    unordered_map<string,BASE>::iterator    itBasesEnd;
    
    //----------------------------------------------------------------
    //Lets go
    //----------------------------------------------------------------
    
    pAlgorithm=(CAlgorithm*)pUserData;
    
    iMinBaseScore=pAlgorithm->iMinBaseScore;
    
    pthread_mutex_lock(&pAlgorithm->mutex);
    
    for(;;)
    {
        //----------------------------------------------------------------
        //Get job from queue
        //----------------------------------------------------------------
        
        if(pAlgorithm->jobs.size()==0)
        {
            pthread_mutex_unlock(&pAlgorithm->mutex);
            pthread_exit(NULL);
        }
        
        job=pAlgorithm->jobs.front();
        pAlgorithm->jobs.pop_front();
        
        pthread_mutex_unlock(&pAlgorithm->mutex);
        
        //----------------------------------------------------------------
        //Check if target is in vcf file
        //----------------------------------------------------------------
        
        if(pAlgorithm->vcfFile.mEntries.count(job.sTarget)==0)
        {
            cerr << "Warning: No entries in vcf for target: " << job.sTarget << " (skipping)" << endl;
            continue;
        }
        
        pTargetEntries=&pAlgorithm->vcfFile.mEntries.at(job.sTarget);
        itTargetEntriesEnd=pTargetEntries->end();
        
        //----------------------------------------------------------------
        //Open bam file
        //----------------------------------------------------------------
        
        sBamFileName=pAlgorithm->vBamFileNames[job.iBamFile];
        
        if(bamFile.Open(sBamFileName,pAlgorithm->iMinAlignmentScore,pAlgorithm->bCountDuplicates)==false)
        {
            cerr << "Warning: Could not open bam file: " << sBamFileName << " (skipping)" << endl;
            continue;
        }
        
        //----------------------------------------------------------------
        //Set target region
        //----------------------------------------------------------------
        
        if(bamFile.SetRegion(job.sTarget)==false)
        {
            cerr << "Warning: Could not set target: " << job.sTarget << " (skipping)" << endl;
            continue;
        }
        
        //----------------------------------------------------------------
        //Pileup target
        //----------------------------------------------------------------
        
        pPileupData=bamFile.GetPileupData();
        
        while(bamFile.PileupRegion())
        {
            //----------------------------------------------------------------
            //Check if pileup position is in vcf file
            //----------------------------------------------------------------
            
            if((itTargetEntries=pTargetEntries->find(pPileupData->iTargetPos))==itTargetEntriesEnd)
            {
                continue;
            }
            
            //----------------------------------------------------------------
            //Iterate over alignments
            //----------------------------------------------------------------

            nAlignments=pPileupData->nAlignments;
            pAlignments=pPileupData->pAlignments;

            bases.clear();

            for(iAlignment=0;iAlignment<nAlignments;iAlignment++,pAlignments++)
            {
                //----------------------------------------------------------------
                //Skip alignment if it is a deletion on this pileup position
                //----------------------------------------------------------------

                if(pAlignments->is_del)     
                {
                    continue;
                }

                //----------------------------------------------------------------
                //Skip alignment if base on position is poor
                //----------------------------------------------------------------

                base.iBase=bam1_seqi(bam1_seq(pAlignments->b),pAlignments->qpos);
                base.iBaseScore=bam1_qual(pAlignments->b)[pAlignments->qpos];

                if(base.iBaseScore<iMinBaseScore)
                {
                    continue;
                }

                //----------------------------------------------------------------
                //Count fragments 
                //----------------------------------------------------------------

                sQueryName=string(bam1_qname(pAlignments->b));

                if((itBases=bases.find(sQueryName))==bases.end())
                {
                    bases[sQueryName]=base;
                }

                else
                {
                    if(itBases->second.iBaseScore<base.iBaseScore)
                    {
                        bases[sQueryName]=base;
                    }
                }
            }
            
            //----------------------------------------------------------------
            //Iterate over fragments and store results
            //----------------------------------------------------------------
            
            iAlt=itTargetEntries->second.iAlt;
            result={0,0};

            itBasesEnd=bases.end();
                
            for(itBases=bases.begin();itBases!=itBasesEnd;itBases++)
            {
                result.iCov+=bIsBase[itBases->second.iBase];    //Only count true bases
                result.iAltCov+=(itBases->second.iBase==iAlt);
            }
            
            ((RESULTS_ENTRY*)itTargetEntries->second.pUserData)[job.iBamFile]=result;
        }

        //----------------------------------------------------------------
        //Done
        //----------------------------------------------------------------
        
        pthread_mutex_lock(&pAlgorithm->mutex);
        
        pAlgorithm->nJobsDone++;
    }    
}
//----------------------------------------------------------------
bool CAlgorithm::InitJobs(const COptions & options)
{
    int             nSNPs;
    int             nSamples;
    int             nBamFiles;
    int             iBamFile;
    int             nTargets;
    int             iTarget;
    
    RESULTS_ENTRY * pResults;
    
    JOB             job;
    
    string          sBamFileName;
    
    CBamFile        bamFile;
    
    map<string,map<int,VCF_ENTRY> >::iterator   itChrom;
    map<string,map<int,VCF_ENTRY> >::iterator   itChromEnd;
    map<int,VCF_ENTRY>::iterator                itPos;
    map<int,VCF_ENTRY>::iterator                itPosEnd;
    
    //----------------------------------------------------------------
    //Init object members
    //----------------------------------------------------------------
    
    bCountDuplicates=options.bCountDuplicates;
        
    iMinAlignmentScore=options.iMinAlignmentScore;
    iMinBaseScore=options.iMinBaseScore;
    nJobsDone=0;
    
    vBamFileNames=options.vBamFileNames;
   
    //----------------------------------------------------------------
    //Open vcf file
    //----------------------------------------------------------------
    
    cerr << "Info: Open vcf file" << endl;
    
    if(vcfFile.Open(options.sSNPFileName,options.dGMAFThreshold)==false)
    {
        return false;
    }
    
    
    //----------------------------------------------------------------
    //Create results matrix
    //----------------------------------------------------------------

    cerr << "Info: Create results matrix" << endl;    
    
    nSNPs=vcfFile.nEntries;
    nSamples=options.vSampleNames.size();
    
    if(this->pResults)
    {
        delete [] this->pResults;
        this->pResults=NULL;
    }
    
    this->pResults=pResults=new RESULTS_ENTRY[nSNPs*nSamples]();
    
    itChromEnd=vcfFile.mEntries.end();
    
    for(itChrom=vcfFile.mEntries.begin();itChrom!=itChromEnd;itChrom++)
    {
        itPosEnd=itChrom->second.end();

        for(itPos=itChrom->second.begin();itPos!=itPosEnd;itPos++)
        {
            itPos->second.pUserData=pResults;
            pResults+=nSamples;
        }
    }
    
    //----------------------------------------------------------------
    //Create jobs
    //----------------------------------------------------------------
    
    cerr << "Info: Create jobs" << endl;
    
    nBamFiles=options.vBamFileNames.size();
    
    for(iBamFile=0;iBamFile<nBamFiles;iBamFile++)
    {
        job.iBamFile=iBamFile;
        
        sBamFileName=options.vBamFileNames[iBamFile];
        
        if(bamFile.Open(sBamFileName,options.iMinAlignmentScore,options.bCountDuplicates)==false)
        {
            cerr << "Warning: Could not open bam file: " << sBamFileName << " (skipping)" << endl;
            continue;
        }
        
        nTargets=bamFile.GetTargetCount();
        
        for(iTarget=0;iTarget<nTargets;iTarget++)
        {
            job.sTarget=string(bamFile.GetTargetName(iTarget));
            
            if(vcfFile.mEntries.count(job.sTarget)==0)
            {
                cerr << "Warning: No entries in vcf for target: " << job.sTarget << " (skipping)" << endl;
                continue;
            }
            
            jobs.push_back(job);
        }
    }
    
    if(jobs.size()==0)
    {
        cerr << "Error: No jobs to run! Please check that target names in vcf and bam files agree" << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return true;
}
//----------------------------------------------------------------
void CAlgorithm::RunJobs(const COptions & options)
{
    int         nJobs;
    int         iThread;
    void *      pThreadReturn;
    
    pthread_t * pThreads;    
    
    nJobs=jobs.size();
    
    //----------------------------------------------------------------
    //Create threads
    //----------------------------------------------------------------
    
    cerr << "Info: Create worker threads" << endl;
    
    pThreads=new pthread_t[options.nThreads];

    for(iThread=0;iThread<options.nThreads;iThread++)
    {
        pthread_create(&pThreads[iThread],NULL,Worker,(void*)this);
    }
    
    //----------------------------------------------------------------
    //Wait for threads to finish
    //----------------------------------------------------------------
    
    cerr << "Info: Wait for worker threads to finish" << endl;
    
    for(iThread=0;iThread<options.nThreads;iThread++)
    {
        while(pthread_tryjoin_np(pThreads[iThread],&pThreadReturn))
        {
            cerr << cProgressBar[(nJobsDone*100)/nJobs] << '\r' << flush;
            sleep(1);
        }
    }
    
    cerr << cProgressBar[(nJobsDone*100)/nJobs] << endl;
    
    delete [] pThreads; //Clean up threads
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
}
//----------------------------------------------------------------
bool CAlgorithm::WriteResults(const COptions & options)
{
    int             nSamples;
    int             iSample;
    int             iSNP;
    
    RESULTS_ENTRY * pResults;
    
    ofstream        covFile;
    ofstream        bafFile;

    map<string,map<int,VCF_ENTRY> >::iterator   itChrom;
    map<string,map<int,VCF_ENTRY> >::iterator   itChromEnd;
    map<int,VCF_ENTRY>::iterator                itPos;
    map<int,VCF_ENTRY>::iterator                itPosEnd;
    
    //----------------------------------------------------------------
    //Create output files
    //----------------------------------------------------------------

    cerr << "Info: Write results to files" << endl;    
    
    covFile.open((options.sOutputPrefix+string("_cov.txt")).c_str(),ofstream::out|ofstream::trunc);
    
    if(covFile.good()==false)
    {
        cerr << "Error: Could not create cov output file" << endl;
        return false;
    }
    
    bafFile.open((options.sOutputPrefix+string("_baf.txt")).c_str(),ofstream::out|ofstream::trunc);
    
    if(bafFile.good()==false)
    {
        cerr << "Error: Could not create baf output file" << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Write header
    //----------------------------------------------------------------
    
    covFile << "\tchrs\tpos";
    bafFile << "\tchrs\tpos";
    
    nSamples=options.vSampleNames.size();
    
    for(iSample=0;iSample<nSamples;iSample++)
    {
        covFile << '\t' << options.vSampleNames[iSample];
        bafFile << '\t' << options.vSampleNames[iSample];
    }

    covFile << '\n';
    bafFile << '\n';

    //----------------------------------------------------------------
    //Write data
    //----------------------------------------------------------------
    
    iSNP=1;
    
    itChromEnd=vcfFile.mEntries.end();
    
    for(itChrom=vcfFile.mEntries.begin();itChrom!=itChromEnd;itChrom++)
    {
        itPosEnd=itChrom->second.end();

        for(itPos=itChrom->second.begin();itPos!=itPosEnd;itPos++)
        {
            covFile << "SNP" << iSNP << '\t' << itChrom->first << '\t' << itPos->first+1;
            bafFile << "SNP" << iSNP << '\t' << itChrom->first << '\t' << itPos->first+1;
            
            pResults=(RESULTS_ENTRY*)itPos->second.pUserData;
        
            for(iSample=0;iSample<nSamples;iSample++)
            {
                if(pResults[iSample].iCov)  //Arrays are initialized to zero!
                {
                    covFile << '\t' << pResults[iSample].iCov;
                    bafFile << '\t' << double(pResults[iSample].iAltCov) / double(pResults[iSample].iCov);
                }
            }
            
            covFile << '\n';
            bafFile << '\n';
            
            iSNP++;
        }
    }
    
    //----------------------------------------------------------------
    //Something went wrong during writing?
    //----------------------------------------------------------------
    
    if(covFile.good()==false)
    {
        cerr << "Error: Could not write cov output file" << endl;
        return false;
    }
    
    if(bafFile.good()==false)
    {
        cerr << "Error: Could not write baf output file" << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return true;
}
//----------------------------------------------------------------
CAlgorithm::CAlgorithm(void)
{
    pResults=NULL;
    pthread_mutex_init(&mutex,NULL);
}
//----------------------------------------------------------------        
CAlgorithm::~CAlgorithm(void)
{
    if(pResults)
    {
        delete [] pResults;
        pResults=NULL;
    }
            
    pthread_mutex_destroy(&mutex);
}
//----------------------------------------------------------------        
bool CAlgorithm::Run(const COptions & options)
{
    //----------------------------------------------------------------
    //Init jobs
    //----------------------------------------------------------------

    if(InitJobs(options)==false)
    {
        return false;
    }
    
    //----------------------------------------------------------------
    //Run jobs
    //----------------------------------------------------------------
    
    RunJobs(options);
    
    //----------------------------------------------------------------
    //Write results 
    //----------------------------------------------------------------
    
    if(WriteResults(options)==false)
    {
        return false;
    }
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return true;
}
//----------------------------------------------------------------

