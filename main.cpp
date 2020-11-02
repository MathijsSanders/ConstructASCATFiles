//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : construct_ascat_files
//----------------------------------------------------------------
#include <iostream>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include "Algorithm.h"
//----------------------------------------------------------------
#define BAM_FILES           'b'
#define SAMPLE_NAMES        's'
#define SNP_FILE            'S'
#define OUTPUT_PREFIX       'o'
#define MIN_ALIGNMENT_SCORE 'a'
#define MIN_BASE_SCORE      'B'
#define THREADS             't'
#define GMAF_THRESHOLD      'g'
#define COUNT_DUPLICATES    'c'
#define HELP                'h'
#define SHORT_OPTIONS       "b:s:S:f:o:a:B:t:g:ch"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
struct option longOptions[] =
{
    {"bam-files",required_argument,NULL,BAM_FILES},
    {"sample-names",required_argument,NULL,SAMPLE_NAMES},
    {"snp-file",required_argument,NULL,SNP_FILE},
    {"output-prefix",required_argument,NULL,OUTPUT_PREFIX},
    {"min-alignment-score",required_argument,NULL,MIN_ALIGNMENT_SCORE},
    {"min-base-score",required_argument,NULL,MIN_BASE_SCORE},
    {"threads",required_argument,NULL,THREADS},
    {"gmaf-threshold",required_argument,NULL,GMAF_THRESHOLD},
    {"count-duplicates",required_argument,NULL,COUNT_DUPLICATES},
    {"help",no_argument,NULL,HELP},
    {0, 0, 0, 0}
};    
//----------------------------------------------------------------
vector<string> & Tokenize(const string & sString,const char * pDelimiter)
{
    int                     iStrLen;
    char *                  pString;
    char *                  pStringTemp;
    static vector<string>   vTokens;
    
    vTokens.clear();
    
    if((iStrLen=sString.size())>0)
    {
        pString=pStringTemp=(char*)memcpy(new char[iStrLen+1],sString.c_str(),iStrLen+1);

        while(pStringTemp)
        {
            vTokens.push_back(string(strsep(&pStringTemp,pDelimiter)));
        }
    
        delete [] pString;
    }
    
    return vTokens;
}
//----------------------------------------------------------------
int main(int argc, char** argv)
{
    bool        bShowHelp;

    int         iOption;
    int         iOptionIndex;
    
    COptions    options;
    CAlgorithm  algorithm;
    
    //----------------------------------------------------------------
    //Get options
    //----------------------------------------------------------------
    
    cerr << "Info: Get options" << endl;
    
    bShowHelp=(argc==1);
    
    while((iOption=getopt_long(argc,argv,SHORT_OPTIONS,longOptions,&iOptionIndex))>=0)
    {
        switch(iOption)
        {
            case BAM_FILES:
                
                options.vBamFileNames=Tokenize(string(optarg),",");
                break;
                
            case SAMPLE_NAMES:
                
                options.vSampleNames=Tokenize(string(optarg),",");
                break;
                
            case SNP_FILE:
                
                options.sSNPFileName=string(optarg);
                break;
                
            case OUTPUT_PREFIX:
                
                options.sOutputPrefix=string(optarg);
                break;
                
            case MIN_ALIGNMENT_SCORE:
                
                options.iMinAlignmentScore=atoi(optarg);
                break;
                
            case MIN_BASE_SCORE:
                
                options.iMinBaseScore=atoi(optarg);
                break;
                
            case THREADS:
                
                options.nThreads=atoi(optarg);
                break;
                
            case GMAF_THRESHOLD:
                
                options.dGMAFThreshold=atof(optarg);
                break;
                
            case COUNT_DUPLICATES:
                
                options.bCountDuplicates=true;
                break;
                
            default:
                
                bShowHelp=true;
                break;
        }
    }
    
    //----------------------------------------------------------------
    //Show help
    //----------------------------------------------------------------
    
    if(bShowHelp)
    {
        cerr << "construct_ascat_files [options]"                                                           << endl;
        cerr                                                                                                << endl;
        cerr << "-b --bam-files <text>          bam files (required comma separated)"                       << endl;
        cerr << "-s --sample-names <text>       sample names (required comma separated)"                    << endl;
        cerr << "-S --snp-file <text>           single vcf file with SNPs (required)"                       << endl;
        cerr << "-o --output-prefix <text>      output prefix (required)"                                   << endl;
        cerr << "-a --min-alignment-score <int> min alignment score (optional default = 40)"                << endl;                
        cerr << "-B --min-base-score <int>      min base score (optional default = 30)"                     << endl;
        cerr << "-t --threads <int>             number of threads to use (optional default = 1)"            << endl;
        cerr << "-g --gmaf-threshold <float>    genomic minor allele frequency (optional default = 0.01)"   << endl;
        cerr << "-c --count-duplicates          if specified duplicates are used (optional)"                << endl;
        cerr                                                                                                << endl;
        cerr << "(Number and order of bam files and sample names must agree!)"                              << endl;
        cerr                                                                                                << endl;
        
        return 1;
    }
    
    //----------------------------------------------------------------
    //Check options
    //----------------------------------------------------------------
    
    cerr << "Info: Check options" << endl;
    
    if(options.vBamFileNames.empty())
    {
        cerr << "Error: Please specify at least one bam file (--help)" << endl;
        return 1;
    }

    if(options.vSampleNames.size()!=options.vBamFileNames.size())
    {
        cerr << "Error: Number and order of bam files and sample names must agree (--help)" << endl;
        return 1;
    }
    
    if(options.sSNPFileName.empty())
    {
        cerr << "Error: Please specify a vcf file with SNPs (--help)" << endl;
        return 1;
    }
    
    if(options.sOutputPrefix.empty())
    {
        cerr << "Error: Please specify an output prefix (--help)" << endl;
        return 1;
    }
    
    if(options.nThreads<1)
    {
        options.nThreads=1;
    }
    
    if(options.dGMAFThreshold<0.0)
    {
        options.dGMAFThreshold=0.0;
    }
    
    if(options.dGMAFThreshold>1.0)
    {
        options.dGMAFThreshold=1.0;
    }
    
    //----------------------------------------------------------------
    //Lets go
    //----------------------------------------------------------------
    
    cerr << "Info: Run algorithm" << endl;
    
    if(algorithm.Run(options)==false)
    {
        return 1;
    }
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    cerr << "Info: Done" << endl;
    
    return 0;
}
//----------------------------------------------------------------

