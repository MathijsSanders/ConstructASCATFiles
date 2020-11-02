//----------------------------------------------------------------
// Name        : VCFFile.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "VCFFile.h"
//----------------------------------------------------------------
#define A   1
#define C   2
#define G   4
#define T   8
#define N   15
//----------------------------------------------------------------
const int iAscii2Index[256]=
{
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,    
    N,A,N,C,N,N,N,G,N,N,N,N,N,N,N,N,
    N,N,N,N,T,N,N,N,N,N,N,N,N,N,N,N,
    N,A,N,C,N,N,N,G,N,N,N,N,N,N,N,N,
    N,N,N,N,T,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,    
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N
};
//----------------------------------------------------------------
bool CVCFFile::Open(const string & sVCFFileName,double dGMAFThreshold=0.01)
{
    int             iRef;
    int             iAlt;

    double          dGMAF;

    char *          pLine;
    char *          pGMAF;
    char *          pChrom;
    char *          pPos;
    char *          pRef;
    char *          pAlt;
    
    string          sLine;
    ifstream        sVCFFile;
    
    //----------------------------------------------------------------
    //Clear previous hash map if any
    //----------------------------------------------------------------
    
    nEntries=0;
    mEntries.clear();
    
    //----------------------------------------------------------------
    //Open vcf file
    //----------------------------------------------------------------
    
    sVCFFile.open(sVCFFileName.c_str(),ifstream::in);
    
    if(sVCFFile.good()==false)
    {
        cerr << "Error: Could not open the vcf file!" << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Read first line from vcf file (fist file should contain vcf file format)
    //----------------------------------------------------------------
    
    getline(sVCFFile,sLine);
    
    if((sVCFFile.good()==false)||(sLine.find_first_of(string("##fileformat=VCFv4"),0)==string::npos))
    {
        cerr << "Error: Specified file not a valid vcf file!" << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Skip header lines
    //----------------------------------------------------------------
    
    for(getline(sVCFFile,sLine);sVCFFile.good() && sLine[0]=='#';getline(sVCFFile,sLine));
    
    //----------------------------------------------------------------
    //Read through the file
    //----------------------------------------------------------------
    
    for(;sVCFFile.good();getline(sVCFFile,sLine))
    {
        //----------------------------------------------------------------
        //Skip empty line
        //----------------------------------------------------------------
        
        if(sLine.empty())
        {
            continue;
        }
        
        //----------------------------------------------------------------
        //Tokenize
        //----------------------------------------------------------------
        
        pLine=(char*)sLine.c_str();
        
        if(dGMAFThreshold>0.0)
        {
            if((pGMAF=strstr(pLine,"GMAF="))==NULL)
            {
                continue;
            }
        
            dGMAF=atof(pGMAF+5);
        
            if(dGMAF<dGMAFThreshold)
            {
                continue;
            }
        }

        pChrom=strsep(&pLine,"\t");
        pPos=strsep(&pLine,"\t");
        strsep(&pLine,"\t");
        pRef=strsep(&pLine,"\t");
        pAlt=strsep(&pLine,"\t");
        
        //----------------------------------------------------------------
        //Not enough columns?
        //----------------------------------------------------------------
        
        if(pAlt==NULL)
        {
            cerr << "Error: Not enough columns in VCF file!" << endl;
            return false;
        }
        
        //----------------------------------------------------------------
        //Add entry 
        //----------------------------------------------------------------
        
        if(strlen(pRef)==1 && strlen(pAlt)==1 && (iRef=iAscii2Index[int(pRef[0])])!=N && (iAlt=iAscii2Index[int(pAlt[0])])!=N)
        {
            nEntries++;
            mEntries[string(pChrom)][atoi(pPos)-1]={iRef,iAlt,NULL};
        }
    }
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return true;
}
//----------------------------------------------------------------


