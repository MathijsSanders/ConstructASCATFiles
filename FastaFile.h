//----------------------------------------------------------------
#ifndef FastaFileH
#define	FastaFileH
//----------------------------------------------------------------
#include <map>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
typedef struct
{
    int     iTargetLength;
    char *  pTarget;
    
}TARGET;
//----------------------------------------------------------------
class CFastaFile
{
    private:
        
    public: 

        char *              pSequence;
        map<string,TARGET>  targets;
        
        CFastaFile(void);
        ~CFastaFile(void);
        
        bool Open(const string & sFileName);
        void Clear(void);
};
//----------------------------------------------------------------
#endif