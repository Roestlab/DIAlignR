#include "affinealignobj.h"


std::ostream& operator<<(std::ostream& out, const TracebackType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(SS);
        PROCESS_VAL(DM);
        PROCESS_VAL(DA);
        PROCESS_VAL(DB);
        PROCESS_VAL(TM);
        PROCESS_VAL(TA);
        PROCESS_VAL(TB);
        PROCESS_VAL(LM);
        PROCESS_VAL(LA);
        PROCESS_VAL(LB);
    }
#undef PROCESS_VAL
    return out << s;
}
