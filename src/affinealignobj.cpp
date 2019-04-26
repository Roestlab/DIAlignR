#include "affinealignobj.h"

// This function overloads << to display TracebackType.
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

// This function converts TracebackType Enum to characters.
std::vector<char> EnumToChar(std::vector<TracebackType> v) {
  std::vector<char> nv(v.begin(), v.end());
  for (auto& it : nv) it += 48;
  //std::vector<char> nv(v.size());
  //std::transform(v.begin(), v.end(), nv.begin(), [](TracebackType c) -> char {return (c + 48); });
  return nv;
}
