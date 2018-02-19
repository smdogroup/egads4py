#ifndef SURREALS_H
#define SURREALS_H

#define SURREAL_TRAD

#if defined(SURREAL_TRAD)
#include "SurrealS_Trad.h"
#elif defined(SURREAL_LAZY)
#include "SurrealS_Lazy.h"
#else
#error "Please define SURREAL_TRAD or SURREAL_LAZY"
#endif

#undef SURREAL_TRAD

#endif // SURREALS_H
