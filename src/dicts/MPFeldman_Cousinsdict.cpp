// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdIMPFeldman_Cousinsdict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "MPFeldman_Cousins.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_MPFeldman_Cousins(void *p = 0);
   static void *newArray_MPFeldman_Cousins(Long_t size, void *p);
   static void delete_MPFeldman_Cousins(void *p);
   static void deleteArray_MPFeldman_Cousins(void *p);
   static void destruct_MPFeldman_Cousins(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MPFeldman_Cousins*)
   {
      ::MPFeldman_Cousins *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MPFeldman_Cousins >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MPFeldman_Cousins", ::MPFeldman_Cousins::Class_Version(), "MPFeldman_Cousins.hh", 21,
                  typeid(::MPFeldman_Cousins), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MPFeldman_Cousins::Dictionary, isa_proxy, 4,
                  sizeof(::MPFeldman_Cousins) );
      instance.SetNew(&new_MPFeldman_Cousins);
      instance.SetNewArray(&newArray_MPFeldman_Cousins);
      instance.SetDelete(&delete_MPFeldman_Cousins);
      instance.SetDeleteArray(&deleteArray_MPFeldman_Cousins);
      instance.SetDestructor(&destruct_MPFeldman_Cousins);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MPFeldman_Cousins*)
   {
      return GenerateInitInstanceLocal((::MPFeldman_Cousins*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MPFeldman_Cousins*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MPFeldman_Cousins::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MPFeldman_Cousins::Class_Name()
{
   return "MPFeldman_Cousins";
}

//______________________________________________________________________________
const char *MPFeldman_Cousins::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MPFeldman_Cousins*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MPFeldman_Cousins::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MPFeldman_Cousins*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MPFeldman_Cousins::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MPFeldman_Cousins*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MPFeldman_Cousins::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MPFeldman_Cousins*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MPFeldman_Cousins::Streamer(TBuffer &R__b)
{
   // Stream an object of class MPFeldman_Cousins.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MPFeldman_Cousins::Class(),this);
   } else {
      R__b.WriteClassBuffer(MPFeldman_Cousins::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MPFeldman_Cousins(void *p) {
      return  p ? new(p) ::MPFeldman_Cousins : new ::MPFeldman_Cousins;
   }
   static void *newArray_MPFeldman_Cousins(Long_t nElements, void *p) {
      return p ? new(p) ::MPFeldman_Cousins[nElements] : new ::MPFeldman_Cousins[nElements];
   }
   // Wrapper around operator delete
   static void delete_MPFeldman_Cousins(void *p) {
      delete ((::MPFeldman_Cousins*)p);
   }
   static void deleteArray_MPFeldman_Cousins(void *p) {
      delete [] ((::MPFeldman_Cousins*)p);
   }
   static void destruct_MPFeldman_Cousins(void *p) {
      typedef ::MPFeldman_Cousins current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MPFeldman_Cousins

namespace ROOT {
   static TClass *vectorlEbeltgR_Dictionary();
   static void vectorlEbeltgR_TClassManip(TClass*);
   static void *new_vectorlEbeltgR(void *p = 0);
   static void *newArray_vectorlEbeltgR(Long_t size, void *p);
   static void delete_vectorlEbeltgR(void *p);
   static void deleteArray_vectorlEbeltgR(void *p);
   static void destruct_vectorlEbeltgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<belt>*)
   {
      vector<belt> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<belt>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<belt>", -2, "vector", 216,
                  typeid(vector<belt>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEbeltgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<belt>) );
      instance.SetNew(&new_vectorlEbeltgR);
      instance.SetNewArray(&newArray_vectorlEbeltgR);
      instance.SetDelete(&delete_vectorlEbeltgR);
      instance.SetDeleteArray(&deleteArray_vectorlEbeltgR);
      instance.SetDestructor(&destruct_vectorlEbeltgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<belt> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<belt>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEbeltgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<belt>*)0x0)->GetClass();
      vectorlEbeltgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEbeltgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEbeltgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<belt> : new vector<belt>;
   }
   static void *newArray_vectorlEbeltgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<belt>[nElements] : new vector<belt>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEbeltgR(void *p) {
      delete ((vector<belt>*)p);
   }
   static void deleteArray_vectorlEbeltgR(void *p) {
      delete [] ((vector<belt>*)p);
   }
   static void destruct_vectorlEbeltgR(void *p) {
      typedef vector<belt> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<belt>

namespace {
  void TriggerDictionaryInitialization_MPFeldman_Cousinsdict_Impl() {
    static const char* headers[] = {
"MPFeldman_Cousins.hh",
0
    };
    static const char* includePaths[] = {
"/home/shoram/Work/ROOT/install/include",
"/home/shoram/Work/Diploma_Thesis/MPFC/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MPFeldman_Cousinsdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$MPFeldman_Cousins.hh")))  MPFeldman_Cousins;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MPFeldman_Cousinsdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "MPFeldman_Cousins.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"MPFeldman_Cousins", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MPFeldman_Cousinsdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MPFeldman_Cousinsdict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MPFeldman_Cousinsdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MPFeldman_Cousinsdict() {
  TriggerDictionaryInitialization_MPFeldman_Cousinsdict_Impl();
}
