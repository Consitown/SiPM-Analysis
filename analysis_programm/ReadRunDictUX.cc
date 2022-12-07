// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME ReadRunDictUX
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

// Header files passed as explicit arguments
#include "ReadRun.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_ReadRun(void *p = nullptr);
   static void *newArray_ReadRun(Long_t size, void *p);
   static void delete_ReadRun(void *p);
   static void deleteArray_ReadRun(void *p);
   static void destruct_ReadRun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ReadRun*)
   {
      ::ReadRun *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ReadRun >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ReadRun", ::ReadRun::Class_Version(), "ReadRun.h", 59,
                  typeid(::ReadRun), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ReadRun::Dictionary, isa_proxy, 4,
                  sizeof(::ReadRun) );
      instance.SetNew(&new_ReadRun);
      instance.SetNewArray(&newArray_ReadRun);
      instance.SetDelete(&delete_ReadRun);
      instance.SetDeleteArray(&deleteArray_ReadRun);
      instance.SetDestructor(&destruct_ReadRun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ReadRun*)
   {
      return GenerateInitInstanceLocal((::ReadRun*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ReadRun*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *Fitf_Dictionary();
   static void Fitf_TClassManip(TClass*);
   static void *new_Fitf(void *p = nullptr);
   static void *newArray_Fitf(Long_t size, void *p);
   static void delete_Fitf(void *p);
   static void deleteArray_Fitf(void *p);
   static void destruct_Fitf(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitf*)
   {
      ::Fitf *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Fitf));
      static ::ROOT::TGenericClassInfo 
         instance("Fitf", "ReadRun.h", 303,
                  typeid(::Fitf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Fitf_Dictionary, isa_proxy, 4,
                  sizeof(::Fitf) );
      instance.SetNew(&new_Fitf);
      instance.SetNewArray(&newArray_Fitf);
      instance.SetDelete(&delete_Fitf);
      instance.SetDeleteArray(&deleteArray_Fitf);
      instance.SetDestructor(&destruct_Fitf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitf*)
   {
      return GenerateInitInstanceLocal((::Fitf*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitf*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Fitf_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Fitf*)nullptr)->GetClass();
      Fitf_TClassManip(theClass);
   return theClass;
   }

   static void Fitf_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Fitf_biased_Dictionary();
   static void Fitf_biased_TClassManip(TClass*);
   static void *new_Fitf_biased(void *p = nullptr);
   static void *newArray_Fitf_biased(Long_t size, void *p);
   static void delete_Fitf_biased(void *p);
   static void deleteArray_Fitf_biased(void *p);
   static void destruct_Fitf_biased(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitf_biased*)
   {
      ::Fitf_biased *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Fitf_biased));
      static ::ROOT::TGenericClassInfo 
         instance("Fitf_biased", "ReadRun.h", 421,
                  typeid(::Fitf_biased), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Fitf_biased_Dictionary, isa_proxy, 4,
                  sizeof(::Fitf_biased) );
      instance.SetNew(&new_Fitf_biased);
      instance.SetNewArray(&newArray_Fitf_biased);
      instance.SetDelete(&delete_Fitf_biased);
      instance.SetDeleteArray(&deleteArray_Fitf_biased);
      instance.SetDestructor(&destruct_Fitf_biased);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitf_biased*)
   {
      return GenerateInitInstanceLocal((::Fitf_biased*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitf_biased*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Fitf_biased_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Fitf_biased*)nullptr)->GetClass();
      Fitf_biased_TClassManip(theClass);
   return theClass;
   }

   static void Fitf_biased_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Fitf_PMT_Dictionary();
   static void Fitf_PMT_TClassManip(TClass*);
   static void *new_Fitf_PMT(void *p = nullptr);
   static void *newArray_Fitf_PMT(Long_t size, void *p);
   static void delete_Fitf_PMT(void *p);
   static void deleteArray_Fitf_PMT(void *p);
   static void destruct_Fitf_PMT(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitf_PMT*)
   {
      ::Fitf_PMT *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Fitf_PMT));
      static ::ROOT::TGenericClassInfo 
         instance("Fitf_PMT", "ReadRun.h", 473,
                  typeid(::Fitf_PMT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Fitf_PMT_Dictionary, isa_proxy, 4,
                  sizeof(::Fitf_PMT) );
      instance.SetNew(&new_Fitf_PMT);
      instance.SetNewArray(&newArray_Fitf_PMT);
      instance.SetDelete(&delete_Fitf_PMT);
      instance.SetDeleteArray(&deleteArray_Fitf_PMT);
      instance.SetDestructor(&destruct_Fitf_PMT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitf_PMT*)
   {
      return GenerateInitInstanceLocal((::Fitf_PMT*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitf_PMT*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Fitf_PMT_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Fitf_PMT*)nullptr)->GetClass();
      Fitf_PMT_TClassManip(theClass);
   return theClass;
   }

   static void Fitf_PMT_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Fitf_PMT_pedestal_Dictionary();
   static void Fitf_PMT_pedestal_TClassManip(TClass*);
   static void *new_Fitf_PMT_pedestal(void *p = nullptr);
   static void *newArray_Fitf_PMT_pedestal(Long_t size, void *p);
   static void delete_Fitf_PMT_pedestal(void *p);
   static void deleteArray_Fitf_PMT_pedestal(void *p);
   static void destruct_Fitf_PMT_pedestal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitf_PMT_pedestal*)
   {
      ::Fitf_PMT_pedestal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Fitf_PMT_pedestal));
      static ::ROOT::TGenericClassInfo 
         instance("Fitf_PMT_pedestal", "ReadRun.h", 521,
                  typeid(::Fitf_PMT_pedestal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Fitf_PMT_pedestal_Dictionary, isa_proxy, 4,
                  sizeof(::Fitf_PMT_pedestal) );
      instance.SetNew(&new_Fitf_PMT_pedestal);
      instance.SetNewArray(&newArray_Fitf_PMT_pedestal);
      instance.SetDelete(&delete_Fitf_PMT_pedestal);
      instance.SetDeleteArray(&deleteArray_Fitf_PMT_pedestal);
      instance.SetDestructor(&destruct_Fitf_PMT_pedestal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitf_PMT_pedestal*)
   {
      return GenerateInitInstanceLocal((::Fitf_PMT_pedestal*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitf_PMT_pedestal*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Fitf_PMT_pedestal_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Fitf_PMT_pedestal*)nullptr)->GetClass();
      Fitf_PMT_pedestal_TClassManip(theClass);
   return theClass;
   }

   static void Fitf_PMT_pedestal_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Fitf_PMT_ideal_Dictionary();
   static void Fitf_PMT_ideal_TClassManip(TClass*);
   static void *new_Fitf_PMT_ideal(void *p = nullptr);
   static void *newArray_Fitf_PMT_ideal(Long_t size, void *p);
   static void delete_Fitf_PMT_ideal(void *p);
   static void deleteArray_Fitf_PMT_ideal(void *p);
   static void destruct_Fitf_PMT_ideal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitf_PMT_ideal*)
   {
      ::Fitf_PMT_ideal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Fitf_PMT_ideal));
      static ::ROOT::TGenericClassInfo 
         instance("Fitf_PMT_ideal", "ReadRun.h", 572,
                  typeid(::Fitf_PMT_ideal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Fitf_PMT_ideal_Dictionary, isa_proxy, 4,
                  sizeof(::Fitf_PMT_ideal) );
      instance.SetNew(&new_Fitf_PMT_ideal);
      instance.SetNewArray(&newArray_Fitf_PMT_ideal);
      instance.SetDelete(&delete_Fitf_PMT_ideal);
      instance.SetDeleteArray(&deleteArray_Fitf_PMT_ideal);
      instance.SetDestructor(&destruct_Fitf_PMT_ideal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitf_PMT_ideal*)
   {
      return GenerateInitInstanceLocal((::Fitf_PMT_ideal*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitf_PMT_ideal*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Fitf_PMT_ideal_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Fitf_PMT_ideal*)nullptr)->GetClass();
      Fitf_PMT_ideal_TClassManip(theClass);
   return theClass;
   }

   static void Fitf_PMT_ideal_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Fitf_langaus_Dictionary();
   static void Fitf_langaus_TClassManip(TClass*);
   static void *new_Fitf_langaus(void *p = nullptr);
   static void *newArray_Fitf_langaus(Long_t size, void *p);
   static void delete_Fitf_langaus(void *p);
   static void deleteArray_Fitf_langaus(void *p);
   static void destruct_Fitf_langaus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitf_langaus*)
   {
      ::Fitf_langaus *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Fitf_langaus));
      static ::ROOT::TGenericClassInfo 
         instance("Fitf_langaus", "ReadRun.h", 609,
                  typeid(::Fitf_langaus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Fitf_langaus_Dictionary, isa_proxy, 4,
                  sizeof(::Fitf_langaus) );
      instance.SetNew(&new_Fitf_langaus);
      instance.SetNewArray(&newArray_Fitf_langaus);
      instance.SetDelete(&delete_Fitf_langaus);
      instance.SetDeleteArray(&deleteArray_Fitf_langaus);
      instance.SetDestructor(&destruct_Fitf_langaus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitf_langaus*)
   {
      return GenerateInitInstanceLocal((::Fitf_langaus*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitf_langaus*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Fitf_langaus_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Fitf_langaus*)nullptr)->GetClass();
      Fitf_langaus_TClassManip(theClass);
   return theClass;
   }

   static void Fitf_langaus_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ReadRun::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ReadRun::Class_Name()
{
   return "ReadRun";
}

//______________________________________________________________________________
const char *ReadRun::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ReadRun*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ReadRun::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ReadRun*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ReadRun::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ReadRun*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ReadRun::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ReadRun*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ReadRun::Streamer(TBuffer &R__b)
{
   // Stream an object of class ReadRun.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ReadRun::Class(),this);
   } else {
      R__b.WriteClassBuffer(ReadRun::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ReadRun(void *p) {
      return  p ? new(p) ::ReadRun : new ::ReadRun;
   }
   static void *newArray_ReadRun(Long_t nElements, void *p) {
      return p ? new(p) ::ReadRun[nElements] : new ::ReadRun[nElements];
   }
   // Wrapper around operator delete
   static void delete_ReadRun(void *p) {
      delete ((::ReadRun*)p);
   }
   static void deleteArray_ReadRun(void *p) {
      delete [] ((::ReadRun*)p);
   }
   static void destruct_ReadRun(void *p) {
      typedef ::ReadRun current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ReadRun

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitf(void *p) {
      return  p ? new(p) ::Fitf : new ::Fitf;
   }
   static void *newArray_Fitf(Long_t nElements, void *p) {
      return p ? new(p) ::Fitf[nElements] : new ::Fitf[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitf(void *p) {
      delete ((::Fitf*)p);
   }
   static void deleteArray_Fitf(void *p) {
      delete [] ((::Fitf*)p);
   }
   static void destruct_Fitf(void *p) {
      typedef ::Fitf current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Fitf

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitf_biased(void *p) {
      return  p ? new(p) ::Fitf_biased : new ::Fitf_biased;
   }
   static void *newArray_Fitf_biased(Long_t nElements, void *p) {
      return p ? new(p) ::Fitf_biased[nElements] : new ::Fitf_biased[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitf_biased(void *p) {
      delete ((::Fitf_biased*)p);
   }
   static void deleteArray_Fitf_biased(void *p) {
      delete [] ((::Fitf_biased*)p);
   }
   static void destruct_Fitf_biased(void *p) {
      typedef ::Fitf_biased current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Fitf_biased

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitf_PMT(void *p) {
      return  p ? new(p) ::Fitf_PMT : new ::Fitf_PMT;
   }
   static void *newArray_Fitf_PMT(Long_t nElements, void *p) {
      return p ? new(p) ::Fitf_PMT[nElements] : new ::Fitf_PMT[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitf_PMT(void *p) {
      delete ((::Fitf_PMT*)p);
   }
   static void deleteArray_Fitf_PMT(void *p) {
      delete [] ((::Fitf_PMT*)p);
   }
   static void destruct_Fitf_PMT(void *p) {
      typedef ::Fitf_PMT current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Fitf_PMT

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitf_PMT_pedestal(void *p) {
      return  p ? new(p) ::Fitf_PMT_pedestal : new ::Fitf_PMT_pedestal;
   }
   static void *newArray_Fitf_PMT_pedestal(Long_t nElements, void *p) {
      return p ? new(p) ::Fitf_PMT_pedestal[nElements] : new ::Fitf_PMT_pedestal[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitf_PMT_pedestal(void *p) {
      delete ((::Fitf_PMT_pedestal*)p);
   }
   static void deleteArray_Fitf_PMT_pedestal(void *p) {
      delete [] ((::Fitf_PMT_pedestal*)p);
   }
   static void destruct_Fitf_PMT_pedestal(void *p) {
      typedef ::Fitf_PMT_pedestal current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Fitf_PMT_pedestal

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitf_PMT_ideal(void *p) {
      return  p ? new(p) ::Fitf_PMT_ideal : new ::Fitf_PMT_ideal;
   }
   static void *newArray_Fitf_PMT_ideal(Long_t nElements, void *p) {
      return p ? new(p) ::Fitf_PMT_ideal[nElements] : new ::Fitf_PMT_ideal[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitf_PMT_ideal(void *p) {
      delete ((::Fitf_PMT_ideal*)p);
   }
   static void deleteArray_Fitf_PMT_ideal(void *p) {
      delete [] ((::Fitf_PMT_ideal*)p);
   }
   static void destruct_Fitf_PMT_ideal(void *p) {
      typedef ::Fitf_PMT_ideal current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Fitf_PMT_ideal

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitf_langaus(void *p) {
      return  p ? new(p) ::Fitf_langaus : new ::Fitf_langaus;
   }
   static void *newArray_Fitf_langaus(Long_t nElements, void *p) {
      return p ? new(p) ::Fitf_langaus[nElements] : new ::Fitf_langaus[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitf_langaus(void *p) {
      delete ((::Fitf_langaus*)p);
   }
   static void deleteArray_Fitf_langaus(void *p) {
      delete [] ((::Fitf_langaus*)p);
   }
   static void destruct_Fitf_langaus(void *p) {
      typedef ::Fitf_langaus current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Fitf_langaus

namespace ROOT {
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary();
   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p);
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<float> >*)
   {
      vector<vector<float> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<float> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<float> >", -2, "vector", 386,
                  typeid(vector<vector<float> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEfloatgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<float> >) );
      instance.SetNew(&new_vectorlEvectorlEfloatgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEfloatgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEfloatgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<float> > >()));

      ::ROOT::AddClassAlternate("vector<vector<float> >","std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<float> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<float> >*)nullptr)->GetClass();
      vectorlEvectorlEfloatgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> > : new vector<vector<float> >;
   }
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> >[nElements] : new vector<vector<float> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete ((vector<vector<float> >*)p);
   }
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete [] ((vector<vector<float> >*)p);
   }
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p) {
      typedef vector<vector<float> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<float> >

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "vector", 386,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));

      ::ROOT::AddClassAlternate("vector<unsigned int>","std::vector<unsigned int, std::allocator<unsigned int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned int>*)nullptr)->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete ((vector<unsigned int>*)p);
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] ((vector<unsigned int>*)p);
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 386,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)nullptr)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = nullptr);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 386,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)nullptr)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEboolgR_Dictionary();
   static void vectorlEboolgR_TClassManip(TClass*);
   static void *new_vectorlEboolgR(void *p = nullptr);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 592,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));

      ::ROOT::AddClassAlternate("vector<bool>","std::vector<bool, std::allocator<bool> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<bool>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<bool>*)nullptr)->GetClass();
      vectorlEboolgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEboolgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEboolgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool> : new vector<bool>;
   }
   static void *newArray_vectorlEboolgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool>[nElements] : new vector<bool>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEboolgR(void *p) {
      delete ((vector<bool>*)p);
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] ((vector<bool>*)p);
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace ROOT {
   static TClass *vectorlETFitResultPtrgR_Dictionary();
   static void vectorlETFitResultPtrgR_TClassManip(TClass*);
   static void *new_vectorlETFitResultPtrgR(void *p = nullptr);
   static void *newArray_vectorlETFitResultPtrgR(Long_t size, void *p);
   static void delete_vectorlETFitResultPtrgR(void *p);
   static void deleteArray_vectorlETFitResultPtrgR(void *p);
   static void destruct_vectorlETFitResultPtrgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TFitResultPtr>*)
   {
      vector<TFitResultPtr> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TFitResultPtr>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TFitResultPtr>", -2, "vector", 386,
                  typeid(vector<TFitResultPtr>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETFitResultPtrgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TFitResultPtr>) );
      instance.SetNew(&new_vectorlETFitResultPtrgR);
      instance.SetNewArray(&newArray_vectorlETFitResultPtrgR);
      instance.SetDelete(&delete_vectorlETFitResultPtrgR);
      instance.SetDeleteArray(&deleteArray_vectorlETFitResultPtrgR);
      instance.SetDestructor(&destruct_vectorlETFitResultPtrgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TFitResultPtr> >()));

      ::ROOT::AddClassAlternate("vector<TFitResultPtr>","std::vector<TFitResultPtr, std::allocator<TFitResultPtr> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TFitResultPtr>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETFitResultPtrgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TFitResultPtr>*)nullptr)->GetClass();
      vectorlETFitResultPtrgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETFitResultPtrgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETFitResultPtrgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TFitResultPtr> : new vector<TFitResultPtr>;
   }
   static void *newArray_vectorlETFitResultPtrgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TFitResultPtr>[nElements] : new vector<TFitResultPtr>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETFitResultPtrgR(void *p) {
      delete ((vector<TFitResultPtr>*)p);
   }
   static void deleteArray_vectorlETFitResultPtrgR(void *p) {
      delete [] ((vector<TFitResultPtr>*)p);
   }
   static void destruct_vectorlETFitResultPtrgR(void *p) {
      typedef vector<TFitResultPtr> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TFitResultPtr>

namespace {
  void TriggerDictionaryInitialization_ReadRunDictUX_Impl() {
    static const char* headers[] = {
"ReadRun.h",
nullptr
    };
    static const char* includePaths[] = {
"/opt/root/include/",
"/mnt/d/Work_SHK_Bachelor/analysis_programm/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "ReadRunDictUX dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  ReadRun;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  Fitf;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  Fitf_biased;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  Fitf_PMT;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  Fitf_PMT_pedestal;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  Fitf_PMT_ideal;
class __attribute__((annotate("$clingAutoload$ReadRun.h")))  Fitf_langaus;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "ReadRunDictUX dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "ReadRun.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Fitf", payloadCode, "@",
"Fitf_PMT", payloadCode, "@",
"Fitf_PMT_ideal", payloadCode, "@",
"Fitf_PMT_pedestal", payloadCode, "@",
"Fitf_biased", payloadCode, "@",
"Fitf_langaus", payloadCode, "@",
"ReadRun", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ReadRunDictUX",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ReadRunDictUX_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ReadRunDictUX_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ReadRunDictUX() {
  TriggerDictionaryInitialization_ReadRunDictUX_Impl();
}
