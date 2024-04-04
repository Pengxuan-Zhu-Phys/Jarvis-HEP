/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/FastJetFinder.h"
#include "modules/FastJetGridMedianEstimator.h"
#include "modules/RunPUPPI.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class FastJetFinder+;
#pragma link C++ class FastJetGridMedianEstimator+;
#pragma link C++ class RunPUPPI+;

#endif
// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdImodulesdIFastJetDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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
#include "modules/FastJetFinder.h"
#include "modules/FastJetGridMedianEstimator.h"
#include "modules/RunPUPPI.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_FastJetFinder(void *p = nullptr);
   static void *newArray_FastJetFinder(Long_t size, void *p);
   static void delete_FastJetFinder(void *p);
   static void deleteArray_FastJetFinder(void *p);
   static void destruct_FastJetFinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastJetFinder*)
   {
      ::FastJetFinder *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastJetFinder >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FastJetFinder", ::FastJetFinder::Class_Version(), "modules/FastJetFinder.h", 51,
                  typeid(::FastJetFinder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastJetFinder::Dictionary, isa_proxy, 4,
                  sizeof(::FastJetFinder) );
      instance.SetNew(&new_FastJetFinder);
      instance.SetNewArray(&newArray_FastJetFinder);
      instance.SetDelete(&delete_FastJetFinder);
      instance.SetDeleteArray(&deleteArray_FastJetFinder);
      instance.SetDestructor(&destruct_FastJetFinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastJetFinder*)
   {
      return GenerateInitInstanceLocal(static_cast<::FastJetFinder*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FastJetFinder*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FastJetGridMedianEstimator(void *p = nullptr);
   static void *newArray_FastJetGridMedianEstimator(Long_t size, void *p);
   static void delete_FastJetGridMedianEstimator(void *p);
   static void deleteArray_FastJetGridMedianEstimator(void *p);
   static void destruct_FastJetGridMedianEstimator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastJetGridMedianEstimator*)
   {
      ::FastJetGridMedianEstimator *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastJetGridMedianEstimator >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FastJetGridMedianEstimator", ::FastJetGridMedianEstimator::Class_Version(), "modules/FastJetGridMedianEstimator.h", 41,
                  typeid(::FastJetGridMedianEstimator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastJetGridMedianEstimator::Dictionary, isa_proxy, 4,
                  sizeof(::FastJetGridMedianEstimator) );
      instance.SetNew(&new_FastJetGridMedianEstimator);
      instance.SetNewArray(&newArray_FastJetGridMedianEstimator);
      instance.SetDelete(&delete_FastJetGridMedianEstimator);
      instance.SetDeleteArray(&deleteArray_FastJetGridMedianEstimator);
      instance.SetDestructor(&destruct_FastJetGridMedianEstimator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastJetGridMedianEstimator*)
   {
      return GenerateInitInstanceLocal(static_cast<::FastJetGridMedianEstimator*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FastJetGridMedianEstimator*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RunPUPPI(void *p = nullptr);
   static void *newArray_RunPUPPI(Long_t size, void *p);
   static void delete_RunPUPPI(void *p);
   static void deleteArray_RunPUPPI(void *p);
   static void destruct_RunPUPPI(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RunPUPPI*)
   {
      ::RunPUPPI *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RunPUPPI >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("RunPUPPI", ::RunPUPPI::Class_Version(), "modules/RunPUPPI.h", 11,
                  typeid(::RunPUPPI), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RunPUPPI::Dictionary, isa_proxy, 4,
                  sizeof(::RunPUPPI) );
      instance.SetNew(&new_RunPUPPI);
      instance.SetNewArray(&newArray_RunPUPPI);
      instance.SetDelete(&delete_RunPUPPI);
      instance.SetDeleteArray(&deleteArray_RunPUPPI);
      instance.SetDestructor(&destruct_RunPUPPI);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RunPUPPI*)
   {
      return GenerateInitInstanceLocal(static_cast<::RunPUPPI*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::RunPUPPI*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr FastJetFinder::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FastJetFinder::Class_Name()
{
   return "FastJetFinder";
}

//______________________________________________________________________________
const char *FastJetFinder::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastJetFinder*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FastJetFinder::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastJetFinder*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastJetFinder::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastJetFinder*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastJetFinder::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastJetFinder*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastJetGridMedianEstimator::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FastJetGridMedianEstimator::Class_Name()
{
   return "FastJetGridMedianEstimator";
}

//______________________________________________________________________________
const char *FastJetGridMedianEstimator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastJetGridMedianEstimator*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FastJetGridMedianEstimator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastJetGridMedianEstimator*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastJetGridMedianEstimator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastJetGridMedianEstimator*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastJetGridMedianEstimator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastJetGridMedianEstimator*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RunPUPPI::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *RunPUPPI::Class_Name()
{
   return "RunPUPPI";
}

//______________________________________________________________________________
const char *RunPUPPI::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RunPUPPI*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int RunPUPPI::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RunPUPPI*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RunPUPPI::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RunPUPPI*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RunPUPPI::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RunPUPPI*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void FastJetFinder::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastJetFinder.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FastJetFinder::Class(),this);
   } else {
      R__b.WriteClassBuffer(FastJetFinder::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastJetFinder(void *p) {
      return  p ? new(p) ::FastJetFinder : new ::FastJetFinder;
   }
   static void *newArray_FastJetFinder(Long_t nElements, void *p) {
      return p ? new(p) ::FastJetFinder[nElements] : new ::FastJetFinder[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastJetFinder(void *p) {
      delete (static_cast<::FastJetFinder*>(p));
   }
   static void deleteArray_FastJetFinder(void *p) {
      delete [] (static_cast<::FastJetFinder*>(p));
   }
   static void destruct_FastJetFinder(void *p) {
      typedef ::FastJetFinder current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FastJetFinder

//______________________________________________________________________________
void FastJetGridMedianEstimator::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastJetGridMedianEstimator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FastJetGridMedianEstimator::Class(),this);
   } else {
      R__b.WriteClassBuffer(FastJetGridMedianEstimator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastJetGridMedianEstimator(void *p) {
      return  p ? new(p) ::FastJetGridMedianEstimator : new ::FastJetGridMedianEstimator;
   }
   static void *newArray_FastJetGridMedianEstimator(Long_t nElements, void *p) {
      return p ? new(p) ::FastJetGridMedianEstimator[nElements] : new ::FastJetGridMedianEstimator[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastJetGridMedianEstimator(void *p) {
      delete (static_cast<::FastJetGridMedianEstimator*>(p));
   }
   static void deleteArray_FastJetGridMedianEstimator(void *p) {
      delete [] (static_cast<::FastJetGridMedianEstimator*>(p));
   }
   static void destruct_FastJetGridMedianEstimator(void *p) {
      typedef ::FastJetGridMedianEstimator current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FastJetGridMedianEstimator

//______________________________________________________________________________
void RunPUPPI::Streamer(TBuffer &R__b)
{
   // Stream an object of class RunPUPPI.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RunPUPPI::Class(),this);
   } else {
      R__b.WriteClassBuffer(RunPUPPI::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RunPUPPI(void *p) {
      return  p ? new(p) ::RunPUPPI : new ::RunPUPPI;
   }
   static void *newArray_RunPUPPI(Long_t nElements, void *p) {
      return p ? new(p) ::RunPUPPI[nElements] : new ::RunPUPPI[nElements];
   }
   // Wrapper around operator delete
   static void delete_RunPUPPI(void *p) {
      delete (static_cast<::RunPUPPI*>(p));
   }
   static void deleteArray_RunPUPPI(void *p) {
      delete [] (static_cast<::RunPUPPI*>(p));
   }
   static void destruct_RunPUPPI(void *p) {
      typedef ::RunPUPPI current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::RunPUPPI

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
         instance("vector<int>", -2, "vector", 389,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr))->GetClass();
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
      delete (static_cast<vector<int>*>(p));
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] (static_cast<vector<int>*>(p));
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      (static_cast<current_t*>(p))->~current_t();
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
         instance("vector<float>", -2, "vector", 389,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<float>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<float>*>(nullptr))->GetClass();
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
      delete (static_cast<vector<float>*>(p));
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] (static_cast<vector<float>*>(p));
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      (static_cast<current_t*>(p))->~current_t();
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
         instance("vector<bool>", -2, "vector", 596,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<bool>","std::vector<bool, std::allocator<bool> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<bool>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<bool>*>(nullptr))->GetClass();
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
      delete (static_cast<vector<bool>*>(p));
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] (static_cast<vector<bool>*>(p));
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace {
  void TriggerDictionaryInitialization_FastJetDict_Impl() {
    static const char* headers[] = {
nullptr
    };
    static const char* includePaths[] = {
"external",
"/home/buding/root/include/",
"/home/buding/Jarvis-HEP/External/Library/Delphes/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "FastJetDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$modules/FastJetFinder.h")))  FastJetFinder;
class __attribute__((annotate("$clingAutoload$modules/FastJetGridMedianEstimator.h")))  FastJetGridMedianEstimator;
class __attribute__((annotate("$clingAutoload$modules/RunPUPPI.h")))  RunPUPPI;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "FastJetDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/FastJetFinder.h"
#include "modules/FastJetGridMedianEstimator.h"
#include "modules/RunPUPPI.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class FastJetFinder+;
#pragma link C++ class FastJetGridMedianEstimator+;
#pragma link C++ class RunPUPPI+;

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"FastJetFinder", payloadCode, "@",
"FastJetGridMedianEstimator", payloadCode, "@",
"RunPUPPI", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("FastJetDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_FastJetDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_FastJetDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_FastJetDict() {
  TriggerDictionaryInitialization_FastJetDict_Impl();
}
