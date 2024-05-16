
/** \class ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTask.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ExRootTreeReader+;
#pragma link C++ class ExRootTreeBranch+;
#pragma link C++ class ExRootTreeWriter+;
#pragma link C++ class ExRootResult+;
#pragma link C++ class ExRootClassifier+;
#pragma link C++ class ExRootFilter+;

#pragma link C++ class ExRootProgressBar+;
#pragma link C++ class ExRootConfReader+;
#pragma link C++ class ExRootConfParam+;
#pragma link C++ class ExRootTask+;

#pragma link C++ function HistStyle;
#pragma link C++ function FillChain;

#endif

// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdIexternaldIExRootAnalysisdIExRootAnalysisDict
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
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTask.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_ExRootTreeReader(void *p = nullptr);
   static void *newArray_ExRootTreeReader(Long_t size, void *p);
   static void delete_ExRootTreeReader(void *p);
   static void deleteArray_ExRootTreeReader(void *p);
   static void destruct_ExRootTreeReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTreeReader*)
   {
      ::ExRootTreeReader *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootTreeReader >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTreeReader", ::ExRootTreeReader::Class_Version(), "ExRootAnalysis/ExRootTreeReader.h", 19,
                  typeid(::ExRootTreeReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ExRootTreeReader::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTreeReader) );
      instance.SetNew(&new_ExRootTreeReader);
      instance.SetNewArray(&newArray_ExRootTreeReader);
      instance.SetDelete(&delete_ExRootTreeReader);
      instance.SetDeleteArray(&deleteArray_ExRootTreeReader);
      instance.SetDestructor(&destruct_ExRootTreeReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTreeReader*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootTreeReader*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootTreeReader*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ExRootTreeWriter(void *p = nullptr);
   static void *newArray_ExRootTreeWriter(Long_t size, void *p);
   static void delete_ExRootTreeWriter(void *p);
   static void deleteArray_ExRootTreeWriter(void *p);
   static void destruct_ExRootTreeWriter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTreeWriter*)
   {
      ::ExRootTreeWriter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootTreeWriter >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTreeWriter", ::ExRootTreeWriter::Class_Version(), "ExRootAnalysis/ExRootTreeWriter.h", 21,
                  typeid(::ExRootTreeWriter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ExRootTreeWriter::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTreeWriter) );
      instance.SetNew(&new_ExRootTreeWriter);
      instance.SetNewArray(&newArray_ExRootTreeWriter);
      instance.SetDelete(&delete_ExRootTreeWriter);
      instance.SetDeleteArray(&deleteArray_ExRootTreeWriter);
      instance.SetDestructor(&destruct_ExRootTreeWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTreeWriter*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootTreeWriter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootTreeWriter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootTreeBranch_Dictionary();
   static void ExRootTreeBranch_TClassManip(TClass*);
   static void delete_ExRootTreeBranch(void *p);
   static void deleteArray_ExRootTreeBranch(void *p);
   static void destruct_ExRootTreeBranch(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTreeBranch*)
   {
      ::ExRootTreeBranch *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootTreeBranch));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTreeBranch", "ExRootAnalysis/ExRootTreeBranch.h", 18,
                  typeid(::ExRootTreeBranch), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExRootTreeBranch_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTreeBranch) );
      instance.SetDelete(&delete_ExRootTreeBranch);
      instance.SetDeleteArray(&deleteArray_ExRootTreeBranch);
      instance.SetDestructor(&destruct_ExRootTreeBranch);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTreeBranch*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootTreeBranch*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootTreeBranch*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootTreeBranch_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ExRootTreeBranch*>(nullptr))->GetClass();
      ExRootTreeBranch_TClassManip(theClass);
   return theClass;
   }

   static void ExRootTreeBranch_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootResult_Dictionary();
   static void ExRootResult_TClassManip(TClass*);
   static void *new_ExRootResult(void *p = nullptr);
   static void *newArray_ExRootResult(Long_t size, void *p);
   static void delete_ExRootResult(void *p);
   static void deleteArray_ExRootResult(void *p);
   static void destruct_ExRootResult(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootResult*)
   {
      ::ExRootResult *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootResult));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootResult", "ExRootAnalysis/ExRootResult.h", 21,
                  typeid(::ExRootResult), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExRootResult_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootResult) );
      instance.SetNew(&new_ExRootResult);
      instance.SetNewArray(&newArray_ExRootResult);
      instance.SetDelete(&delete_ExRootResult);
      instance.SetDeleteArray(&deleteArray_ExRootResult);
      instance.SetDestructor(&destruct_ExRootResult);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootResult*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootResult*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootResult*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootResult_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ExRootResult*>(nullptr))->GetClass();
      ExRootResult_TClassManip(theClass);
   return theClass;
   }

   static void ExRootResult_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootClassifier_Dictionary();
   static void ExRootClassifier_TClassManip(TClass*);
   static void delete_ExRootClassifier(void *p);
   static void deleteArray_ExRootClassifier(void *p);
   static void destruct_ExRootClassifier(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootClassifier*)
   {
      ::ExRootClassifier *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootClassifier));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootClassifier", "ExRootAnalysis/ExRootClassifier.h", 8,
                  typeid(::ExRootClassifier), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExRootClassifier_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootClassifier) );
      instance.SetDelete(&delete_ExRootClassifier);
      instance.SetDeleteArray(&deleteArray_ExRootClassifier);
      instance.SetDestructor(&destruct_ExRootClassifier);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootClassifier*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootClassifier*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootClassifier*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootClassifier_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ExRootClassifier*>(nullptr))->GetClass();
      ExRootClassifier_TClassManip(theClass);
   return theClass;
   }

   static void ExRootClassifier_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootFilter_Dictionary();
   static void ExRootFilter_TClassManip(TClass*);
   static void delete_ExRootFilter(void *p);
   static void deleteArray_ExRootFilter(void *p);
   static void destruct_ExRootFilter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootFilter*)
   {
      ::ExRootFilter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootFilter));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootFilter", "ExRootAnalysis/ExRootFilter.h", 13,
                  typeid(::ExRootFilter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExRootFilter_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootFilter) );
      instance.SetDelete(&delete_ExRootFilter);
      instance.SetDeleteArray(&deleteArray_ExRootFilter);
      instance.SetDestructor(&destruct_ExRootFilter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootFilter*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootFilter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootFilter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootFilter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ExRootFilter*>(nullptr))->GetClass();
      ExRootFilter_TClassManip(theClass);
   return theClass;
   }

   static void ExRootFilter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootProgressBar_Dictionary();
   static void ExRootProgressBar_TClassManip(TClass*);
   static void delete_ExRootProgressBar(void *p);
   static void deleteArray_ExRootProgressBar(void *p);
   static void destruct_ExRootProgressBar(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootProgressBar*)
   {
      ::ExRootProgressBar *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootProgressBar));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootProgressBar", "ExRootAnalysis/ExRootProgressBar.h", 6,
                  typeid(::ExRootProgressBar), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExRootProgressBar_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootProgressBar) );
      instance.SetDelete(&delete_ExRootProgressBar);
      instance.SetDeleteArray(&deleteArray_ExRootProgressBar);
      instance.SetDestructor(&destruct_ExRootProgressBar);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootProgressBar*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootProgressBar*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootProgressBar*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootProgressBar_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ExRootProgressBar*>(nullptr))->GetClass();
      ExRootProgressBar_TClassManip(theClass);
   return theClass;
   }

   static void ExRootProgressBar_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootConfParam_Dictionary();
   static void ExRootConfParam_TClassManip(TClass*);
   static void *new_ExRootConfParam(void *p = nullptr);
   static void *newArray_ExRootConfParam(Long_t size, void *p);
   static void delete_ExRootConfParam(void *p);
   static void deleteArray_ExRootConfParam(void *p);
   static void destruct_ExRootConfParam(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootConfParam*)
   {
      ::ExRootConfParam *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootConfParam));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootConfParam", "ExRootAnalysis/ExRootConfReader.h", 20,
                  typeid(::ExRootConfParam), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExRootConfParam_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootConfParam) );
      instance.SetNew(&new_ExRootConfParam);
      instance.SetNewArray(&newArray_ExRootConfParam);
      instance.SetDelete(&delete_ExRootConfParam);
      instance.SetDeleteArray(&deleteArray_ExRootConfParam);
      instance.SetDestructor(&destruct_ExRootConfParam);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootConfParam*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootConfParam*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootConfParam*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootConfParam_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ExRootConfParam*>(nullptr))->GetClass();
      ExRootConfParam_TClassManip(theClass);
   return theClass;
   }

   static void ExRootConfParam_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_ExRootConfReader(void *p = nullptr);
   static void *newArray_ExRootConfReader(Long_t size, void *p);
   static void delete_ExRootConfReader(void *p);
   static void deleteArray_ExRootConfReader(void *p);
   static void destruct_ExRootConfReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootConfReader*)
   {
      ::ExRootConfReader *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootConfReader >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootConfReader", ::ExRootConfReader::Class_Version(), "ExRootAnalysis/ExRootConfReader.h", 42,
                  typeid(::ExRootConfReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ExRootConfReader::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootConfReader) );
      instance.SetNew(&new_ExRootConfReader);
      instance.SetNewArray(&newArray_ExRootConfReader);
      instance.SetDelete(&delete_ExRootConfReader);
      instance.SetDeleteArray(&deleteArray_ExRootConfReader);
      instance.SetDestructor(&destruct_ExRootConfReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootConfReader*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootConfReader*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootConfReader*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ExRootTask(void *p = nullptr);
   static void *newArray_ExRootTask(Long_t size, void *p);
   static void delete_ExRootTask(void *p);
   static void deleteArray_ExRootTask(void *p);
   static void destruct_ExRootTask(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTask*)
   {
      ::ExRootTask *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootTask >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTask", ::ExRootTask::Class_Version(), "ExRootAnalysis/ExRootTask.h", 19,
                  typeid(::ExRootTask), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ExRootTask::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTask) );
      instance.SetNew(&new_ExRootTask);
      instance.SetNewArray(&newArray_ExRootTask);
      instance.SetDelete(&delete_ExRootTask);
      instance.SetDeleteArray(&deleteArray_ExRootTask);
      instance.SetDestructor(&destruct_ExRootTask);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTask*)
   {
      return GenerateInitInstanceLocal(static_cast<::ExRootTask*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ExRootTask*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ExRootTreeReader::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootTreeReader::Class_Name()
{
   return "ExRootTreeReader";
}

//______________________________________________________________________________
const char *ExRootTreeReader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootTreeReader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootTreeReader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootTreeReader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExRootTreeWriter::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootTreeWriter::Class_Name()
{
   return "ExRootTreeWriter";
}

//______________________________________________________________________________
const char *ExRootTreeWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootTreeWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootTreeWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootTreeWriter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExRootConfReader::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootConfReader::Class_Name()
{
   return "ExRootConfReader";
}

//______________________________________________________________________________
const char *ExRootConfReader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootConfReader*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootConfReader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootConfReader*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootConfReader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootConfReader*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootConfReader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootConfReader*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExRootTask::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootTask::Class_Name()
{
   return "ExRootTask";
}

//______________________________________________________________________________
const char *ExRootTask::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTask*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootTask::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTask*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootTask::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTask*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootTask::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTask*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ExRootTreeReader::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootTreeReader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootTreeReader::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootTreeReader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootTreeReader(void *p) {
      return  p ? new(p) ::ExRootTreeReader : new ::ExRootTreeReader;
   }
   static void *newArray_ExRootTreeReader(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootTreeReader[nElements] : new ::ExRootTreeReader[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootTreeReader(void *p) {
      delete (static_cast<::ExRootTreeReader*>(p));
   }
   static void deleteArray_ExRootTreeReader(void *p) {
      delete [] (static_cast<::ExRootTreeReader*>(p));
   }
   static void destruct_ExRootTreeReader(void *p) {
      typedef ::ExRootTreeReader current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTreeReader

//______________________________________________________________________________
void ExRootTreeWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootTreeWriter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootTreeWriter::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootTreeWriter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootTreeWriter(void *p) {
      return  p ? new(p) ::ExRootTreeWriter : new ::ExRootTreeWriter;
   }
   static void *newArray_ExRootTreeWriter(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootTreeWriter[nElements] : new ::ExRootTreeWriter[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootTreeWriter(void *p) {
      delete (static_cast<::ExRootTreeWriter*>(p));
   }
   static void deleteArray_ExRootTreeWriter(void *p) {
      delete [] (static_cast<::ExRootTreeWriter*>(p));
   }
   static void destruct_ExRootTreeWriter(void *p) {
      typedef ::ExRootTreeWriter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTreeWriter

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootTreeBranch(void *p) {
      delete (static_cast<::ExRootTreeBranch*>(p));
   }
   static void deleteArray_ExRootTreeBranch(void *p) {
      delete [] (static_cast<::ExRootTreeBranch*>(p));
   }
   static void destruct_ExRootTreeBranch(void *p) {
      typedef ::ExRootTreeBranch current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTreeBranch

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootResult(void *p) {
      return  p ? new(p) ::ExRootResult : new ::ExRootResult;
   }
   static void *newArray_ExRootResult(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootResult[nElements] : new ::ExRootResult[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootResult(void *p) {
      delete (static_cast<::ExRootResult*>(p));
   }
   static void deleteArray_ExRootResult(void *p) {
      delete [] (static_cast<::ExRootResult*>(p));
   }
   static void destruct_ExRootResult(void *p) {
      typedef ::ExRootResult current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootResult

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootClassifier(void *p) {
      delete (static_cast<::ExRootClassifier*>(p));
   }
   static void deleteArray_ExRootClassifier(void *p) {
      delete [] (static_cast<::ExRootClassifier*>(p));
   }
   static void destruct_ExRootClassifier(void *p) {
      typedef ::ExRootClassifier current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootClassifier

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootFilter(void *p) {
      delete (static_cast<::ExRootFilter*>(p));
   }
   static void deleteArray_ExRootFilter(void *p) {
      delete [] (static_cast<::ExRootFilter*>(p));
   }
   static void destruct_ExRootFilter(void *p) {
      typedef ::ExRootFilter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootFilter

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootProgressBar(void *p) {
      delete (static_cast<::ExRootProgressBar*>(p));
   }
   static void deleteArray_ExRootProgressBar(void *p) {
      delete [] (static_cast<::ExRootProgressBar*>(p));
   }
   static void destruct_ExRootProgressBar(void *p) {
      typedef ::ExRootProgressBar current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootProgressBar

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootConfParam(void *p) {
      return  p ? new(p) ::ExRootConfParam : new ::ExRootConfParam;
   }
   static void *newArray_ExRootConfParam(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootConfParam[nElements] : new ::ExRootConfParam[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootConfParam(void *p) {
      delete (static_cast<::ExRootConfParam*>(p));
   }
   static void deleteArray_ExRootConfParam(void *p) {
      delete [] (static_cast<::ExRootConfParam*>(p));
   }
   static void destruct_ExRootConfParam(void *p) {
      typedef ::ExRootConfParam current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootConfParam

//______________________________________________________________________________
void ExRootConfReader::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootConfReader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootConfReader::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootConfReader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootConfReader(void *p) {
      return  p ? new(p) ::ExRootConfReader : new ::ExRootConfReader;
   }
   static void *newArray_ExRootConfReader(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootConfReader[nElements] : new ::ExRootConfReader[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootConfReader(void *p) {
      delete (static_cast<::ExRootConfReader*>(p));
   }
   static void deleteArray_ExRootConfReader(void *p) {
      delete [] (static_cast<::ExRootConfReader*>(p));
   }
   static void destruct_ExRootConfReader(void *p) {
      typedef ::ExRootConfReader current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootConfReader

//______________________________________________________________________________
void ExRootTask::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootTask.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootTask::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootTask::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootTask(void *p) {
      return  p ? new(p) ::ExRootTask : new ::ExRootTask;
   }
   static void *newArray_ExRootTask(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootTask[nElements] : new ::ExRootTask[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootTask(void *p) {
      delete (static_cast<::ExRootTask*>(p));
   }
   static void deleteArray_ExRootTask(void *p) {
      delete [] (static_cast<::ExRootTask*>(p));
   }
   static void destruct_ExRootTask(void *p) {
      typedef ::ExRootTask current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTask

namespace {
  void TriggerDictionaryInitialization_ExRootAnalysisDict_Impl() {
    static const char* headers[] = {
nullptr
    };
    static const char* includePaths[] = {
"external",
"/opt/homebrew/Cellar/root/6.30.04/include/root",
"/Users/buding/Workshop/Jarvis/External/Library/Delphes/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "ExRootAnalysisDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootTreeReader.h")))  ExRootTreeReader;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootTreeWriter.h")))  ExRootTreeWriter;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootTreeBranch.h")))  ExRootTreeBranch;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootResult.h")))  ExRootResult;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootClassifier.h")))  ExRootClassifier;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootFilter.h")))  ExRootFilter;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootProgressBar.h")))  ExRootProgressBar;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootConfReader.h")))  ExRootConfParam;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootConfReader.h")))  ExRootConfReader;
class __attribute__((annotate("$clingAutoload$ExRootAnalysis/ExRootTask.h")))  ExRootTask;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "ExRootAnalysisDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers

/** \class ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTask.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ExRootTreeReader+;
#pragma link C++ class ExRootTreeBranch+;
#pragma link C++ class ExRootTreeWriter+;
#pragma link C++ class ExRootResult+;
#pragma link C++ class ExRootClassifier+;
#pragma link C++ class ExRootFilter+;

#pragma link C++ class ExRootProgressBar+;
#pragma link C++ class ExRootConfReader+;
#pragma link C++ class ExRootConfParam+;
#pragma link C++ class ExRootTask+;

#pragma link C++ function HistStyle;
#pragma link C++ function FillChain;

#endif


#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ExRootClassifier", payloadCode, "@",
"ExRootConfParam", payloadCode, "@",
"ExRootConfReader", payloadCode, "@",
"ExRootFilter", payloadCode, "@",
"ExRootProgressBar", payloadCode, "@",
"ExRootResult", payloadCode, "@",
"ExRootTask", payloadCode, "@",
"ExRootTreeBranch", payloadCode, "@",
"ExRootTreeReader", payloadCode, "@",
"ExRootTreeWriter", payloadCode, "@",
"FillChain", payloadCode, "@",
"HistStyle", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ExRootAnalysisDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ExRootAnalysisDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ExRootAnalysisDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ExRootAnalysisDict() {
  TriggerDictionaryInitialization_ExRootAnalysisDict_Impl();
}
