/*!
 * \file   QwTrackedRNTupleModel.h
 * \brief  Wrapper around ROOT::RNTupleModel that tracks field names
 *
 * This wrapper intercepts MakeField calls and records field names alongside
 * the field pointers. This enables per-detector RNTuple functionality by
 * allowing us to group fields by detector after model construction.
 */

#ifndef __QWTRACKEDRNTUPLEMODEL__
#define __QWTRACKEDRNTUPLEMODEL__

#ifdef HAS_RNTUPLE_SUPPORT

#include <ROOT/RNTupleModel.hxx>
#include <string>
#include <vector>
#include <memory>
#include "QwLog.h"

/**
 * \class QwTrackedRNTupleModel
 * \brief Wrapper around ROOT::RNTupleModel that tracks field names
 * 
 * This wrapper intercepts MakeField calls and records field names alongside
 * the field pointers. This enables per-detector RNTuple functionality by
 * allowing us to group fields by detector after model construction.
 */
class QwTrackedRNTupleModel {
public:
  QwTrackedRNTupleModel() {
    fModel = ROOT::RNTupleModel::Create();
  }
  
  /// Create a field and track its name
  template<typename T>
  std::shared_ptr<T> MakeField(const std::string& field_name) {
    auto field_ptr = fModel->MakeField<T>(field_name);
    fFieldNames.push_back(field_name);
    return field_ptr;
  }
  
  /// Overload for TString (commonly used in JAPAN)
  template<typename T>
  std::shared_ptr<T> MakeField(const TString& field_name) {
    return MakeField<T>(std::string(field_name.Data()));
  }
  
  /// Overload for const char* (commonly used)
  template<typename T>
  std::shared_ptr<T> MakeField(const char* field_name) {
    return MakeField<T>(std::string(field_name));
  }
  
  /// Get the underlying ROOT model
  std::unique_ptr<ROOT::RNTupleModel>& GetModel() { return fModel; }
  
  /// Get the default entry from the model
  ROOT::REntry& GetDefaultEntry() { return fModel->GetDefaultEntry(); }
  
  /// Release ownership of the model (for transferring to RNTupleWriter)
  std::unique_ptr<ROOT::RNTupleModel> ReleaseModel() { return std::move(fModel); }
  
  /// Get the tracked field names
  const std::vector<std::string>& GetFieldNames() const { return fFieldNames; }
  
private:
  std::unique_ptr<ROOT::RNTupleModel> fModel;
  std::vector<std::string> fFieldNames;
};

#endif // HAS_RNTUPLE_SUPPORT

#endif // __QWTRACKEDRNTUPLEMODEL__
