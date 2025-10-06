#ifndef __QWPERDETECTORRNTUPLES__
#define __QWPERDETECTORRNTUPLES__

#ifdef HAS_RNTUPLE_SUPPORT

// System headers
#include <map>
#include <string>
#include <memory>
#include <vector>

// ROOT headers
#include "TFile.h"
#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RField.hxx"
#include "ROOT/RNTupleWriter.hxx"

// Qweak headers
#include "QwLog.h"
#include "QwOptions.h"

/**
 * \class QwPerDetectorRNTuples
 * \brief Manages multiple RNTuples, one per detector/branch
 * 
 * Instead of having one large RNTuple with ~18,000 fields (650 branches x 28 fields each),
 * this class creates a separate RNTuple for each detector/branch, each containing only ~28 fields.
 * This should significantly improve write performance by reducing the row size.
 * 
 * Key benefits:
 * - Smaller row size per RNTuple (~28 fields vs ~18,000 fields)
 * - Better columnar compression (each detector's data is separate)
 * - More efficient I/O patterns
 * - Easier to parallelize writes in the future
 */
class QwPerDetectorRNTuples {

public:
  /// Constructor
  QwPerDetectorRNTuples(const std::string& base_name, const std::string& description);
  
  /// Destructor
  ~QwPerDetectorRNTuples();
  
  /// Initialize with TFile
  void Initialize(TFile* file);
  
  /// Register a detector with its fields
  void RegisterDetector(const std::string& detector_name, 
                       const std::vector<std::string>& field_names);
  
  /// Fill data for a specific detector
  void FillDetector(const std::string& detector_name, 
                   const std::vector<Double_t>& values);
  
  /// Commit all RNTuple clusters (called periodically for batching)
  void CommitAllClusters();
  
  /// Close all RNTuples and finalize
  void Close();
  
  /// Set cluster size for all RNTuples
  void SetClusterSize(UInt_t cluster_size);
  
  /// Enable/disable batching
  void SetBatching(Bool_t enable);
  
  /// Get number of registered detectors
  size_t GetNumDetectors() const { return fDetectorWriters.size(); }
  
  /// Check if a detector is registered
  bool HasDetector(const std::string& detector_name) const;
  
  /// Get statistics
  UInt_t GetEventsInCurrentCluster() const { return fEventsInCurrentCluster; }
  UInt_t GetClusterSize() const { return fClusterSize; }
  
  /// Extract detector name from field name (e.g., "bpm1c10WS_hw_sum" -> "bpm1c10WS")
  static std::string ExtractDetectorName(const std::string& field_name);
  
  /// Group fields by detector name
  static std::map<std::string, std::vector<std::string>> GroupFieldsByDetector(
    const std::vector<std::string>& field_names);
  
private:
  
  /// Structure to hold per-detector RNTuple data
  struct DetectorRNTuple {
    std::string name;
    std::unique_ptr<ROOT::RNTupleModel> model;
    std::unique_ptr<ROOT::RNTupleWriter> writer;
    std::vector<std::shared_ptr<Double_t>> field_ptrs;
    std::vector<std::string> field_names;
  };
  
  /// Base name for RNTuples (e.g., "evts" or "muls")
  std::string fBaseName;
  
  /// Description
  std::string fDescription;
  
  /// TFile pointer
  TFile* fFile;
  
  /// Map of detector name to DetectorRNTuple structure
  std::map<std::string, DetectorRNTuple> fDetectorWriters;
  
  /// Cluster size for batching
  UInt_t fClusterSize;
  
  /// Events accumulated in current cluster
  UInt_t fEventsInCurrentCluster;
  
  /// Enable batching
  Bool_t fEnableBatching;
  
  /// Has been initialized
  Bool_t fInitialized;
  
  /// Has been closed
  Bool_t fClosed;
};

#endif // HAS_RNTUPLE_SUPPORT

#endif // __QWPERDETECTORRNTUPLES__
