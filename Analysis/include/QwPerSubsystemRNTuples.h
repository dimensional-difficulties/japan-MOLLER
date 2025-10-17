/**
 * QwPerSubsystemRNTuples.h
 *
 * Per-Subsystem RNTuple Manager
 * 
 * Groups detectors by subsystem (BPM, BCM, MainDetector, etc.) to reduce
 * the number of RNTuples from 430+ (per-detector) to ~10-15 (per-subsystem).
 * This dramatically reduces finalization time from ~15 minutes to <1 minute.
 *
 * Subsystems:
 * - BPM: Beam Position Monitors (bpm*, qwk_*bpm*)
 * - BCM: Beam Current Monitors (bcm*, qwk_bcm*)
 * - MainDetector: Main detector arrays (md*, pmtled*, qpd*)
 * - Combined: Combined devices (bpm_target, bcm_target, etc.)
 * - Energy: Energy calculators
 * - Metadata: Event metadata (CodaEventNumber, patterns, etc.)
 */

#ifndef QW_PER_SUBSYSTEM_RNTUPLES_H
#define QW_PER_SUBSYSTEM_RNTUPLES_H

#ifdef HAS_RNTUPLE_SUPPORT

#include "Rtypes.h"
#include "TFile.h"
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>

#include <string>
#include <vector>
#include <map>
#include <memory>

#include "QwLog.h"

/**
 * QwPerSubsystemRNTuples
 * 
 * Manages multiple RNTuples grouped by detector subsystem.
 * Each subsystem (BPM, BCM, MainDetector, etc.) gets one RNTuple
 * containing all fields from all detectors in that subsystem.
 */
class QwPerSubsystemRNTuples {
public:
  /**
   * Constructor
   * @param base_name Base name for RNTuples (e.g., "evts")
   * @param description Description for the RNTuple collection
   */
  QwPerSubsystemRNTuples(const std::string& base_name, const std::string& description);
  
  /**
   * Destructor
   */
  virtual ~QwPerSubsystemRNTuples();
  
  /**
   * Initialize with TFile
   * @param file ROOT file to write RNTuples to
   */
  void Initialize(TFile* file);
  
  /**
   * Register fields grouped by subsystem
   * @param field_names List of all field names to register
   */
  void RegisterFields(const std::vector<std::string>& field_names);
  
  /**
   * Fill data for all fields (called once per event)
   * @param values Map of field_name -> value
   */
  void Fill(const std::map<std::string, Double_t>& values);
  
  /**
   * Commit clusters for batching (called periodically)
   */
  void CommitAllClusters();
  
  /**
   * Close all RNTuples and finalize
   */
  void Close();
  
  /**
   * Set cluster size for batching
   */
  void SetClusterSize(UInt_t cluster_size);
  
  /**
   * Enable/disable batching
   */
  void SetBatching(Bool_t enable);
  
  /**
   * Extract subsystem name from field name
   * @param field_name Full field name (e.g., "bpm1c10WS_hw_sum")
   * @return Subsystem name (e.g., "BPM")
   */
  static std::string ExtractSubsystemName(const std::string& field_name);
  
  /**
   * Group fields by subsystem
   * @param field_names List of all field names
   * @return Map of subsystem_name -> list of field names in that subsystem
   */
  static std::map<std::string, std::vector<std::string>> 
    GroupFieldsBySubsystem(const std::vector<std::string>& field_names);

private:
  /**
   * SubsystemRNTuple - structure holding one subsystem's RNTuple data
   */
  struct SubsystemRNTuple {
    std::string name;                                      // Subsystem name
    std::unique_ptr<ROOT::RNTupleModel> model;            // RNTuple data model
    std::unique_ptr<ROOT::RNTupleWriter> writer;          // RNTuple writer
    std::vector<std::shared_ptr<Double_t>> field_ptrs;    // Pointers to fields
    std::vector<std::string> field_names;                  // Names of fields
    std::map<std::string, size_t> field_index_map;        // field_name -> index in field_ptrs
  };
  
  std::string fBaseName;                                   // Base name for RNTuples
  std::string fDescription;                                // Description
  TFile* fFile;                                            // ROOT file
  
  std::map<std::string, SubsystemRNTuple> fSubsystemWriters;  // subsystem_name -> RNTuple data
  
  UInt_t fClusterSize;                                     // Events per cluster
  UInt_t fEventsInCurrentCluster;                          // Event counter for batching
  Bool_t fEnableBatching;                                  // Whether batching is enabled
  Bool_t fInitialized;                                     // Initialization flag
  Bool_t fClosed;                                          // Close flag
  
  ClassDef(QwPerSubsystemRNTuples, 0);
};

#endif // HAS_RNTUPLE_SUPPORT

#endif // QW_PER_SUBSYSTEM_RNTUPLES_H
