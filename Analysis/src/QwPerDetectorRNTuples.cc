#ifdef HAS_RNTUPLE_SUPPORT

#include "QwPerDetectorRNTuples.h"
#include <stdexcept>
#include <chrono>

/**
 * Constructor
 */
QwPerDetectorRNTuples::QwPerDetectorRNTuples(const std::string& base_name, const std::string& description)
  : fBaseName(base_name),
    fDescription(description),
    fFile(nullptr),
    fClusterSize(25000),
    fEventsInCurrentCluster(0),
    fEnableBatching(kTRUE),
    fInitialized(kFALSE),
    fClosed(kFALSE)
{
  QwMessage << "QwPerDetectorRNTuples: Creating manager for '" << fBaseName << "'" << QwLog::endl;
}

/**
 * Destructor
 */
QwPerDetectorRNTuples::~QwPerDetectorRNTuples()
{
  if (!fClosed) {
    Close();
  }
}

/**
 * Initialize with TFile
 */
void QwPerDetectorRNTuples::Initialize(TFile* file)
{
  if (fInitialized) {
    QwWarning << "QwPerDetectorRNTuples: Already initialized!" << QwLog::endl;
    return;
  }
  
  if (!file) {
    QwError << "QwPerDetectorRNTuples: NULL TFile pointer!" << QwLog::endl;
    return;
  }
  
  fFile = file;
  fInitialized = kTRUE;
  
  QwMessage << "QwPerDetectorRNTuples: Initialized with file " << file->GetName() << QwLog::endl;
}

/**
 * Register a detector with its fields
 */
void QwPerDetectorRNTuples::RegisterDetector(const std::string& detector_name, 
                                            const std::vector<std::string>& field_names)
{
  if (!fInitialized) {
    QwError << "QwPerDetectorRNTuples: Not initialized! Call Initialize() first." << QwLog::endl;
    return;
  }
  
  if (fDetectorWriters.find(detector_name) != fDetectorWriters.end()) {
    QwWarning << "QwPerDetectorRNTuples: Detector '" << detector_name << "' already registered!" << QwLog::endl;
    return;
  }
  
  // Create detector RNTuple structure
  DetectorRNTuple det;
  det.name = detector_name;
  det.field_names = field_names;
  det.model = ROOT::RNTupleModel::Create();
  
  // Create fields in the model
  for (const auto& field_name : field_names) {
    auto field_ptr = det.model->MakeField<Double_t>(field_name);
    det.field_ptrs.push_back(field_ptr);
  }
  
  // Create the RNTuple name: baseName_detectorName (e.g., "evts_bpm1c10WS")
  std::string rntuple_name = fBaseName + "_" + detector_name;
  
  try {
    // Create the writer
    det.writer = ROOT::RNTupleWriter::Append(std::move(det.model), rntuple_name, *fFile);
    
    // Store in map
    fDetectorWriters[detector_name] = std::move(det);
    
  } catch (const std::exception& e) {
    QwError << "QwPerDetectorRNTuples: Failed to create RNTuple for detector '" 
           << detector_name << "': " << e.what() << QwLog::endl;
  }
}

/**
 * Fill data for a specific detector
 */
void QwPerDetectorRNTuples::FillDetector(const std::string& detector_name, 
                                        const std::vector<Double_t>& values)
{
  auto it = fDetectorWriters.find(detector_name);
  if (it == fDetectorWriters.end()) {
    return;
  }
  
  DetectorRNTuple& det = it->second;
  
  // Check that value count matches field count
  if (values.size() != det.field_ptrs.size()) {
    QwError << "QwPerDetectorRNTuples: Value count mismatch for detector '" 
           << detector_name << "': expected " << det.field_ptrs.size() 
           << " but got " << values.size() << QwLog::endl;
    return;
  }
  
  // Copy values to field pointers
  for (size_t i = 0; i < values.size(); ++i) {
    if (det.field_ptrs[i]) {
      *(det.field_ptrs[i]) = values[i];
    }
  }
  
  // Fill the RNTuple
  if (det.writer) {
    det.writer->Fill();
  }
}

/**
 * Commit all RNTuple clusters (called periodically for batching)
 */
void QwPerDetectorRNTuples::CommitAllClusters()
{
  if (!fEnableBatching) {
    return;  // Batching is disabled
  }
  
  fEventsInCurrentCluster++;
  
  // Commit clusters when we reach the cluster size
  if (fEventsInCurrentCluster >= fClusterSize) {
    for (auto& pair : fDetectorWriters) {
      DetectorRNTuple& det = pair.second;
      if (det.writer) {
        det.writer->CommitCluster();
      }
    }
    
    fEventsInCurrentCluster = 0;
  }
}

/**
 * Close all RNTuples and finalize
 */
void QwPerDetectorRNTuples::Close()
{
  if (fClosed) {
    return;
  }
  
  QwMessage << "QwPerDetectorRNTuples::Close: Starting to close " << fDetectorWriters.size() 
           << " detector RNTuples for '" << fBaseName << "'" << QwLog::endl;
  
  // Commit any remaining events in current cluster
  if (fEnableBatching && fEventsInCurrentCluster > 0) {
    QwMessage << "QwPerDetectorRNTuples::Close: Committing final cluster with " 
             << fEventsInCurrentCluster << " events" << QwLog::endl;
    
    size_t count = 0;
    for (auto& pair : fDetectorWriters) {
      DetectorRNTuple& det = pair.second;
      if (det.writer) {
        det.writer->CommitCluster();
        count++;
        if (count % 100 == 0) {
          QwMessage << "QwPerDetectorRNTuples::Close: Committed cluster for " << count 
                   << "/" << fDetectorWriters.size() << " detectors" << QwLog::endl;
        }
      }
    }
    QwMessage << "QwPerDetectorRNTuples::Close: Finished committing clusters for " << count 
             << " detectors" << QwLog::endl;
  }
  
  // Close all writers
  QwMessage << "QwPerDetectorRNTuples::Close: Starting to finalize writers..." << QwLog::endl;
  size_t count = 0;
  for (auto& pair : fDetectorWriters) {
    DetectorRNTuple& det = pair.second;
    if (det.writer) {
      auto start_time = std::chrono::high_resolution_clock::now();
      
      QwMessage << "QwPerDetectorRNTuples::Close: Finalizing writer " << (count+1) 
               << "/" << fDetectorWriters.size() << " (" << pair.first << ")..." << QwLog::endl;
      
      det.writer.reset();  // Calls destructor which finalizes the RNTuple
      
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
      
      QwMessage << "QwPerDetectorRNTuples::Close: Writer " << (count+1) 
               << " (" << pair.first << ") finalized in " << duration.count() << " ms" << QwLog::endl;
      
      count++;
      if (count % 100 == 0) {
        QwMessage << "QwPerDetectorRNTuples::Close: Finalized " << count 
                 << "/" << fDetectorWriters.size() << " writers" << QwLog::endl;
      }
    }
  }
  
  fClosed = kTRUE;
  
  QwMessage << "QwPerDetectorRNTuples: Successfully closed " << fDetectorWriters.size() 
           << " detector RNTuples for '" << fBaseName << "'" << QwLog::endl;
}

/**
 * Set cluster size for all RNTuples
 */
void QwPerDetectorRNTuples::SetClusterSize(UInt_t cluster_size)
{
  fClusterSize = cluster_size;
  QwMessage << "QwPerDetectorRNTuples: Set cluster size to " << cluster_size << QwLog::endl;
}

/**
 * Enable/disable batching
 */
void QwPerDetectorRNTuples::SetBatching(Bool_t enable)
{
  fEnableBatching = enable;
  if (enable) {
    QwMessage << "QwPerDetectorRNTuples: Batching enabled with cluster size " 
             << fClusterSize << QwLog::endl;
  } else {
    QwMessage << "QwPerDetectorRNTuples: Batching disabled (less efficient)" << QwLog::endl;
  }
}

/**
 * Check if a detector is registered
 */
bool QwPerDetectorRNTuples::HasDetector(const std::string& detector_name) const
{
  return fDetectorWriters.find(detector_name) != fDetectorWriters.end();
}

/**
 * Extract detector name from field name
 * In TTrees, each branch is a detector, and leaves are fields within that detector.
 * In RNTuples, we flatten to: detectorname_fieldname
 * 
 * Examples:
 *   "bpm1c10WS_hw_sum" -> "bpm1c10WS"
 *   "bpm1c10WS_block0" -> "bpm1c10WS"
 *   "bpm1c10XP_sumsq0_high" -> "bpm1c10XP"
 *   "CodaEventNumber" -> "metadata" (no underscore, special case)
 */
std::string QwPerDetectorRNTuples::ExtractDetectorName(const std::string& field_name)
{
  // Special cases - metadata fields that don't belong to a detector
  if (field_name.find("Coda") == 0 || 
      field_name.find("pattern") == 0 ||
      field_name.find("time") == 0 ||
      field_name.find("delayed") == 0 ||
      field_name.find("reported") == 0 ||
      field_name.find("event_") == 0 ||
      field_name.find("input_") == 0 ||
      field_name.find("output_") == 0 ||
      field_name.find("mps_") == 0 ||
      field_name.find("pat_") == 0 ||
      field_name.find("cleandata") == 0 ||
      field_name.find("scandata") == 0) {
    return "metadata";
  }
  
  // Find the FIRST underscore - everything before it is the detector name
  // This matches TTree structure where branch name = detector, leaf name = field
  size_t first_underscore = field_name.find('_');
  
  if (first_underscore == std::string::npos) {
    // No underscore - this is likely a simple field name, treat as its own detector
    return field_name;
  }
  
  // Everything before the first underscore is the detector name
  return field_name.substr(0, first_underscore);
}

/**
 * Group fields by detector name
 */
std::map<std::string, std::vector<std::string>> QwPerDetectorRNTuples::GroupFieldsByDetector(
  const std::vector<std::string>& field_names)
{
  std::map<std::string, std::vector<std::string>> grouped;
  
  for (const auto& field_name : field_names) {
    std::string detector = ExtractDetectorName(field_name);
    grouped[detector].push_back(field_name);
  }
  
  return grouped;
}

#endif // HAS_RNTUPLE_SUPPORT
