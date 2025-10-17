#ifdef HAS_RNTUPLE_SUPPORT

#include "QwPerSubsystemRNTuples.h"
#include <ROOT/RNTupleWriter.hxx>
#include <stdexcept>
#include <chrono>
#include <algorithm>
#include <cctype>

/**
 * Constructor
 */
QwPerSubsystemRNTuples::QwPerSubsystemRNTuples(const std::string& base_name, const std::string& description)
  : fBaseName(base_name),
    fDescription(description),
    fFile(nullptr),
    fClusterSize(25000),
    fEventsInCurrentCluster(0),
    fEnableBatching(kTRUE),
    fInitialized(kFALSE),
    fClosed(kFALSE)
{
  QwMessage << "QwPerSubsystemRNTuples: Creating manager for '" << fBaseName << "'" << QwLog::endl;
}

/**
 * Destructor
 */
QwPerSubsystemRNTuples::~QwPerSubsystemRNTuples()
{
  if (!fClosed) {
    Close();
  }
}

/**
 * Initialize with TFile
 */
void QwPerSubsystemRNTuples::Initialize(TFile* file)
{
  if (fInitialized) {
    QwWarning << "QwPerSubsystemRNTuples: Already initialized!" << QwLog::endl;
    return;
  }
  
  if (!file) {
    QwError << "QwPerSubsystemRNTuples: NULL TFile pointer!" << QwLog::endl;
    return;
  }
  
  fFile = file;
  fInitialized = kTRUE;
  
  QwMessage << "QwPerSubsystemRNTuples: Initialized with file " << file->GetName() << QwLog::endl;
}

/**
 * Extract subsystem name from field name
 * 
 * Logic:
 * - bpm*, qwk_*bpm* -> BPM
 * - bcm*, qwk_bcm* -> BCM
 * - md*, pmtled*, qpd* -> MainDetector
 * - *_target, *_energy -> Combined
 * - Coda*, pattern*, time*, event_*, etc. -> Metadata
 * - Everything else -> Other
 */
std::string QwPerSubsystemRNTuples::ExtractSubsystemName(const std::string& field_name)
{
  // Convert to lowercase for comparison
  std::string lower_name = field_name;
  std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  
  // Metadata fields
  if (lower_name.find("coda") == 0 || 
      lower_name.find("pattern") == 0 ||
      lower_name.find("time") == 0 ||
      lower_name.find("delayed") == 0 ||
      lower_name.find("reported") == 0 ||
      lower_name.find("event_") == 0 ||
      lower_name.find("input_") == 0 ||
      lower_name.find("output_") == 0 ||
      lower_name.find("mps_") == 0 ||
      lower_name.find("pat_") == 0 ||
      lower_name.find("cleandata") == 0 ||
      lower_name.find("scandata") == 0 ||
      lower_name.find("error") == 0 ||
      lower_name.find("sequence") == 0 ||
      lower_name.find("yield") == 0 ||
      lower_name.find("asym") == 0 ||
      lower_name.find("diff") == 0) {
    return "Metadata";
  }
  
  // Combined/derived quantities
  if (lower_name.find("_target") != std::string::npos ||
      lower_name.find("_energy") != std::string::npos ||
      lower_name.find("_targ") != std::string::npos ||
      lower_name.find("combo") == 0) {
    return "Combined";
  }
  
  // BPM (Beam Position Monitors)
  if (lower_name.find("bpm") != std::string::npos ||
      lower_name.find("qwk_0") == 0 ||  // QWK_0R06, QWK_0L06, etc.
      lower_name.find("qwk_1") == 0) {  // QWK_1C12, QWK_1H04, etc.
    return "BPM";
  }
  
  // BCM (Beam Current Monitors)
  if (lower_name.find("bcm") != std::string::npos ||
      lower_name.find("unser") != std::string::npos) {
    return "BCM";
  }
  
  // Main Detector
  if (lower_name.find("md") == 0 ||  // md detectors
      lower_name.find("pmtled") != std::string::npos ||
      lower_name.find("qpd") != std::string::npos ||
      lower_name.find("det") == 0 ||
      lower_name.find("usl") != std::string::npos ||  // upstream left
      lower_name.find("usr") != std::string::npos ||  // upstream right
      lower_name.find("dsl") != std::string::npos ||  // downstream left
      lower_name.find("dsr") != std::string::npos) {  // downstream right
    return "MainDetector";
  }
  
  // Lumi monitors
  if (lower_name.find("lumi") != std::string::npos ||
      lower_name.find("sam") != std::string::npos) {
    return "Luminosity";
  }
  
  // Slow controls / EPICS
  if (lower_name.find("ioc") == 0 ||
      lower_name.find("hac") == 0 ||
      lower_name.find("ips") == 0) {
    return "SlowControls";
  }
  
  // Default: put in "Other" subsystem
  return "Other";
}

/**
 * Group fields by subsystem
 */
std::map<std::string, std::vector<std::string>> 
QwPerSubsystemRNTuples::GroupFieldsBySubsystem(const std::vector<std::string>& field_names)
{
  std::map<std::string, std::vector<std::string>> grouped;
  
  for (const auto& field_name : field_names) {
    std::string subsystem = ExtractSubsystemName(field_name);
    grouped[subsystem].push_back(field_name);
  }
  
  return grouped;
}

/**
 * Register fields grouped by subsystem
 */
void QwPerSubsystemRNTuples::RegisterFields(const std::vector<std::string>& field_names)
{
  if (!fInitialized) {
    QwError << "QwPerSubsystemRNTuples: Not initialized! Call Initialize() first." << QwLog::endl;
    return;
  }
  
  // Group fields by subsystem
  auto grouped_fields = GroupFieldsBySubsystem(field_names);
  
  QwMessage << "QwPerSubsystemRNTuples: Registering " << field_names.size() 
           << " fields across " << grouped_fields.size() << " subsystems" << QwLog::endl;
  
  // Create an RNTuple for each subsystem
  for (const auto& pair : grouped_fields) {
    const std::string& subsystem_name = pair.first;
    const std::vector<std::string>& subsystem_fields = pair.second;
    
    if (fSubsystemWriters.find(subsystem_name) != fSubsystemWriters.end()) {
      QwWarning << "QwPerSubsystemRNTuples: Subsystem '" << subsystem_name 
               << "' already registered!" << QwLog::endl;
      continue;
    }
    
    QwMessage << "QwPerSubsystemRNTuples: Creating RNTuple for subsystem '" << subsystem_name 
             << "' with " << subsystem_fields.size() << " fields" << QwLog::endl;
    
    // Create subsystem RNTuple structure
    SubsystemRNTuple subsys;
    subsys.name = subsystem_name;
    subsys.field_names = subsystem_fields;
    subsys.model = ROOT::RNTupleModel::Create();
    
    // Create fields in the model
    size_t field_idx = 0;
    for (const auto& field_name : subsystem_fields) {
      auto field_ptr = subsys.model->MakeField<Double_t>(field_name);
      subsys.field_ptrs.push_back(field_ptr);
      subsys.field_index_map[field_name] = field_idx++;
    }
    
    // Create the RNTuple name: baseName_subsystemName (e.g., "evts_BPM")
    std::string rntuple_name = fBaseName + "_" + subsystem_name;
    
    try {
      // Create the writer
      subsys.writer = ROOT::RNTupleWriter::Append(std::move(subsys.model), rntuple_name, *fFile);
      
      // Store in map
      fSubsystemWriters[subsystem_name] = std::move(subsys);
      
      QwMessage << "QwPerSubsystemRNTuples: Created RNTuple '" << rntuple_name << "'" << QwLog::endl;
      
    } catch (const std::exception& e) {
      QwError << "QwPerSubsystemRNTuples: Failed to create RNTuple for subsystem '" 
             << subsystem_name << "': " << e.what() << QwLog::endl;
    }
  }
  
  QwMessage << "QwPerSubsystemRNTuples: Successfully created " << fSubsystemWriters.size() 
           << " subsystem RNTuples" << QwLog::endl;
}

/**
 * Fill data for all fields (called once per event)
 */
void QwPerSubsystemRNTuples::Fill(const std::map<std::string, Double_t>& values)
{
  // Group values by subsystem and fill each subsystem's RNTuple
  std::map<std::string, bool> subsystem_filled;
  
  for (const auto& value_pair : values) {
    const std::string& field_name = value_pair.first;
    Double_t value = value_pair.second;
    
    // Determine which subsystem this field belongs to
    std::string subsystem = ExtractSubsystemName(field_name);
    
    auto it = fSubsystemWriters.find(subsystem);
    if (it == fSubsystemWriters.end()) {
      // Subsystem not registered, skip
      continue;
    }
    
    SubsystemRNTuple& subsys = it->second;
    
    // Find the field in this subsystem
    auto field_it = subsys.field_index_map.find(field_name);
    if (field_it == subsys.field_index_map.end()) {
      // Field not in this subsystem, skip
      continue;
    }
    
    size_t field_idx = field_it->second;
    
    // Set the field value
    if (field_idx < subsys.field_ptrs.size() && subsys.field_ptrs[field_idx]) {
      *(subsys.field_ptrs[field_idx]) = value;
    }
    
    // Mark this subsystem as having data for this event
    subsystem_filled[subsystem] = true;
  }
  
  // Fill all subsystems that have data
  for (auto& pair : fSubsystemWriters) {
    const std::string& subsystem_name = pair.first;
    SubsystemRNTuple& subsys = pair.second;
    
    if (subsys.writer && subsystem_filled[subsystem_name]) {
      subsys.writer->Fill();
    }
  }
}

/**
 * Commit all RNTuple clusters (called periodically for batching)
 */
void QwPerSubsystemRNTuples::CommitAllClusters()
{
  if (!fEnableBatching) {
    return;  // Batching is disabled
  }
  
  fEventsInCurrentCluster++;
  
  // Commit clusters when we reach the cluster size
  if (fEventsInCurrentCluster >= fClusterSize) {
    for (auto& pair : fSubsystemWriters) {
      SubsystemRNTuple& subsys = pair.second;
      if (subsys.writer) {
        subsys.writer->CommitCluster();
      }
    }
    
    fEventsInCurrentCluster = 0;
  }
}

/**
 * Close all RNTuples and finalize
 */
void QwPerSubsystemRNTuples::Close()
{
  if (fClosed) {
    return;
  }
  
  QwMessage << "QwPerSubsystemRNTuples::Close: Starting to close " << fSubsystemWriters.size() 
           << " subsystem RNTuples for '" << fBaseName << "'" << QwLog::endl;
  
  // Commit any remaining events in current cluster
  if (fEnableBatching && fEventsInCurrentCluster > 0) {
    QwMessage << "QwPerSubsystemRNTuples::Close: Committing final cluster with " 
             << fEventsInCurrentCluster << " events" << QwLog::endl;
    
    for (auto& pair : fSubsystemWriters) {
      SubsystemRNTuple& subsys = pair.second;
      if (subsys.writer) {
        subsys.writer->CommitCluster();
      }
    }
  }
  
  // Close all writers
  QwMessage << "QwPerSubsystemRNTuples::Close: Starting to finalize writers..." << QwLog::endl;
  auto overall_start = std::chrono::high_resolution_clock::now();
  
  size_t count = 0;
  for (auto& pair : fSubsystemWriters) {
    SubsystemRNTuple& subsys = pair.second;
    if (subsys.writer) {
      auto start_time = std::chrono::high_resolution_clock::now();
      
      QwMessage << "QwPerSubsystemRNTuples::Close: Finalizing subsystem " << (count+1) 
               << "/" << fSubsystemWriters.size() << " (" << pair.first 
               << " with " << subsys.field_names.size() << " fields)..." << QwLog::endl;
      
      subsys.writer.reset();  // Calls destructor which finalizes the RNTuple
      
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
      
      QwMessage << "QwPerSubsystemRNTuples::Close: Subsystem " << (count+1) 
               << " (" << pair.first << ") finalized in " << duration.count() << " ms" << QwLog::endl;
      
      count++;
    }
  }
  
  auto overall_end = std::chrono::high_resolution_clock::now();
  auto overall_duration = std::chrono::duration_cast<std::chrono::milliseconds>(overall_end - overall_start);
  
  fClosed = kTRUE;
  
  QwMessage << "QwPerSubsystemRNTuples: Successfully closed " << fSubsystemWriters.size() 
           << " subsystem RNTuples for '" << fBaseName << "' in " 
           << overall_duration.count() << " ms total" << QwLog::endl;
}

/**
 * Set cluster size for all RNTuples
 */
void QwPerSubsystemRNTuples::SetClusterSize(UInt_t cluster_size)
{
  fClusterSize = cluster_size;
  QwMessage << "QwPerSubsystemRNTuples: Set cluster size to " << cluster_size << QwLog::endl;
}

/**
 * Enable/disable batching
 */
void QwPerSubsystemRNTuples::SetBatching(Bool_t enable)
{
  fEnableBatching = enable;
  if (enable) {
    QwMessage << "QwPerSubsystemRNTuples: Batching enabled with cluster size " 
             << fClusterSize << QwLog::endl;
  } else {
    QwMessage << "QwPerSubsystemRNTuples: Batching disabled (less efficient)" << QwLog::endl;
  }
}

#endif // HAS_RNTUPLE_SUPPORT
