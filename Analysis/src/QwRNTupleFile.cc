#include "QwRNTupleFile.h"
#include "QwLog.h"
#include "QwOptions.h"

// Static member initialization
std::string QwRNTupleFile::fDefaultRootFileDir = "";
std::string QwRNTupleFile::fDefaultRootFileStem = "";

QwRNTupleFile::QwRNTupleFile(const TString& run_label)
: fRootFile(nullptr), fOwnsFile(kTRUE), fRunLabel(run_label.Data()), fUpdateInterval(0) {
  // Constructor implementation
}

QwRNTupleFile::QwRNTupleFile(const TString& run_label, TFile* existing_file)
: fRootFile(existing_file), fOwnsFile(kFALSE), fRunLabel(run_label.Data()), fUpdateInterval(0) {
  // Constructor that uses an existing file for same-file storage
  if (fRootFile) {
    QwMessage << "QwRNTupleFile: Using existing ROOT file: " << fRootFile->GetName() << QwLog::endl;
  }
}

QwRNTupleFile::~QwRNTupleFile() {
  // Clean up RNTuples
  for (auto& pair : fRNTupleByName) {
    for (auto* ntuple : pair.second) {
      delete ntuple;
    }
  }
  
  // Close file only if we own it
  if (fRootFile && fOwnsFile) {
    Close();
  } else if (fRootFile) {
    // For shared files, just clear our pointer but don't close
    QwMessage << "QwRNTupleFile: Not closing shared ROOT file in destructor" << QwLog::endl;
    fRootFile = nullptr;
  }
}

void QwRNTupleFile::DefineOptions(QwOptions &options) {
  // Define command line options for RNTuple files
  options.AddOptions()("QwRNTupleFile.disable-rntuple", 
                      po::value<std::vector<std::string>>()->multitoken(),
                      "Disable specific RNTuples by name");
  
  options.AddOptions()("QwRNTupleFile.update-interval", 
                      po::value<int>()->default_value(0),
                      "Update interval for histogram updates");
  
  options.AddOptions()("QwRNTupleFile.root-file-dir", 
                      po::value<std::string>(),
                      "Directory for ROOT files");
  
  options.AddOptions()("QwRNTupleFile.root-file-stem", 
                      po::value<std::string>(),
                      "Stem for ROOT file names");
}

void QwRNTupleFile::ProcessOptions(QwOptions &options) {
  // Process disabled RNTuples
  if (options.HasValue("QwRNTupleFile.disable-rntuple")) {
    std::vector<std::string> disabled = 
      options.GetValueVector<std::string>("QwRNTupleFile.disable-rntuple");
    for (const auto& name : disabled) {
      fDisabledRNTuples.insert(name);
      QwMessage << "RNTuple " << name << " disabled" << QwLog::endl;
    }
  }
  
  // Process update interval
  if (options.HasValue("QwRNTupleFile.update-interval")) {
    fUpdateInterval = options.GetValue<int>("QwRNTupleFile.update-interval");
  }
  
  // Process file directory and stem
  if (options.HasValue("QwRNTupleFile.root-file-dir")) {
    fDefaultRootFileDir = options.GetValue<std::string>("QwRNTupleFile.root-file-dir");
  }
  
  if (options.HasValue("QwRNTupleFile.root-file-stem")) {
    fDefaultRootFileStem = options.GetValue<std::string>("QwRNTupleFile.root-file-stem");
  }
  
  // Create the ROOT file only if not provided (separate files mode)
  if (!fRootFile) {
    std::string filename = fDefaultRootFileDir;
    if (!filename.empty() && filename.back() != '/') {
      filename += "/";
    }
    filename += fDefaultRootFileStem;
    if (!filename.empty()) {
      filename += "_";
    }
    filename += fRunLabel + "_rntuple.root";
    
    fRootFile = TFile::Open(filename.c_str(), "RECREATE");
    fOwnsFile = kTRUE; // We created this file, so we own it
    
    if (!fRootFile || fRootFile->IsZombie()) {
      QwError << "Failed to create ROOT file: " << filename << QwLog::endl;
      if (fRootFile) {
        delete fRootFile;
        fRootFile = nullptr;
      }
    } else {
      QwMessage << "Created RNTuple ROOT file: " << filename << QwLog::endl;
    }
  }
}

void QwRNTupleFile::CommitRNTuples() {
  // Commit all RNTuple writers and clear them to avoid double deletion
  if (fRNTupleWriters.empty()) {
    QwVerbose << "QwRNTupleFile::CommitRNTuples() - No writers to commit" << QwLog::endl;
    return;
  }
  
  for (auto& pair : fRNTupleWriters) {
    if (pair.second) {
      try {
        QwVerbose << "Committing RNTuple writer: " << pair.first << QwLog::endl;
        pair.second->CommitDataset();
      } catch (const std::exception& e) {
        QwError << "Exception committing RNTuple " << pair.first << ": " << e.what() << QwLog::endl;
      }
    }
  }
  // Clear the writers map to prevent further access
  fRNTupleWriters.clear();
}
