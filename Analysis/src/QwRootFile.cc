#include "QwRootFile.h"
#include "QwRNTupleFile.h"
#include "QwRunCondition.h"
#include "TH1.h"

#include <unistd.h>
#include <cstdio>

std::string QwRootFile::fDefaultRootFileDir = ".";
std::string QwRootFile::fDefaultRootFileStem = "Qweak_";

const Long64_t QwRootFile::kMaxTreeSize = 100000000000LL;
const Int_t QwRootFile::kMaxMapFileSize = 0x3fffffff; // 1 GiB

const TString QwRootTree::kUnitsName = "ppm/D:ppb/D:um/D:mm/D:mV_uA/D:V_uA/D";
Double_t QwRootTree::kUnitsValue[] = { 1e-6, 1e-9, 1e-3, 1 , 1e-3, 1};

/**
 * Constructor with relative filename
 */
QwRootFile::QwRootFile(const TString& run_label)
  : fRootFile(0), fMakePermanent(0),
    fMapFile(0), fEnableMapFile(kFALSE),
    fUpdateInterval(-1), fUseRNTuple(false), fRNTupleFile(nullptr), 
    fRunLabel(run_label), fOptions(&gQwOptions)
{
  // Process the configuration options
  ProcessOptions(gQwOptions);

#ifdef QW_ENABLE_MAPFILE
  // Check for the memory-mapped file flag
  if (fEnableMapFile) {

    TString mapfilename = "/dev/shm/";

    mapfilename += "/QwMemMapFile.map";

    fMapFile = TMapFile::Create(mapfilename,"UPDATE", kMaxMapFileSize, "RealTime Producer File");

    if (not fMapFile) {
      QwError << "Memory-mapped file " << mapfilename
              << " could not be opened!" << QwLog::endl;
      return;
    }

    QwMessage << "================== RealTime Producer Memory Map File =================" << QwLog::endl;
    fMapFile->Print();
    QwMessage << "======================================================================" << QwLog::endl;
  } else
#endif
  {

    TString rootfilename = fRootFileDir;
    TString hostname = gSystem -> HostName();

    // Use a probably-unique temporary file name.
    pid_t pid = getpid();

    fPermanentName = rootfilename
      + Form("/%s%s.root", fRootFileStem.Data(), run_label.Data());
    if (fUseTemporaryFile){
      rootfilename += Form("/%s%s.%s.%d.root",
			   fRootFileStem.Data(), run_label.Data(),
			   hostname.Data(), pid);
    } else {
      rootfilename = fPermanentName;
    }
    fRootFile = new TFile(rootfilename.Data(), "RECREATE", "myfile1");
    if (! fRootFile) {
      QwError << "ROOT file " << rootfilename
              << " could not be opened!" << QwLog::endl;
      return;
    } else {
      QwMessage << "Opened "<< (fUseTemporaryFile?"temporary ":"")
		<<"rootfile " << rootfilename << QwLog::endl;
    }

    TString run_condition_name = Form("condition_%s", run_label.Data());
    TList *run_cond_list = (TList*) fRootFile -> FindObjectAny(run_condition_name);
    if (not run_cond_list) {
      QwRunCondition run_condition(
          gQwOptions.GetArgc(),
          gQwOptions.GetArgv(),
          run_condition_name
      );

      fRootFile -> WriteObject(
          run_condition.Get(),
          run_condition.GetName()
      );
    }

    fRootFile->SetCompressionLevel(fCompressionLevel);
  }
}


/**
 * Destructor
 */
QwRootFile::~QwRootFile()
{
  // Keep the file on disk if any trees or histograms have been filled.
  // Also respect any other requests to keep the file around.
  if (!fMakePermanent) fMakePermanent = HasAnyFilled();

  // Close the map file
  if (fMapFile) {
    fMapFile->Close();
    // TMapFiles may not be deleted
    fMapFile = 0;
  }

  // Close the ROOT file.
  // Rename if permanence is requested, remove otherwise
  if (fRootFile) {
    TString rootfilename = fRootFile->GetName();

    fRootFile->Close();
    delete fRootFile;
    fRootFile = 0;

    int err;
    const char* action;
    if (fUseTemporaryFile){
      if (fMakePermanent) {
	action = " rename ";
	err = rename( rootfilename.Data(), fPermanentName.Data() );
      } else {
	action = " remove ";
	err = remove( rootfilename.Data() );
      }
      // It'd be proper to "extern int errno" and strerror() here,
      // but that doesn't seem very C++-ish.
      if (err) {
	QwWarning << "Couldn't" << action << rootfilename << QwLog::endl;
      } else {
	QwMessage << "Was able to" << action << rootfilename << QwLog::endl;
	QwMessage << "Root file is " << fPermanentName << QwLog::endl;
      }
    }
  }

  // Clean up RNTuple file wrapper
  if (fRNTupleFile) {
    delete fRNTupleFile;
    fRNTupleFile = nullptr;
  }

  // Delete Qweak ROOT trees
  std::map< const std::string, std::vector<QwRootTree*> >::iterator map_iter;
  std::vector<QwRootTree*>::iterator vec_iter;
  for (map_iter = fTreeByName.begin(); map_iter != fTreeByName.end(); map_iter++) {
    for (vec_iter = map_iter->second.begin(); vec_iter != map_iter->second.end(); vec_iter++) {
      delete *vec_iter;
    }
  }
}

/**
 * Defines configuration options using QwOptions functionality.
 * @param options Options object
 */
void QwRootFile::DefineOptions(QwOptions &options)
{
  // Define the ROOT files directory
  options.AddOptions("Default options")
    ("rootfiles", po::value<std::string>()->default_value(fDefaultRootFileDir),
     "directory of the output ROOT files");

  // Define the ROOT filename stem
  options.AddOptions("Default options")
    ("rootfile-stem", po::value<std::string>()->default_value(fDefaultRootFileStem),
     "stem of the output ROOT filename");

  // Define the memory map option
  options.AddOptions()
    ("enable-mapfile", po::value<bool>()->default_bool_value(false),
     "enable output to memory-mapped file\n(likely requires circular-buffer too)");
  options.AddOptions()
    ("write-temporary-rootfiles", po::value<bool>()->default_bool_value(true),
     "When writing ROOT files, use the PID to create a temporary filename");

  // Define the histogram and tree options
  options.AddOptions("ROOT output options")
    ("disable-tree", po::value<std::vector<std::string>>()->composing(),
     "disable output to tree regex");
  options.AddOptions("ROOT output options")
    ("disable-trees", po::value<bool>()->default_bool_value(false),
     "disable output to all trees");
  options.AddOptions("ROOT output options")
    ("disable-histos", po::value<bool>()->default_bool_value(false),
     "disable output to all histograms");

  // Define the helicity window versus helicity pattern options
  options.AddOptions("ROOT output options")
    ("disable-mps-tree", po::value<bool>()->default_bool_value(false),
     "disable helicity window output");
  options.AddOptions("ROOT output options")
    ("disable-pair-tree", po::value<bool>()->default_bool_value(false),
     "disable helicity pairs output");
  options.AddOptions("ROOT output options")
    ("disable-hel-tree", po::value<bool>()->default_bool_value(false),
     "disable helicity pattern output");
  options.AddOptions("ROOT output options")
    ("disable-burst-tree", po::value<bool>()->default_bool_value(false),
     "disable burst tree");
  options.AddOptions("ROOT output options")
    ("disable-slow-tree", po::value<bool>()->default_bool_value(false),
     "disable slow control tree");

  // Define the tree output prescaling options
  options.AddOptions("ROOT output options")
    ("num-mps-accepted-events", po::value<int>()->default_value(0),
     "number of accepted consecutive MPS events");
  options.AddOptions("ROOT output options")
    ("num-mps-discarded-events", po::value<int>()->default_value(0),
     "number of discarded consecutive MPS events");
  options.AddOptions("ROOT output options")
    ("num-hel-accepted-events", po::value<int>()->default_value(0),
     "number of accepted consecutive pattern events");
  options.AddOptions("ROOT output options")
    ("num-hel-discarded-events", po::value<int>()->default_value(0),
     "number of discarded consecutive pattern events");
  options.AddOptions("ROOT output options")
    ("mapfile-update-interval", po::value<int>()->default_value(-1),
     "Events between a map file update");

  // Define the autoflush and autosave option (default values by ROOT)
  options.AddOptions("ROOT performance options")
    ("autoflush", po::value<int>()->default_value(0),
     "TTree autoflush");
  options.AddOptions("ROOT performance options")
    ("autosave", po::value<int>()->default_value(300000000),
     "TTree autosave");
  options.AddOptions("ROOT performance options")
    ("basket-size", po::value<int>()->default_value(16000),
     "TTree basket size");
  options.AddOptions("ROOT performance options")
    ("circular-buffer", po::value<int>()->default_value(0),
     "TTree circular buffer");
  options.AddOptions("ROOT performance options")
    ("compression-level", po::value<int>()->default_value(1),
     "TFile compression level");

  // Define RNTuple options
  options.AddOptions("RNTuple options")
    ("enable-rntuple", po::value<bool>()->default_bool_value(false),
     "enable RNTuple output format instead of TTree");
  options.AddOptions("RNTuple options")
    ("disable-rntuple", po::value<std::vector<std::string>>()->composing(),
     "disable output to specific RNTuple regex");
  
  // Define RNTuple-specific options via the wrapper class
  QwRNTupleFile::DefineOptions(options);
}


/**
 * Parse the configuration options and store in class fields
 * @param options Options object
 */
void QwRootFile::ProcessOptions(QwOptions &options)
{
  // Option 'rootfiles' to specify ROOT files dir
  fRootFileDir = TString(options.GetValue<std::string>("rootfiles"));

  // Option 'root-stem' to specify ROOT file stem
  fRootFileStem = TString(options.GetValue<std::string>("rootfile-stem"));

  // Option 'mapfile' to enable memory-mapped ROOT file
  fEnableMapFile = options.GetValue<bool>("enable-mapfile");
#ifndef QW_ENABLE_MAPFILE
  if( fEnableMapFile ) {
    QwMessage << QwLog::endl;
    QwWarning << "QwRootFile::ProcessOptions:  "
              << "The 'enable-mapfile' flag is not supported by the ROOT "
                 "version with which this app is built. Disabling it."
              << QwLog::endl;
    fEnableMapFile = false;
  }
#endif
  fUseTemporaryFile = options.GetValue<bool>("write-temporary-rootfiles");

  // Options 'disable-trees' and 'disable-histos' for disabling
  // tree and histogram output
  auto v = options.GetValueVector<std::string>("disable-tree");
  std::for_each(v.begin(), v.end(), [&](const std::string& s){ this->DisableTree(s); });
  if (options.GetValue<bool>("disable-trees"))  DisableTree(".*");
  if (options.GetValue<bool>("disable-histos")) DisableHisto(".*");

  // Options 'disable-mps' and 'disable-hel' for disabling
  // helicity window and helicity pattern output
  if (options.GetValue<bool>("disable-mps-tree"))  DisableTree("^evt$");
  if (options.GetValue<bool>("disable-pair-tree"))  DisableTree("^pr$");
  if (options.GetValue<bool>("disable-hel-tree"))  DisableTree("^mul$");
  if (options.GetValue<bool>("disable-burst-tree"))  DisableTree("^burst$");
  if (options.GetValue<bool>("disable-slow-tree")) DisableTree("^slow$");

  // Options 'num-accepted-events' and 'num-discarded-events' for
  // prescaling of the tree output
  fNumMpsEventsToSave = options.GetValue<int>("num-mps-accepted-events");
  fNumMpsEventsToSkip = options.GetValue<int>("num-mps-discarded-events");
  fNumHelEventsToSave = options.GetValue<int>("num-mps-accepted-events");
  fNumHelEventsToSkip = options.GetValue<int>("num-mps-discarded-events");

  // Update interval for the map file
  fCircularBufferSize = options.GetValue<int>("circular-buffer");
  fUpdateInterval = options.GetValue<int>("mapfile-update-interval");
  fCompressionLevel = options.GetValue<int>("compression-level");
  fBasketSize = options.GetValue<int>("basket-size");

  // Autoflush and autosave
  fAutoFlush = options.GetValue<int>("autoflush");
  if ((ROOT_VERSION_CODE < ROOT_VERSION(5,26,00)) && fAutoFlush != -30000000){
    QwMessage << QwLog::endl;
    QwWarning << "QwRootFile::ProcessOptions:  "
              << "The 'autoflush' flag is not supported by ROOT version "
              << ROOT_RELEASE
              << QwLog::endl;
  }
  fAutoSave  = options.GetValue<int>("autosave");

  // Process RNTuple options
  fUseRNTuple = options.GetValue<bool>("enable-rntuple");
  
  // Process disabled RNTuples if RNTuple mode is enabled
  if (fUseRNTuple) {
    auto disabled_rntuple_v = options.GetValueVector<std::string>("disable-rntuple");
    std::for_each(disabled_rntuple_v.begin(), disabled_rntuple_v.end(), 
                  [&](const std::string& s){ this->DisableRNTuple(s); });
  }
  
  return;
}

/**
 * Determine whether the rootfile object has any non-empty trees or
 * histograms.
 */
Bool_t QwRootFile::HasAnyFilled(void) {
  return this->HasAnyFilled(fRootFile);
}
Bool_t QwRootFile::HasAnyFilled(TDirectory* d) {
  if (!d) return false;
  TList* l = d->GetListOfKeys();

  for( int i=0; i < l->GetEntries(); ++i) {
    const char* name = l->At(i)->GetName();
    TObject* obj = d->FindObjectAny(name);

    // Objects which can't be found don't count.
    if (!obj) continue;

    // Lists of parameter files, map files, and job conditions don't count.
    if ( TString(name).Contains("parameter_file") ) continue;
    if ( TString(name).Contains("mapfile") ) continue;
    if ( TString(name).Contains("_condition") ) continue;
    //  The EPICS tree doesn't count
    if ( TString(name).Contains("slow") ) continue;

    // Recursively check subdirectories.
    if (obj->IsA()->InheritsFrom( "TDirectory" ))
      if (this->HasAnyFilled( (TDirectory*)obj )) return true;

    if (obj->IsA()->InheritsFrom( "TTree" ))
      if ( ((TTree*) obj)->GetEntries() ) return true;

    if (obj->IsA()->InheritsFrom( "TH1" ))
      if ( ((TH1*) obj)->GetEntries() ) return true;
  }

  return false;
}

//=== RNTuple Method Implementations ===

void QwRootFile::NewRNTuple(const std::string& name, const std::string& desc) {
  // Only proceed if RNTuple mode is enabled
  if (!fUseRNTuple) return;
  
  if (IsRNTupleDisabled(name)) return;
  
  // Create RNTuple file wrapper if needed
  if (!fRNTupleFile) {
    fRNTupleFile = new QwRNTupleFile(fRunLabel);
    fRNTupleFile->ProcessOptions(*fOptions);
  }
  
  // Delegate to RNTuple file wrapper
  fRNTupleFile->NewRNTuple(name, desc);
}

Int_t QwRootFile::FillRNTuple(const std::string& name) {
  // Only proceed if RNTuple mode is enabled
  if (!fUseRNTuple || !fRNTupleFile) return 0;
  
  // Delegate to RNTuple file wrapper
  return fRNTupleFile->FillRNTuple(name);
}

Int_t QwRootFile::FillRNTuples() {
  // Only proceed if RNTuple mode is enabled
  if (!fUseRNTuple || !fRNTupleFile) return 0;
  
  // Delegate to RNTuple file wrapper
  return fRNTupleFile->FillRNTuples();
}

void QwRootFile::PrintRNTuples() const {
  // Only proceed if RNTuple mode is enabled
  if (!fUseRNTuple || !fRNTupleFile) {
    QwMessage << "RNTuple mode disabled or not initialized" << QwLog::endl;
    return;
  }
  
  // Delegate to RNTuple file wrapper
  fRNTupleFile->PrintRNTuples();
}
