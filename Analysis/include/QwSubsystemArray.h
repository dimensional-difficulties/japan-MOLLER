/**********************************************************\
* File: QwSubsystemArray.h                                 *
*                                                          *
* Author: P. M. King,  Rakitha Beminiwattha                *
* Time-stamp: <2008-07-22 15:50>                           *
\**********************************************************/

#ifndef __QWSUBSYSTEMARRAY__
#define __QWSUBSYSTEMARRAY__

#include <vector>
#include <map>
#include "Rtypes.h"
#include "TString.h"
#include "TDirectory.h"
#include "TTree.h"


#include <boost/shared_ptr.hpp>
#include <boost/mem_fn.hpp>

// Qweak headers
#include "MQwPublishable.h"
#include "VQwSubsystem.h"
#include "QwOptions.h"

// Forward declarations
class VQwHardwareChannel;
class QwParameterFile;
class QwRNTuple;

///
/// \ingroup QwAnalysis
class QwSubsystemArray:
    public std::vector<boost::shared_ptr<VQwSubsystem>>,
    public MQwPublishable<QwSubsystemArray, VQwSubsystem> {

 private:
  typedef std::vector<boost::shared_ptr<VQwSubsystem> >  SubsysPtrs;
 public:
  using SubsysPtrs::const_iterator;
  using SubsysPtrs::iterator;
  using SubsysPtrs::begin;
  using SubsysPtrs::end;
  using SubsysPtrs::size;
  using SubsysPtrs::empty;

  typedef Bool_t (*CanContainFn)(VQwSubsystem*);

 private:

  /// Private default constructor
  QwSubsystemArray(); // not implemented, will throw linker error on use

 public:

  /// \brief Constructor with options
  QwSubsystemArray(QwOptions& options, CanContainFn myCanContain);
  /// \brief Copy constructor by reference
  QwSubsystemArray(const QwSubsystemArray& source);
  /// \brief Virtual destructor
  virtual ~QwSubsystemArray() { };

  /// \brief Assignment operator
  QwSubsystemArray& operator=(const QwSubsystemArray& value);

  /// \brief Set the internal record of the CODA run number
  void SetCodaRunNumber(UInt_t runnum) { fCodaRunNumber = runnum; };
  /// \brief Set the internal record of the CODA segment number
  void SetCodaSegmentNumber(UInt_t segnum) { fCodaSegmentNumber = segnum; };
  /// \brief Set the internal record of the CODA event number
  void SetCodaEventNumber(UInt_t evtnum) { fCodaEventNumber = evtnum; };
  /// \brief Set the internal record of the CODA event type
  void SetCodaEventType(UInt_t evttype) { fCodaEventType = evttype; };
  /// \brief Get the internal record of the CODA run number
  UInt_t GetCodaRunNumber() const { return fCodaRunNumber; };
  /// \brief Get the internal record of the CODA segment number
  UInt_t GetCodaSegmentNumber() const { return fCodaSegmentNumber; };
  /// \brief Get the internal record of the CODA event number
  UInt_t GetCodaEventNumber() const { return fCodaEventNumber; };
  /// \brief Get the internal record of the CODA event type
  UInt_t GetCodaEventType() const { return fCodaEventType; };

  /// \brief Set the internal record of the CODA event number
  void SetCleanParameters(Double_t cleanparameter[3])
  {
    fCleanParameter[0] = cleanparameter[0];
    fCleanParameter[1] = cleanparameter[1];
    fCleanParameter[2] = cleanparameter[2];
  };

  /// \brief Set event type mask
  void   SetEventTypeMask(const UInt_t mask) { fEventTypeMask = mask; };
  /// \brief Get event type mask
  UInt_t GetEventTypeMask() const { return fEventTypeMask; };
 /// \brief Update the event type mask from the subsystems
  UInt_t UpdateEventTypeMask() {
    for (iterator subsys_iter = begin(); subsys_iter != end(); ++subsys_iter) {
      VQwSubsystem* subsys = dynamic_cast<VQwSubsystem*>(subsys_iter->get());
      fEventTypeMask |= subsys->GetEventTypeMask();
    }
    return fEventTypeMask; 
  };


  /// \brief Set data loaded flag
  void   SetDataLoaded(const Bool_t flag) { fHasDataLoaded = flag; };
  /// \brief Get data loaded flag
  Bool_t HasDataLoaded() const { return fHasDataLoaded; };

  /// \brief Define configuration options for global array
  static void DefineOptions(QwOptions &options);
  /// \brief Process configuration options for the subsystem array itself
  void ProcessOptionsToplevel(QwOptions &options);
  /// \brief Process configuration options for all subsystems in the array
  void ProcessOptionsSubsystems(QwOptions &options);
  /// \brief Process configuration options (default behavior)
  void ProcessOptions(QwOptions &options) { ProcessOptionsSubsystems(options); };
  void LoadAllEventRanges(QwOptions &options);
  
  /// \brief Add the subsystem to this array
  void push_back(VQwSubsystem* subsys);

  /// \brief Get the subsystem with the specified name
  virtual VQwSubsystem* GetSubsystemByName(const TString& name);

  /// \brief Get the list of subsystems of the specified type
  virtual std::vector<VQwSubsystem*> GetSubsystemByType(const std::string& type);

  //each of the methods below will call their counterpart method separately.

  void  ClearEventData();

  /// \brief Process the event buffer for configuration events
  Int_t ProcessConfigurationBuffer(const ROCID_t roc_id, const BankID_t bank_id,
				   UInt_t *buffer, UInt_t num_words);

  /// \brief Process the event buffer for events
  Int_t ProcessEvBuffer(const UInt_t event_type, const ROCID_t roc_id,
			const BankID_t bank_id, UInt_t *buffer,
			UInt_t num_words);

  /// \brief Get the ROCID list
  void GetROCIDList(std::vector<ROCID_t> &list);

  /// \brief Process the decoded data in this event
  void  ProcessEvent();

  /// \brief Perform actions at the end of the event loop
  void  AtEndOfEventLoop();

 public:


  /// \brief Print list of parameter files
  void PrintParamFileList() const;

  /// \brief Get list of parameter files
  TList* GetParamFileNameList(TString name) const;

 public:

  /// \name Object construction and maintenance
  // @{
  /// Construct the objects for this subsystem
  void  ConstructObjects() {
    ConstructObjects((TDirectory*) NULL);
  };
  /// Construct the objects for this subsystem in a folder
  void  ConstructObjects(TDirectory *folder) {
    TString prefix = "";
    ConstructObjects(folder, prefix);
  };
  /// \brief Construct the objects for this subsystem in a folder with a prefix
  void  ConstructObjects(TDirectory *folder, TString &prefix);
  // @}


  /// \name Histogram construction and maintenance
  // @{
  /// Construct the histograms for this subsystem
  void  ConstructHistograms() {
    ConstructHistograms((TDirectory*) NULL);
  };
  /// Construct the histograms for this subsystem in a folder
  void  ConstructHistograms(TDirectory *folder) {
    TString prefix = "";
    ConstructHistograms(folder, prefix);
  };
  /// \brief Construct the histograms for this subsystem in a folder with a prefix
  void  ConstructHistograms(TDirectory *folder, TString &prefix);
  /// \brief Fill the histograms for this subsystem
  void  FillHistograms();
  /// \brief Share the histograms with another subsystem
  void  ShareHistograms(const QwSubsystemArray& source);
  // @}


  /// \name Tree and vector construction and maintenance
  // @{
  /// Construct the tree and vector for this subsystem
  void ConstructBranchAndVector(TTree *tree, std::vector <Double_t> &values) {
    TString tmpstr("");
    ConstructBranchAndVector(tree,tmpstr,values);
  };
  /// \brief Construct a branch and vector for this subsystem with a prefix
  void ConstructBranchAndVector(TTree *tree, TString& prefix, std::vector <Double_t> &values);
  /// \brief Construct a branch for this subsystem with a prefix
  void ConstructBranch(TTree *tree, TString& prefix);
  /// \brief Construct a branch for this subsystem with a prefix after tree leave trimming
  void ConstructBranch(TTree *tree, TString& prefix, QwParameterFile& trim_file);
  /// \brief Fill the vector for this subsystem
  void  FillTreeVector(std::vector<Double_t> &values) const;
  // @}

  /// \name RNTuple field construction and maintenance
  // @{
  /// \brief Construct the RNTuple fields for this subsystem with a prefix
  void ConstructRNTupleFields(QwRNTuple* rntuple, const TString& prefix);
  /// \brief Fill the RNTuple vector for this subsystem
  void FillRNTupleVector(std::vector<Double_t>& values) const;
  // @}


  /// \name Tree construction and maintenance
  /// These functions are not purely virtual, since not every subsystem is
  /// expected to implement them.  They are intended for expert output to
  /// trees.
  // @{
  /// Construct the tree for this subsystem
  void  ConstructTree() {
    ConstructTree((TDirectory*) NULL);
  };
  /// Construct the tree for this subsystem in a folder
  void  ConstructTree(TDirectory *folder) {
    TString prefix = "";
    ConstructTree(folder, prefix);
  };
  /// \brief Construct the tree for this subsystem in a folder with a prefix
  void  ConstructTree(TDirectory *folder, TString &prefix);

  /// \brief Fill the tree for this subsystem
  void  FillTree();
  /// \brief Delete the tree for this subsystem
  void  DeleteTree();
  // @}


  /// \brief Print some information about the subsystem
  void PrintInfo() const;
  
  void push_back(boost::shared_ptr<VQwSubsystem> subsys);

 protected:
  void LoadSubsystemsFromParameterFile(QwParameterFile& detectors);

 public:
  void GetMarkerWordList(const ROCID_t roc_id, const BankID_t bank_id, std::vector<UInt_t>& marker) const;

 protected:
  size_t fTreeArrayIndex;  //! Index of this data element in root tree

 protected:
  UInt_t fCodaRunNumber;     ///< CODA run number as provided by QwEventBuffer
  UInt_t fCodaSegmentNumber; ///< CODA segment number as provided by QwEventBuffer
  UInt_t fCodaEventNumber;   ///< CODA event number as provided by QwEventBuffer
  UInt_t fCodaEventType;     ///< CODA event type as provided by QwEventBuffer

  Double_t fCleanParameter[3];
  UInt_t fEventTypeMask;   ///< Mask of event types
  Bool_t fHasDataLoaded;   ///< Has this array gotten data to be processed?

 protected:
  /// Function to determine which subsystems we can accept
  CanContainFn fnCanContain;

  /// Test whether this subsystem array can contain a particular subsystem
  static Bool_t CanContain(VQwSubsystem* subsys) {
    if (subsys == 0) {
      QwError << "Zero pointer passed!" << QwLog::endl;
    }
    // should never occur
    return kFALSE;
  };

  std::vector< std::pair<UInt_t,UInt_t> > fBadEventRange; 

 private:
  /// Filename of the global detector map
  std::string fSubsystemsMapFile;
  std::vector<std::string> fSubsystemsDisabledByName; ///< List of disabled types
  std::vector<std::string> fSubsystemsDisabledByType; ///< List of disabled names

public:
  // Mock Data Variables
    /// \brief Randomize the data in this event
  void  RandomizeEventData(int helicity = 0, double time = 0.0);

  /// \brief Encode the data in this event
  void  EncodeEventData(std::vector<UInt_t> &buffer);

  double GetWindowPeriod(){return fWindowPeriod;};

protected:
  double fWindowPeriod;
  

}; // class QwSubsystemArray


#endif // __QWSUBSYSTEMARRAY__
