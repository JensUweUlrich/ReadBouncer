#include "gmock/gmock.h"

using namespace interleave; 

class MockIBF : public interleave::IBF {

 public:

 // override only if we overriding a virtual method
 // MOCK_METHOD(return value, methods general name, (parameters), (override)); or (const) instead of override or if both (const, override)
  MOCK_METHOD(std::future< void >, parse_ref_seqs, (SafeQueue< Seqs > &queue_refs, interleave::IBFConfig &config, FilterStats &stats));
  MOCK_METHOD(void, add_sequences_to_filter, (std::vector< std::future< void > >& tasks, IBFConfig &config, uint64_t &binid, SafeQueue< Seqs > &queue_refs));
  MOCK_METHOD(std::vector<std::string>, cutOutNNNs, (std::string& seq, uint64_t seqlen));
  MOCK_METHOD(uint64_t, calculate_filter_size_bits, (IBFConfig& config, uint64_t numberOfBins));
  
  MOCK_METHOD(FilterStats, create_filter, (IBFConfig& config));
  MOCK_METHOD(FilterStats, load_filter, (IBFConfig& config));
  MOCK_METHOD(FilterStats, update_filter, (IBFConfig& config));
  MOCK_METHOD(TIbf, getFilter, ());

};

class MockRead : public interleave::Read {

 public:


 // override only if we overriding a virtual method
 // MOCK_METHOD(return value, methods general name, (parameters), (override)); or (const) instead of override or if both (const, override)

 // Private Methods: 
  MOCK_METHOD(bool, find_matches, (std::vector< TIbf >& filters, ClassifyConfig&  config ));
  MOCK_METHOD(uint64_t, count_matches, (IBFMeta& filter, ClassifyConfig& config));

  MOCK_METHOD(bool, select_matches, (std::vector< uint16_t >& selectedBins,
                                     std::vector< uint16_t >& selectedBinsRev,
                                     TIbf&                    filter,
                                     uint16_t                 threshold));
  MOCK_METHOD(uint64_t, max_matches, (std::vector< uint16_t >& selectedBins, std::vector< uint16_t >& selectedBinsRev,
                TIbf& filter, uint16_t threshold));

  // Public Methods: 

  MOCK_METHOD(uint32_t, getReadLength, ());
  MOCK_METHOD(bool, classify, (std::vector< TIbf >& filters, ClassifyConfig& config));
  MOCK_METHOD(int, classify, (std::vector< IBFMeta >& filters, ClassifyConfig& config));
  MOCK_METHOD((std::pair<int, int>), classify, (std::vector< IBFMeta >& filt1, std::vector< IBFMeta >& filt2, ClassifyConfig& config));



};

//class MockConfig  : public interleave::ClassifyConfig {};
