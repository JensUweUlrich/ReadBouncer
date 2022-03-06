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
