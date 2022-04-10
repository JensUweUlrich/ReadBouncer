/*
 * GuppyBasecaller.hpp
 *
 *  Created on: 11.10.2021
 *      Author: juu
 */

#pragma once

#include <string>
#include <vector>
#include <memory>

namespace guppy_cpp_client
{
	class GuppyCPPClientImpl;

	enum class Result {
		success,                      ///< The method returned successfully.
		no_connection,                ///< The client is not connected, and the method cannot be called unless connected.
		already_connected,            ///< The client is already connected, and the method cannot be called while connected.
		not_ready,                    ///< The client is not able to accept any more reads (usually means the queue is full).
		timed_out,                    ///< The request has timed out.
		invalid_response,             ///< The server has returned an invalid response type for the request that was sent.
		basecall_config_unavailable,  ///< The requested basecall configuration is not available on the server.
		barcode_kit_unavailable,      ///< The requested barcoding configuration is not available on the server.
		align_index_unavailable,      ///< The requested alignment index is not available on the server.
		bed_file_unavailable,         ///< The requested bed file is not available on the server.
		align_index_no_seq,           ///< The requested alignment index does not contain the sequence(s), but full alignment was requested.
		file_load_timeout,            ///< The requested file could not be loaded on the server before timing out.
		failed,                       ///< Used to indicate general failure, for which no specific code is defined.
		bad_request,                  ///< The server could not interpret the request, either malformed or an incompatible protocol version.
		bad_reply,                    ///< The response received from the server was malformed and could not be interpreted.
	};

	enum class Priority : int {
		high = 0,     ///< Usually reserved for read-until applications.
		medium = 1,   ///< Standard priority for MinKNOW live basecalling.
		low = 2,      ///< Standalone applications running concurrently with MinKNOW should use this.
		unknown = -1  ///< This value should never occur. If it does, an error has occurred.
	};

	enum class ConnectionStatus {
		disconnected,              ///< The client is not connected to a server.
		preparing_connection,      ///< Performing tasks that must be done before connection can be established.
		loading_config,            ///< Currently waiting for a requested configuration to be loaded on the server.
		loading_alignment_index,   ///< Currently waiting for an alignment index to be loaded on the server.
		loading_bed_file,          ///< Currently waiting for a BED file to be loaded on the server.
		connecting,                ///< In the process of attempting to connect to a server.
		connected,                 ///< Currently connected to a server.
		error,                     ///< The client is in an error state.
	};

	enum class PassReadResult {
		success,                ///< The read has been successfully passed.
		not_accepting_reads,    ///< Reads are not being accepted. This could be because the client is currently inactive.
								///< (status is 'disconnected' or 'error') or finish() has been called and the client is
								///< accepting no more reads.
		queue_full,             ///< The read could not be accepted at this time as the input queue is full and
								///< no more reads will be accepted until some have been processed.
		bad_read,               ///< The read is malformed and will not be accepted for processing.
	};

	struct GuppyRead
	{
		uint64_t read_tag;
		std::string read_id;
		std::vector<int16_t> signal;
	};

	struct CalledRead
	{
		uint64_t read_tag;
		std::string read_id;
		std::string sequence;
		std::string qstring;
		uint16_t seqlen;
	};

	class GuppyCPPClient
	{
		private:
			std::unique_ptr<GuppyCPPClientImpl> m_impl;

		public:
			GuppyCPPClient(std::string& address, std::string& config);
			/// Destructor.
			~GuppyCPPClient();

			/// Copying of object is not allowed.
			GuppyCPPClient(const GuppyCPPClient&) = delete;

			/// Copying of object is not allowed.
			GuppyCPPClient& operator=(const GuppyCPPClient&) = delete;

			/** Specify a name for this client, this name will be reported in the guppy basecall server log alongside its
			*  server assigned ID
			*  @remarks sent as part of the connection request so this must be set prior to connecting.
			*/
			Result set_name(std::string& name);

			/// Specify the timeout (in milliseconds) to use for all queries to the server.
			Result set_timeout(uint32_t timeout);

			/** Specify a default timeout (in milliseconds) that the pass_read method should use when waiting to queue
			 *  a read to be passed to the caller.
			 *  This can be used instead of providing the optional timeout parameter to @see pass_read.
			 *
			 *  @remarks Default is zero, which means that if pass_read is called without a wait duration it will not wait
			 *  i.e. it will be non-blocking.
			 */
			Result set_default_pass_read_timeout_ms(uint32_t default_pass_read_timeout_ms);

			/** Specify the timeout (in seconds) for connecting to the server, or reconnecting in the event of a lost
			*  connection.
			*  If a server connection cannot be established before timing out the client will become inactive and enter
			*  an error state.
			*
			*  @remark The default is 300 seconds.
			*
			*  @remark Note that this represents a minimum timeout period, since each failed call to the server takes several
			*  seconds, the actual time taken trying to (re)connect will be at least this value.
			*/
			Result set_reconnect_timeout(uint32_t reconnect_timeout);

			/** Specify the maximum number of reads to queue on the client for sending to the server. Reads are sent to
			*  the server asynchronously. Setting this too large will consume a lot of memory on the client (particularly
			*  if reads are very long), but setting it too low may negatively impact performance. If the queue is full
			*  the *see pass_read method will return a PassReadResult::queue_full result.
			*/
			Result set_max_reads_queued(size_t max_reads_queued);

			/** Specify the priority the client should use. Each priority has its own set of queues on the server, and
			 *  higher priority reads will get a larger share of resources if the server is unable to keep up. Typically
			 *  the high priority should be used for Read Until applications, the medium priority for live MinKNOW
			 *  basecalling, and the low priority for standalone basecalling, either through the Guppy C++ application
			 *  or via the python interface.
			 */
			Result set_priority(Priority& priority);

			/** Specify whether Move and Trace tables should be returned by the server. If Move and Trace data are not
			*  going to be used by the client, then disabling this will reduce the interprocess communication overhead.
			*/
			Result set_move_and_trace_enabled(bool enabled);

			/** Specify whether full State data should be returned by the server. If the State data is not going to be
			 *  used by the client, then this should be left disabled. Enabling it will greatly increase interprocess
			 *  communication overhead.
			 */
			Result set_state_data_enabled(bool enabled);

			/// Get the status of the client.
			ConnectionStatus get_status() const;

			/// Get the error message, if any.
			std::string get_error_message() const;

			/** For debugging/diagnostics only. Returns a string describing the internal state of the client.
			 *  @remarks it may be useful to capture this information if diagnosing an issue.
			 */
			std::string get_internal_state() const;

			/** Connect to basecall server.
			 *  @return A code indicating either successful connection, or indicating what went wrong.
			 *
			 *  All results other than success put the object into the 'error' state, get_error_message() should be
			 *  called for further details.
			 *
			 *  @remarks If a bad_request is returned it is most likely this is due to incompatible versions of of the
			 *  messaging protocol used by client and server, get_error_message() will have details of any version issues.
			 *
			 *  @remarks The client will be reset, any reads held due to an earlier connection will be lost.
			 */
			Result connect();

			/** Queue a read to be passed to the caller.
			 *
			*  @param read A Read object containing raw data for a read.
			*
			*  @param wait_duration_ms (optional) Duration to wait to try to pass the read. If the input queue is full this
			*  call will block until either space is available on the queue or the timeout period is reached. A value of
			*  zero is valid and means the call will not block.
			*
			*  @return A code indicating success, or why the attempt was unsuccessful.
			*
			*  The Read object must be properly initialised. This includes the following:
			*  - The read tag must be set to a value unique to this client. The value must fit in a 32 bit unsigned int.
			*  - The priority can be ignored. This will be set automatically based on the priority of the client.
			*  - The datasets must contain the entry dataset::raw_data, which must be a 1d array of signed 16 bit
			*    integers representing the raw signal.
			*  - The metadata must contain the following entries:
			*    - read_data::read_id -> The string formatted read-id for the read. This is mostly for debugging, so the
			*      the exact value passed is not important, nor does it need to be unique. But it should not be empty.
			*    - read_data::daq_offset -> The offset to add to the daq values for conversion to pA.
			*    - read_data::daq_scaling -> The scale factor to multiply the daq values by for conversion to pA.
			*
			*  Any other fields in the metadata will be ignored by the client.
			*
			*  Note that this function uses move semantics to take ownership of the read. After calling,
			*  if the read was accepted, then the passed read object will be blank and the caller will
			*  have taken ownership of its data. If the read is not accepted then the passed read object
			*  will be unaffected.
			*
			*  @@returns Any return code other than success means the read has not been accepted.
			*/
			PassReadResult pass_read(GuppyRead& read, uint32_t wait_duration_ms = 0);


			/** Get completed reads.
			 *  @param[in/out] reads An empty list to be populated with Read objects.
			 *
			 *  If this method is called after the object has gone into the 'error' state, then any reads that
			 *  have come back from the server, but not been retrieved yet, will be provided.
			 */
			void get_completed_reads(std::vector<CalledRead>& calledReads);

			/** Blocks until all reads are finished processing.
			 *  @param[in] timeout The timeout (in seconds) to wait for reads to finish. Zero means no timeout.
			 *  @return A code indicating success, or why the attempt was unsuccessful.
			 *
			 *  @@returns
			 *  - success
			 *  - no_connection : the client is in the 'disconnected' state.
			 *  - failed        : the client is in an error state, call get_error_message() to retrieve further details.
			 *  - not_ready     : the timeout was reached before all reads could be processed. The client is now in the
			 *                    'disconnected' state.
			 */
			Result finish(uint64_t timeout = 0);

			/** Puts the client into the 'disconnected' state. This will disconnect from the server if required and clear any
			 *  error state message.
			 *
			 *  @return returns 'Result.no_connection' if already in the 'disconnected' state otherwise returns 'Result.success'
			 *
			 *  @remark N.B. this will lose any error status information. Any completed reads may still be retrieved by
			 *  calling get_completed_reads.
			 *
			 *  @remark if wishing to reconnect after after an error, first call this method to enter the 'disconnected' state
			 *  prior to calling connect()
			 */
			Result disconnect();
	}
	;

} //namespace

