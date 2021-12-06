/*
 * Acquisition.cpp
 *
 *  Created on: 07.04.2020
 *      Author: jens-uwe.ulrich
 */

#include "Acquisition.hpp"

namespace readuntil
{

    Acquisition::Acquisition(std::shared_ptr<::grpc::Channel> channel)
    {
        stub = AcquisitionService::NewStub(channel);
//        acquisition_logger = spdlog::get("RUClientLog");
    }

    bool Acquisition::hasStarted()
    {
        grpc::Status status;
        WatchForStatusChangeResponse response;

        stream = stub->watch_for_status_change(&context);
        // first write setup the stream
//        acquisition_logger->debug("reading Minknow status");
        while (stream->Read(&response))
        {
            
            if(response.status() == MinknowStatus::PROCESSING)
            {
                
                WatchForStatusChangeRequest request{};
                request.set_stop(true);
                if (!stream->Write(request))
                {
                    throw ReadUntilClientException("Unable to stop status acquisition!");
                }
                stream->WritesDone();
            }
        }
//        acquisition_logger->debug("Minknow status changed to PROCESSING");
	    context.TryCancel();
        status = stream->Finish();

//        std::cout << "Device is performing MUX scan. Wainting to finish..." << std::endl;

//        std::this_thread::sleep_for(std::chrono::seconds(150));

        return true;
    }

    bool Acquisition::isFinished()
    {
        CurrentStatusRequest request;
		CurrentStatusResponse response;
		::grpc::ClientContext context;
		::grpc::Status status = stub->current_status(&context, request, &response);
		if (status.ok())
		{
            if(response.status() == MinknowStatus::FINISHING)
            {
                return true;
            }
            else
            {
                return false;
            }
		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
    }

}

