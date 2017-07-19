#pragma once

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"
	#include "libavformat\avio.h"
	#include "libswscale\swscale.h"
}


const AVCodecID codec_id = AV_CODEC_ID_MPEG4;
const int framerate = 25;
const char filename[] = "C:/DevStuff/thesis/output.mp4";

class Video{
	// Class that copies openGL output into a video file. Uses ffmpeg(libav).
	private:
		AVCodecContext *codec_context;
		AVFormatContext *format_context;
		AVStream *ostream;
		AVFrame *rgbframe, *yuvframe;
		SwsContext *sws;
		FILE *f;
		int current_frame;
		bool finalized;

		int save_packets();
		// Grab packets from codec output and write them to the file. All contexts are expected to
		// have been initialized.

	public:
		Video() :
			codec_context(NULL),
			format_context(NULL),
			rgbframe(NULL),
			current_frame(0),
			finalized(false) {}
		~Video();

		// memory allocations & ffmpeg initialization
		void video_init();

		// Return true if a new frame needs to be rendered and saved given current simulation_time
		inline bool need_new_frame(float simulation_time) {return simulation_time >= current_frame / (float)framerate; }

		// Encodes current openGL output as a new frame if simulation_time warrants new frame. This check is done to ensure
		// that simulation time will be equal to real time when the video is played.
		void encode_frame(float simulation_time);

		// Finalizes the video file and de-allocates everything. Is called by the destructor so it will be moved to
		// private at some point when it serves no debugging purpose to call from outside the class.
		void video_finalize();
};
