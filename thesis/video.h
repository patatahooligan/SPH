#pragma once

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"
	#include "libavformat\avio.h"
	#include "libswscale\swscale.h"
}

class Video{
	private:
		AVCodecContext *codec_context;
		AVFormatContext *format_context;
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
		void video_init();
		void encode_frame(float simulation_time);
		void video_finalize();
};
