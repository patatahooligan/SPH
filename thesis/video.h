#pragma once

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"
	#include "libavformat\avio.h"
}

class Video{
	private:
		AVCodecContext *codec_context;
		AVFormatContext *format_context;
		AVIOContext *io_context;
		AVFrame *frame;
		FILE *f;
		int current_frame;

	public:
		Video() :
			codec_context(NULL),
			format_context(NULL),
			io_context(NULL),
			frame(NULL),
			current_frame(0) {}
		~Video();
		void video_init();
		void video_finalize();
};
