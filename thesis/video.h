#pragma once

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"

class Video{
	private:
		AVCodecContext *context;
		AVFrame *frame;
		FILE *f;
		int current_frame;

	public:
		Video() :
			context(NULL),
			frame(NULL),
			current_frame(0) {}
		~Video();
		void video_init();
		void video_finalize();
};
