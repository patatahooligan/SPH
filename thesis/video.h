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
		AVFrame *frame;
		AVPacket *pkt;
		FILE *f;
		int current_frame;

		int save_packets();
		// Grab packets from codec output and write them to the file. All contexts are expected to
		// have been initialized.

	public:
		Video() :
			codec_context(NULL),
			format_context(NULL),
			frame(NULL),
			pkt(NULL),
			current_frame(0) {}
		~Video();
		void video_init();
		void video_finalize();
};
