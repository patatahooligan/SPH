#pragma once

class Video{
	private:
		AVCodecContext *context;
		AVFrame *frame;
		FILE *f;

	public:
		Video() :
			context(NULL),
			frame(NULL) {}
		~Video();
		void video_init();
		void video_finalize();
};
