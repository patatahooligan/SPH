#pragma once

class Video{
	private:
		AVCodecContext *context;

	public:
		Video() : context(NULL) {};
		~Video();
		void video_init();
		void video_finalize();
};
