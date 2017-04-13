#include "stdafx.h"

#include <stdexcept>
#include <iostream>

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"
}

#include "video.h"
#include "constants.h"

const AVCodecID codec_id = AV_CODEC_ID_H264;
const int framerate = 25;

Video::~Video() {
	// Free all AVCodec structs that have been allocated while initializing.
	// Post a warning if the context is still open because it might mean that the file was
	// not finalized and created.

	if (context) {
		avcodec_free_context(&context);
		std::cerr << "WARNING : A Video object was destroyed with a live AVCodecContext. Possibly a video file was not finalized.";
	}

	av_frame_free(&frame);
}

void Video::video_init() {
	// Initialize libav

	// Find the desired codec
	AVCodec *codec = avcodec_find_encoder(codec_id);
	if (!codec) {
		// If avcodec_find_encoder returned null, the specified codec was not found
		throw std::runtime_error("Codec id not found by avcodec");
	}

	// Initialize a context which serves as an API
	context = avcodec_alloc_context3(codec);
	if (!context) {
		throw std::runtime_error("Could not allocate AVCodecContext");
	}

	// Set video dimensions. They match with the glut output window.
	context->width  = output_width;
	context->height = output_height;
	context->time_base.num = 1;
	context->time_base.den = framerate;
	// IMPORTANT : Not yet sure if context needs more params initialized!

	frame = av_frame_alloc();
	if (!frame) {
		throw std::runtime_error("Could not allocate AVFrame");
	}

	frame->format = context->pix_fmt;
	frame->width  = output_width;
	frame->height = output_height;
}

void Video::video_finalize() {
	// Currently a placeholder
	// TODO : add stuff that save the video file

	// It's important to destroy the context because the destructor checks if it exists
	// and posts a warning to catch situations where a Video object was destroyed without
	// finalizing the video capture process.
	avcodec_free_context(&context);
}