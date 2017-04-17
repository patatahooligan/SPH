#include "stdafx.h"

#include <stdexcept>
#include <iostream>

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"
	#include "libavformat\avio.h"
}

#include "video.h"
#include "constants.h"


const AVCodecID codec_id = AV_CODEC_ID_H264;
const int framerate = 25;


Video::~Video() {
	// Free all AVCodec structs that have been allocated while initializing.
	// Post a warning if the context is still open because it might mean that the file was
	// not finalized and created.

	if (codec_context) {
		avcodec_free_context(&codec_context);
		std::cerr << "WARNING : A Video object was destroyed with a live AVCodecContext. Possibly a video file was not finalized.";
	}

	if (format_context) {
		avformat_free_context(format_context);
		std::cerr << "WARNING : A Video object was destroyed with a live AVFormatContext. Possibly a video file was not finalized.";
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

	// Initialize a codec context
	codec_context = avcodec_alloc_context3(codec);
	if (!codec_context) {
		throw std::runtime_error("Could not allocate AVCodecContext");
	}

	// Initialize a format context
	format_context = avformat_alloc_context();
	if (!format_context) {
		throw std::runtime_error("Could not allocate AVFormatContext");
	}

	// Set video dimensions. They match with the glut output window.
	codec_context->width  = output_width;
	codec_context->height = output_height;
	codec_context->time_base.num = 1;
	codec_context->time_base.den = framerate;
    codec_context->gop_size = 10;
    codec_context->max_b_frames = 1;
    codec_context->pix_fmt = AV_PIX_FMT_YUV420P;
	// IMPORTANT : Not yet sure if context needs more params initialized!

	frame = av_frame_alloc();
	if (!frame) {
		throw std::runtime_error("Could not allocate AVFrame");
	}

	frame->format = codec_context->pix_fmt;
	frame->width  = output_width;
	frame->height = output_height;
}

void Video::video_finalize() {
	// Currently a placeholder
	// TODO : add stuff that save the video file
	AVPacket *pkt = NULL;

	// Send NULL to the context to put it in flush mode.
	avcodec_send_frame(codec_context, NULL);

	// The codec keeps a buffer of packets, so we have to drain it after sending the flush signal.
	// Keep requesting packets until codec returns EOF
	while (avcodec_receive_packet(codec_context, pkt) != AVERROR_EOF) {

	}

	// It's important to destroy the context because the destructor checks if it exists
	// and posts a warning to catch situations where a Video object was destroyed without
	// finalizing the video capture process.
	avcodec_free_context(&codec_context);
	avformat_free_context(format_context);
	format_context = NULL;
}