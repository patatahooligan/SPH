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
const char filename[] = "output.mp4";


int Video::save_packets() {
	// Grab packets from codec output and write them to the file. All contexts are expected to
	// have been initialized.

	// TODO : figure out if the packet timebase needs to be adjusted before muxing with
	// void av_packet_rescale_ts(AVPacket * pkt, AVRational tb_src, AVRational tb_dst)

	// Attempt to grab first packet.
	int err_code = avcodec_receive_packet(codec_context, pkt);

	while (err_code == 0) {			// 0 means a packet was received

		// Send packet to the format context. Give a warning if an error occurs.
		int write_frame_err_code = av_interleaved_write_frame(format_context, pkt);
		if (write_frame_err_code != 0) {
			std::cerr << "av_interleaved_write_frame returned " << write_frame_err_code << std::endl;
		}
		err_code = avcodec_receive_packet(codec_context, pkt);
	}

	// These error codes are expected when the codec needs more frames to encode or has
	// been flushed. Post warning if a different error code occurs.
	if (err_code != AVERROR(EAGAIN) && err_code != AVERROR_EOF) {
		std::cerr << "avcodec_receive_packet returned " << err_code << std::endl;
	}

	// Return the code because it might be useful for the caller to know if the codec was in flush mode.
	return err_code;
}


Video::~Video() {
	// Free all AVCodec structs that have been allocated while initializing.
	// Normally all contexts should be deallocated in video_finalize(), so this code serves to warn of
	// possible bugs and premature Video destruction. It also doesn't allow memory leaks if a Video object
	// is mishandled.

	if (codec_context) {
		avcodec_free_context(&codec_context);
		std::cerr << "WARNING : A Video object was destroyed with a live AVCodecContext. Possibly a video file was not finalized.";
	}

	if (format_context) {
		avformat_free_context(format_context);
		std::cerr << "WARNING : A Video object was destroyed with a live AVFormatContext. Possibly a video file was not finalized.";
	}

	av_frame_free(&frame);

	if (io_context) {
		std::cerr << "WARNING : A Video object was destroyed with a live AVIOContext. Possibly a video file was not finalized.";
		int err_code = avio_close(io_context);
		if (err_code != 0) {
			std::cerr << "avio_close returned " << err_code << std::endl;
		}
	}
}

void Video::video_init() {
	// Initialize libav

	// Initializes libavformat and registers all the muxers, demuxers and protocols.
	av_register_all();

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

	int err_code = avio_open(&io_context, filename, AVIO_FLAG_WRITE);
	if (err_code < 0) {
		std::cerr << "avio_open returned " << err_code << std::endl;
		throw std::runtime_error("Could not allocate avio_open");
	}

	// Set codec parameters. Dimensions match with the glut output window.
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

	// Send NULL to the context to put it in flush mode.
	avcodec_send_frame(codec_context, NULL);

	int err_code = save_packets();
	if (err_code != AVERROR_EOF) {
		throw std::runtime_error("Unexpected error in save_packets");
	}

	// Free all the stuff. The destructor also does this, but only as a fail-safe. 
	avcodec_free_context(&codec_context);
	avformat_free_context(format_context);
	format_context = NULL;
	err_code = avio_close(io_context);
	io_context = NULL;
	if (err_code != 0) {
		std::cerr << "avio_close returned " << err_code << std::endl;
	}
}