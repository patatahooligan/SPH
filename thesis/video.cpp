#include "stdafx.h"

#include <stdexcept>
#include <iostream>
#include <assert.h>

// libav is a pure C project so it needs to be included as such
extern "C" {
	#include "libavcodec\avcodec.h"
	#include "libavformat\avformat.h"
	#include "libavformat\avio.h"
	#include "libavutil\imgutils.h"
	#include "libswscale\swscale.h"
	#include "libavutil\opt.h"
}

#include "video.h"
#include "constants.h"


const AVCodecID codec_id = AV_CODEC_ID_H264;
const int framerate = 25;
const char filename[] = "C:\DevStuff\thesis\output.mp4";


int Video::save_packets() {
	// Grab packets from codec output and write them to the file. All contexts are expected to
	// have been initialized.

	// TODO : figure out if the packet timebase needs to be adjusted before muxing with
	// void av_packet_rescale_ts(AVPacket * pkt, AVRational tb_src, AVRational tb_dst)

	// Attempt to grab first packet.
	AVPacket pkt;
	av_init_packet(&pkt);
	pkt.data = NULL;
	pkt.size = 0;
	int err_code = avcodec_receive_packet(codec_context, &pkt);

	while (err_code == 0) {			// 0 means a packet was received

		// Send packet to the format context. Give a warning if an error occurs.
		int write_frame_err_code = av_interleaved_write_frame(format_context, &pkt);
		if (write_frame_err_code != 0) {
			char error[AV_ERROR_MAX_STRING_SIZE];
			std::cerr << "av_interleaved_write_frame returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, write_frame_err_code) << std::endl;
			return err_code;
		}
		err_code = avcodec_receive_packet(codec_context, &pkt);
	}

	// These error codes are expected when the codec needs more frames to encode or has
	// been flushed. Post warning if a different error code occurs.
	if (err_code != AVERROR(EAGAIN) && err_code != AVERROR_EOF) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cerr << "avcodec_receive_packet returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}

	// Return the code because it might be useful for the caller to know if the codec was in flush mode.
	return err_code;
}


Video::~Video() {
	// Call video_finalize(). This should have been done manually, but in case of crashes or user error, this might save the
	// video file after all. Post a bunch of warnings if things have not been finalized to draw attention to bugs.

	if (!finalized) {
		std::cerr << "WARNING : A Video object is not in a finalized state. The program has either crashed or video.finalize() was not called " << std::endl;
	}

	if (codec_context) {
		std::cerr << "WARNING : A Video object was destroyed with a live AVCodecContext. Possibly a video file was not finalized.";
	}

	if (format_context) {
		std::cerr << "WARNING : A Video object was destroyed with a live AVFormatContext. Possibly a video file was not finalized.";
	}

	if (format_context->pb) {
		std::cerr << "WARNING : A Video object was destroyed with a live AVIOContext. Possibly a video file was not finalized.";
	}

	video_finalize();
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

	// Set codec parameters. Dimensions match with the glut output window.
	codec_context->width  = output_width;
	codec_context->height = output_height;
	codec_context->time_base.num = 1;
	codec_context->time_base.den = framerate;
    codec_context->gop_size = 10;
    codec_context->max_b_frames = 2;
    codec_context->pix_fmt = AV_PIX_FMT_YUV420P;
	codec_context->bit_rate = 2500000;

	// Initialize a format context
	int err_code = avformat_alloc_output_context2(&format_context, NULL, NULL, filename);
	if (!format_context) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cerr << "avformat_alloc_output_context2 returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
		throw std::runtime_error("Could not allocate AVFormatContext");
	}

	ostream = avformat_new_stream(format_context, codec);

	err_code = avformat_write_header(format_context, NULL);
	if (err_code < 0) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cerr << "avformat_write_header returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}

	// Copy codec parameters to output stream
	avcodec_parameters_from_context(ostream->codecpar, codec_context);

	err_code = av_opt_set(codec_context->priv_data, "preset", "slow", 0);
	if (err_code < 0) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cout << "av_opt_set returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}
	err_code = av_opt_set(codec_context->priv_data, "profile", "baseline", AV_OPT_SEARCH_CHILDREN);
	if (err_code < 0) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cout << "av_opt_set returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}

	err_code = avcodec_open2(codec_context, codec, NULL);
	if (err_code) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cout << "av_opt_set returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
		throw std::runtime_error("Unexpected error in video_init when calling avcodec_open2");
	}

	// Allocate frame objects
	rgbframe = av_frame_alloc();
	yuvframe = av_frame_alloc();
	if (!rgbframe) {
		throw std::runtime_error("Could not allocate AVFrame");
	}

	// Set their parameters
	rgbframe->format = AV_PIX_FMT_RGB24;
	rgbframe->width  = output_width;
	rgbframe->height = output_height;
	rgbframe->format = codec_context->pix_fmt;
	rgbframe->width  = output_width;
	rgbframe->height = output_height;

	yuvframe->format = AV_PIX_FMT_YUV420P;
	yuvframe->width  = output_width;
	yuvframe->height = output_height;
	yuvframe->format = codec_context->pix_fmt;
	yuvframe->width  = output_width;
	yuvframe->height = output_height;

	// Allocate an SwScontext which will be used to convert RGB24->YUV420P
	sws = sws_getContext (
		output_height, output_width, AV_PIX_FMT_RGB24,
		output_height, output_width, AV_PIX_FMT_YUV420P,
		0, NULL, NULL, NULL);
	if (!sws) {
		throw std::runtime_error("Could not allocate SwsContext");
	}

	// Allocate space for the image data within the frame objects
	err_code = av_image_alloc(rgbframe->data, rgbframe->linesize, output_width, output_height, codec_context->pix_fmt, 32);
	if (err_code < 0) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cout << "av_image_alloc returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}
	err_code = av_image_alloc(yuvframe->data, yuvframe->linesize, output_width, output_height, codec_context->pix_fmt, 32);
	if (err_code < 0) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cout << "av_image_alloc returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}
}

void Video::encode_frame(float simulation_time) {
	if (!sws)

	// Check if enough time has passed to warrant a frame.
	if (simulation_time <= current_frame * framerate) return;

	// Grab color info from the render buffer
	GLubyte* gl_image = new GLubyte[output_width*output_height*3];
	glReadPixels(0, 0, output_width, output_height, GL_RGB, GL_UNSIGNED_BYTE, gl_image);

	// Move the data to the frame struct. The reason we don't just pass frame->data to glReadPixels is
	// that there may be padding around the image data meaning that possibly linesize!= ouput_width so
	// glReadPixels would not move the data to the proper location.
	for (unsigned y = 0; y < output_height; y++) {
		for (unsigned x = 0; x < output_width; x++) {
			for (unsigned color = 0; color < 3; color++) {
				rgbframe->data[0][y * rgbframe->linesize[0] + 3 * x + color] = gl_image[y * 3 * output_width + 3 * x + color];
			}
		}
	}

	delete [] gl_image;
	gl_image = NULL;

	// Convert data from rgbframe to yuvframe
	int err_code = sws_scale(sws, rgbframe->data, rgbframe->linesize, 0, output_height, yuvframe->data, yuvframe->linesize);
	assert(err_code == output_height);

	// In case the simulation step is larger than the fps (unlikely), send the frame multiple times.
	do {
		yuvframe->pts = current_frame;
		int errcode = avcodec_send_frame(codec_context, yuvframe);
		if (errcode) {
			char error[AV_ERROR_MAX_STRING_SIZE];
			std::cerr << "avcodec_send_frame returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, errcode) << std::endl;
		}
		current_frame++;
	} while (simulation_time > current_frame*framerate);

	err_code = save_packets();

	if (err_code != AVERROR(EAGAIN)) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cerr << "save_packets returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
		throw std::runtime_error("Unexpected error in save_packets");
	}
}

void Video::video_finalize() {
	// Currently a placeholder
	// TODO : add stuff that save the video file

	// Fail safe if already finalized
	if (finalized) return;
	
	// Send NULL to the context to put it in flush mode.
	avcodec_send_frame(codec_context, NULL);

	int err_code = save_packets();
	if (err_code != AVERROR_EOF) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cerr << "save_packets returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}

	// Free all the stuff.
	
	avcodec_free_context(&codec_context);

	err_code = avio_close(format_context->pb);
	format_context->pb = NULL;
	if (err_code != 0) {
		char error[AV_ERROR_MAX_STRING_SIZE];
		std::cerr << "avio_close returned " << av_make_error_string(error, AV_ERROR_MAX_STRING_SIZE, err_code) << std::endl;
	}

	avformat_free_context(format_context);
	format_context = NULL;
	finalized = true;
}