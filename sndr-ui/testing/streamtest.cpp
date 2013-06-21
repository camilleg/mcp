/*
    Synthesize a video signal, and stream it over SDP
    ( http://en.wikipedia.org/wiki/Session_Description_Protocol ).

    Usage:
        ./streamtest &

    cat > myStream.sdp
    v=0
    o=- 0 0 IN IP4 127.0.0.1
    s=No Name
    c=IN IP4 127.0.0.1
    t=0 0
    a=tool:libavformat 55.2.100
    m=video 49990 RTP/AVP 96
    a=rtpmap:96 H264/90000
    a=fmtp:96 packetization-mode=1
    ^D
    # The values 127.0.0.1 and 49990 must agree with argv[1] and argv[2].

    To play what streamtest streams:
	ffplay -i myStream.sdp
	vlc myStream.sdp	# causes errors
*/

#include <cstdio>
#include <stdint.h>
#include <unistd.h>
#include <cassert>

// Workaround for libavutil/common.h:170:47: error: ‘UINT64_C’ was not declared in this scope
#ifndef INT64_C
#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)
#endif
extern "C" {
#include <x264.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
}

const int WIDTH = 320;
const int HEIGHT = 240;
const int FPS = 30;
const int BITRATE = 400000;


struct AVFormatContext* avctx;
struct x264_t* encoder;
struct SwsContext* imgctx;
AVCodecContext* pCodecCtx = NULL;

void create_sample_picture(x264_picture_t* picDst)
{
  // Generate a test signal, in imgctx's AVPixelFormat AV_PIX_FMT_RGB24,
  // described in libavutil/pixfmt.h and http://libav.org/doxygen/master/pixfmt_8h.html .
  static double t = 0.0;
  t += 0.3;
  static uint8_t rgbStore[WIDTH * HEIGHT * 3];
  for (int y=0;y<HEIGHT; ++y) {
    for (int x=0; x<WIDTH; ++x) {
      uint8_t* rgb = rgbStore + ((y*WIDTH + x) *3);
      const double i = x/double(WIDTH);
      const double j = y/double(HEIGHT);
      rgb[0] = (y<HEIGHT/4) ? 255*i * drand48() : 0;
      rgb[1] = (x<WIDTH/2) ? 90 : 0;
      rgb[2] = 255*j*(0.5+0.5*sin(t));
    }
  }
  AVFrame* picSrc = avcodec_alloc_frame();
  avpicture_fill((AVPicture *)picSrc, rgbStore, AV_PIX_FMT_RGB24, WIDTH, HEIGHT);

#ifdef flip_frame_vertically
  picSrc->data[0] += picSrc->linesize[0] * (HEIGHT-1);
  picSrc->linesize[0] *= -1;
#endif

  // Create the frame *picDst for sws_scale() to stuff.
  x264_picture_alloc(picDst, X264_CSP_I420, WIDTH, HEIGHT);

  // Scale the image.
  // http://ffmpeg.org/doxygen/trunk/group__lsws.html#gae531c9754c9205d90ad6800015046d74
  const int hSlice = sws_scale(imgctx,
    picSrc->data, picSrc->linesize, 0, HEIGHT,
    picDst->img.plane, picDst->img.i_stride);
  assert(hSlice == HEIGHT);
  av_free(picSrc);
}

int encode_frame(x264_picture_t* picSrc, x264_nal_t** nals)
{
  // Encode a frame into a sequence of NAL units.

#ifdef doesnt_work
  // Even when commenting out where AVPacket sets p.pts = AV_NOPTS_VALUE,
  // this didn't remove the warning:
  //     [rtp @ 0x269c0a0] Encoder did not produce proper pts, making some up.
  static int pts = 0;
  ((AVFrame*)picSrc)->pts = ++pts; // c->frame_number; //;;;; (1.0/FPS) * ...;
#endif

  x264_picture_t picDst;
  int num_nals;
  const int frame_size = x264_encoder_encode(encoder, nals, &num_nals, picSrc, &picDst);
  if (frame_size >= 0 && num_nals < 0)
    printf("invalid frame size: %d\n", num_nals);

  // Ignore bad frames (huh?).
  if (frame_size < 0)
    return frame_size;

  return num_nals;
}

// Initalize and send a packet.
void stream_frame(uint8_t* payload, int size)
{
  static AVPacket p;
  av_init_packet(&p);
  p.data = payload;
  p.size = size;
  p.stream_index = 0;
  p.flags = AV_PKT_FLAG_KEY;
  p.pts = AV_NOPTS_VALUE;
  p.dts = AV_NOPTS_VALUE;
  av_interleaved_write_frame(avctx, &p);
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
      printf("usage: streamtest dottedquad port\n");
      return -1;
    }

    // Initalize ffmpeg.
    av_register_all();
    avformat_network_init();

    // Initialize image scaler
    imgctx = sws_getContext(WIDTH, HEIGHT, AV_PIX_FMT_RGB24,       // in
                            WIDTH, HEIGHT, AV_PIX_FMT_YUV420P,     // out
                            SWS_FAST_BILINEAR, NULL, NULL, NULL);  // scaling method

    // Initalize encoder.
    x264_param_t param;
    x264_param_default_preset(&param, "ultrafast", "zerolatency");
    param.i_threads = 3;
    param.i_width = WIDTH;
    param.i_height = HEIGHT;
    param.i_fps_num = FPS;
    param.i_fps_den = 1;
    param.i_keyint_max = FPS;
    param.b_intra_refresh = 0;
    param.rc.i_rc_method = X264_RC_ABR; // enable param.rc.i_bitrate
    param.rc.i_bitrate = BITRATE;
    param.b_repeat_headers = 1; // repeat headers or write just once
    param.b_annexb = 1;         // place start codes (1) or sizes (0)
    x264_param_apply_profile(&param, "high");
    encoder = x264_encoder_open(&param);

    // x264_encoder_headers can now be used, but it has had no effect

    // Initialize streaming context.
    // Todo: add standard error handling.
    avctx = avformat_alloc_context();
    struct AVOutputFormat* fmt = av_guess_format("rtp", NULL, NULL);
    avctx->oformat = fmt;

    const int rtpPort = strtol(argv[2], NULL, 10);
    if (rtpPort <= 0 || rtpPort > 65535) {
      printf("streamtest warning: unexpected port number '%s'\n", argv[2]);
    }
    snprintf(avctx->filename, sizeof(avctx->filename), "rtp://%s:%s", argv[1], argv[2] /* an int, e.g. 49990 */ );
    if (avio_open(&avctx->pb, avctx->filename, AVIO_FLAG_WRITE) < 0)
    {
        perror("avio_open failed");
        return 1;
    }
    struct AVStream* stream = avformat_new_stream(avctx, NULL /* pCodecCtx->codec isn't available yet */);

    // Initalize codec.
    pCodecCtx = stream->codec;
    pCodecCtx->codec_id = CODEC_ID_H264;
    pCodecCtx->codec_type = AVMEDIA_TYPE_VIDEO;
    pCodecCtx->flags = CODEC_FLAG_GLOBAL_HEADER;
    pCodecCtx->width = WIDTH;
    pCodecCtx->height = HEIGHT;
    pCodecCtx->time_base.den = FPS;
    pCodecCtx->time_base.num = 1;
    pCodecCtx->gop_size = FPS;
    pCodecCtx->bit_rate = BITRATE;
    // avctx->flags = AVFMT_FLAG_RTP_HINT;
    // Use avctx->extradata?

    // Write header.
    avformat_write_header(avctx, NULL);

    printf("Emitting frames.\n");
    for (long frame=0;; frame++)
    {
        // Create frame.
        x264_picture_t pic;
        create_sample_picture(&pic);

        // Encode frame into a sequence of Network Abstraction Layer (NAL) units.
        x264_nal_t* nals;
        const int num_nals = encode_frame(&pic, &nals);
        if (num_nals < 0) {
            printf("invalid frame size: %d\n", num_nals);
	    continue;
	}

        // Emit the NALs.
        for (int i = 0; i < num_nals; i++)
        {
            stream_frame(nals[i].p_payload, nals[i].i_payload);
        }

        // Free frame.
        x264_picture_clean(&pic);

        usleep(2000); // imprecise CPU throttle
    }
    return 0;
}
