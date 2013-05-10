// Workaround for libavutil/common.h:170:47: error: ‘UINT64_C’ was not declared in this scope
#ifndef INT64_C
#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)
#endif
extern "C" {
#include <libavcodec/avcodec.h>		// apt-get install libavcodec-dev
#include <libavformat/avformat.h>	// apt-get install libavformat-dev
}

#include <cassert>
#include "array.h"

// todo: encapsulate into a class all globals and statics in this .cpp file
AVFormatContext *pFormatContext = NULL;
AVCodecContext *pDecoderContext = NULL;
int iAudioStream = -1;
AVRational timebase;

int open_input_file(const char *filename)
{
    int ret;
    if ((ret = avformat_open_input(&pFormatContext, filename, NULL, NULL)) < 0) {
        av_log(NULL, AV_LOG_ERROR, "Cannot open input file\n");
        return ret;
    }
    if ((ret = avformat_find_stream_info(pFormatContext, NULL)) < 0) {
        av_log(NULL, AV_LOG_ERROR, "Cannot find stream information\n");
        return ret;
    }

    AVCodec *dec = NULL;
    ret = av_find_best_stream(pFormatContext, AVMEDIA_TYPE_AUDIO, -1, -1, &dec, 0);
    if (ret < 0) {
        av_log(NULL, AV_LOG_ERROR, "Cannot find a audio stream in the input file\n");
        return ret;
    }
    iAudioStream = ret;
    pDecoderContext = pFormatContext->streams[iAudioStream]->codec;
    timebase        = pFormatContext->streams[iAudioStream]->time_base;

    if ((ret = avcodec_open2(pDecoderContext, dec, NULL)) < 0) {
        av_log(NULL, AV_LOG_ERROR, "Cannot open audio decoder\n");
        return ret;
    }
    return 0;
}

AVFrame *src_frame = NULL;

int ffmpeg_init()
{
    src_frame = av_frame_alloc();
    if (!src_frame) {
        perror("Could not allocate frame");
        return 1;
    }
    avcodec_register_all();
    av_register_all();
    return 0;
}

int ret_opened = 0;
bool hit_eof = false;

bool ffmpeg_open(const std::string& filename)
{
  hit_eof = false;
  return (ret_opened = open_input_file(filename.c_str())) >= 0;
}

bool ffmpeg_eof() { return hit_eof; }

double ffmpeg_samplerate() { return pDecoderContext->sample_rate; }

void debug_frame(const AVFrame *frame)
{
  const int n = frame->nb_samples * av_get_channel_layout_nb_channels(av_frame_get_channel_layout(frame));
  static int a = 0; a += n;
  int64_t dur = av_frame_get_pkt_duration(frame); // in AVStream->time_base units
  double seconds = (dur * timebase.num) / double(timebase.den);
  static double s = 0; s += seconds;
  printf("another %d samples, %f seconds, respective totals %d %f\n", n, seconds, a, s);
  const int m = pDecoderContext->frame_size;
  printf("%d samples per channel per audio frame, at %f SR, == %f seconds.\n", m, ffmpeg_samplerate(), m/ffmpeg_samplerate());
}

// ffmpeg can report the number of frames only by explicitly counting them.
// To avoid a separate pass through a file (which doesn't generalize to streamed input),
// it can estimate it: http://superuser.com/questions/84631/how-do-i-get-the-number-of-frames-in-a-video-on-the-linux-command-line

bool ffmpeg_read_frame(uint16_t*& pw, size_t& cw)
{
  if (hit_eof || ret_opened < 0)
      return false;
  int ret;
  AVPacket packet;
again:
  do {
    ret = av_read_frame(pFormatContext, &packet);
    if (ret < 0) {
      hit_eof = true;
      if (pDecoderContext)
	avcodec_close(pDecoderContext);
      avformat_close_input(&pFormatContext);
      if (ret < 0 && ret != AVERROR_EOF) {
	char buf[1024];
	av_strerror(ret, buf, sizeof(buf));
	fprintf(stderr, "Error occurred: %s\n", buf);
      }
      return false;
    }
  } while (packet.stream_index != iAudioStream);

  avcodec_get_frame_defaults(src_frame);
  int got_frame = 0;
  ret = avcodec_decode_audio4(pDecoderContext, src_frame, &got_frame, &packet);
  if (ret < 0) {
    av_log(NULL, AV_LOG_ERROR, "Error decoding audio\n");
    return false;
  }

  if (got_frame) {
    debug_frame(src_frame);
    const size_t csamples = src_frame->nb_samples *
      av_get_channel_layout_nb_channels(av_frame_get_channel_layout(src_frame));
    if (csamples == 0) {
      printf("Skipping frame with 0 samples.\n");;;;
      goto again;
    }
    cw = csamples;
    static array<uint16_t> buf;
    buf.grow(cw);
    const uint16_t* pSrc = (const uint16_t*)src_frame->data[0];
    std::copy(pSrc, pSrc+cw, buf+0); //;;;; extra array-copy.  Nicer would be to av_free_packet "after return."
    pw = buf+0;
  }
  av_free_packet(&packet);
  return true;
}

// for .h file: template <class T> bool ffmpeg_read( array<T> &p)
bool ffmpeg_read(array<double> &dst)
{
  if (hit_eof || ret_opened < 0)
    return false;

  static array<double> cache;
  static size_t icache = 0;

  // Drain cache into dst before reading more frames.
  const size_t cwDst = dst.size();
  size_t i = std::min(cwDst, icache);
  std::copy(cache+0, cache+i, (double*)dst);

  if (i < icache) {
    // Reconstruct incompletely drained cache.
    // todo: for speed, replace memmove'd array with fifo or circular buffer.
    memmove(cache, cache+i, (icache-i) * sizeof(cache[0]));
    icache -= i;
  }

  // Continue filling dst[i..cwDst].
  uint16_t* pw = NULL;
  size_t cw = 0;
  while (i < cwDst) {
    if (!ffmpeg_read_frame(pw, cw)) {
      // Failed to read another frame from file.
      return false;
    }
    // Fill dst[i..] from pw[0..cw].
    const size_t iw = std::min(cwDst-i, cw);
    for (size_t j=0; j<iw; ++j)
      dst[i+j] = double(pw[iw]) / SHRT_MAX;
    i += iw;
    // We can't std::copy(pw, pw+iw, (double*)dst), because we also convert from short to double.
    if (iw < cw) {
      assert(i == cwDst); // Filled dst.
      // Save the unused samples pw[iw..cw] into local cache.
      const size_t icacheNew = icache + cw-iw;
      cache.grow(icacheNew);
      std::copy(pw+iw, pw+cw, cache+icache);
      icache = icacheNew;
    }
  }
  return true;
}
