// Workaround for libavutil/common.h:170:47: error: ‘UINT64_C’ was not declared in this scope
#ifndef INT64_C
#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)
#endif

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
}

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

void emit_frame(const AVFrame *frame)
{
    const int n = frame->nb_samples * av_get_channel_layout_nb_channels(av_frame_get_channel_layout(frame));
    const uint16_t *p     = (uint16_t*)frame->data[0];
    const uint16_t *p_end = p + n;

#if 0
    while (p < p_end) {
        fputc(*p    & 0xff, stdout);
        fputc(*p>>8 & 0xff, stdout);
        p++;
    }
#else
    static int a = 0; a += n;
    int64_t dur = av_frame_get_pkt_duration(frame); // in AVStream->time_base units
    double seconds = (double(dur) * timebase.num) / timebase.den;
    static double s = 0; s += seconds;
    printf("another %d samples, %f seconds, respective totals %d %f\n", n, seconds, a, s);
#endif
    fflush(stdout);
}

int main(int argc, char **argv)
{
    AVFrame *src_frame = av_frame_alloc();
    if (!src_frame) {
        perror("Could not allocate frame");
        return 1;
    }

    avcodec_register_all();
    av_register_all();

    AVPacket packet;
    int ret;
    if ((ret = open_input_file(argv[1])) >= 0)
    {
	/* read all packets */
	while ((ret = av_read_frame(pFormatContext, &packet)) >= 0) {
	  if (packet.stream_index != iAudioStream)
	      continue;
	  avcodec_get_frame_defaults(src_frame);
	  int got_frame = 0;
	  ret = avcodec_decode_audio4(pDecoderContext, src_frame, &got_frame, &packet);
	  if (ret < 0) {
	      av_log(NULL, AV_LOG_ERROR, "Error decoding audio\n");
	      continue;
	  }
	  if (got_frame) {
	      emit_frame(src_frame);
	  }
	  av_free_packet(&packet);
      }
    }
    if (pDecoderContext)
        avcodec_close(pDecoderContext);
    avformat_close_input(&pFormatContext);
    // av_frame_free(&src_frame); // This causes a segfault "double free or corruption (!prev)."
    if (ret < 0 && ret != AVERROR_EOF) {
        char buf[1024];
        av_strerror(ret, buf, sizeof(buf));
        fprintf(stderr, "Error occurred: %s\n", buf);
        return 1;
    }
    return 0;
}
