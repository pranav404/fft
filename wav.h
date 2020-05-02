#pragma once
struct wavheader {
	unsigned char riff[4];//1-4
	unsigned int size;
	unsigned char wave[4];
	unsigned char fmt[4];
	unsigned int length;
	unsigned int format;
	unsigned int no_of_channels;
	unsigned int sample_rate;
	unsigned int byte_rate;
	unsigned int align;
	unsigned int bits_per_sample;
	unsigned char data_chunk_marker[4];
	unsigned int size_of_data;


};