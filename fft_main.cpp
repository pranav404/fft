#pragma warning(disable:4996)
#include<iostream>
#include<conio.h>
#include<cstdlib>
#include<string>
#include "wav.h"
#include<math.h>
using namespace std;

void errorDisplay(const char* msg) {
	cout << msg << endl;

}
void fft(float* ipr, float tr[], float ti[], double* ptr)
{
	float br[128], bi[128];
	float r, r1, im, im1;
	int i, k, j, bits;
	for (i = 0; i < 128; i++)               //Performing bit reversal
	{
		k = i;
		for (bits = j = 0; bits < 7; bits++)
		{
			j = (j << 1) | (k & 1);
			k >>= 1; 

		}
		br[j] = *(ipr + i);
		bi[j] = 0.0;

	}
	//First stage
	for (k = 0; k < 127; k += 2)
	{
		r = br[k];
		im = bi[k];
		r1 = br[k + 1];
		im1 = bi[k + 1];
		br[k] = (r + r1);
		bi[k] = (im + im1);
		br[k + 1] = (r - r1);
		bi[k + 1] = (im - im1);

	}

	//Second stage
	for (i = 0; i < 127; i += 4)
	{
		for (k = 0; k < 2; k++)
		{
			r = br[k + i];
			im = bi[k + i];
			r1 = br[k + i + 2];
			im1 = bi[k + i + 2];
			br[k + i] = (r + ((tr[32 * k] * r1) - (ti[32 * k] * im1)));
			bi[k + i] = (im + ((tr[32 * k] * im1) + (ti[32 * k] * r1)));
			br[k + i + 2] = (r - ((tr[32 * k] * r1) - (ti[32 * k] * im1)));
			bi[k + i + 2] = (im - ((tr[32 * k] * im1) + (ti[32 * k] * r1)));
		}
	}


	//third stage
	for (i = 0; i < 127; i += 8)
	{
		for (k = 0; k < 4; k++)
		{
			r = br[k + i];
			im = bi[k + i];
			r1 = br[k + i + 4];
			im1 = bi[k + i + 4];
			br[k + i] = (r + ((tr[16 * k] * r1) - (ti[16 * k] * im1)));
			bi[k + i] = (im + ((tr[16 * k] * im1) + (ti[16 * k] * r1)));
			br[k + i + 4] = (r - ((tr[16 * k] * r1) - (ti[16 * k] * im1)));
			bi[k + i + 4] = (im - ((tr[16 * k] * im1) + (ti[16 * k] * r1)));
		}
	}
	//fourth stage:
	for (i = 0; i < 127; i += 16)
	{
		for (k = 0; k < 8; k++)
		{
			r = br[k + i];
			im = bi[k + i];
			r1 = br[k + i + 8];
			im1 = bi[k + i + 8];
			br[k + i] = (r + ((tr[8 * k] * r1) - (ti[8 * k] * im1)));
			bi[k + i] = (im + ((tr[8 * k] * im1) + (ti[8 * k] * r1)));
			br[k + i + 8] = (r - ((tr[8 * k] * r1) - (ti[8 * k] * im1)));
			bi[k + i + 8] = (im - ((tr[8 * k] * im1) + (ti[8 * k] * r1)));
		}
	}
	//Fifth Stage
	for (i = 0; i < 127; i += 32)
	{
		for (k = 0; k < 16; k++)
		{
			r = br[k + i];
			im = bi[k + i];
			r1 = br[k + i + 16];
			im1 = bi[k + i + 16];
			br[k + i] = (r + ((tr[4 * k] * r1) - (ti[4 * k] * im1)));
			bi[k + i] = (im + ((tr[4 * k] * im1) + (ti[4 * k] * r1)));
			br[k + i + 16] = (r - ((tr[4 * k] * r1) - (ti[4 * k] * im1)));
			bi[k + i + 16] = (im - ((tr[4 * k] * im1) + (ti[4 * k] * r1)));
		}
	}
	//6th stage
	for (i = 0; i < 127; i += 64)
	{
		for (k = 0; k < 32; k++)
		{
			r = br[k + i];
			im = bi[k + i];
			r1 = br[k + i + 32];
			im1 = bi[k + i + 32];
			br[k + i] = (r + ((tr[2 * k] * r1) - (ti[2 * k] * im1)));
			bi[k + i] = (im + ((tr[2 * k] * im1) + (ti[2 * k] * r1)));
			br[k + i + 32] = (r - ((tr[2 * k] * r1) - (ti[2 * k] * im1)));
			bi[k + i + 32] = (im - ((tr[2 * k] * im1) + (ti[2 * k] * r1)));
		}
	}
	//Seventh Stage
	for (i = 0; i < 127; i += 128)
	{
		for (k = 0; k < 64; k++)
		{
			r = br[k + i];
			im = bi[k + i];
			r1 = br[k + i + 64];
			im1 = bi[k + i + 64];
			br[k + i] = (r + ((tr[k] * r1) - (ti[k] * im1)));
			bi[k + i] = (im + ((tr[k] * im1) + (ti[k] * r1)));
			br[k + i + 64] = (r - ((tr[k] * r1) - (ti[k] * im1)));
			bi[k + i + 64] = (im - ((tr[k] * im1) + (ti[k] * r1)));
		}
	}
	for (i = 0; i < 128; i++)
	{
		*(ptr + i) = (double)sqrt((br[i] * br[i]) + (bi[i] * bi[i]));
	}


}
void fft256(float* ipr, float tr[], float ti[], double* op)
{
	float br[256], bi[256], r, i, r1, i1;
	int k, j, temp, bits;
	//performing bit reversal
	for (k = 0; k < 256; k++)
	{
		temp = k;
		j = 0;
		for (bits = 0; bits < 8; bits++)
		{
			j = (j << 1) | (temp & 1);
			temp >>= 1;
		}
		br[j] = *(ipr + k);
		bi[k] = 0;
	}
	//first stage
	for (k = 0; k < 256; k += 2)
	{
		r = br[k];
		i = bi[k];
		r1 = br[k + 1];
		i1 = bi[k + 1];
		br[k] = r + r1;
		bi[k] = i + i1;
		br[k + 1] = r - r1;
		bi[k + 1] = i - i1;
	}
	//second stage
	for (k = 0; k < 256; k+=4)
	{
		for (j = 0; j < 2; j++)
		{
			r = br[k+j];
			i = bi[k+j];
			r1 = br[k + j + 2];
			i1 = bi[k + j + 2];
			br[k + j] = r + ((r1 * tr[64 * j]) - (i1 * ti[64 * j]));
			bi[k + j] = i + ((r1 * ti[64 * j]) + (i1 * tr[64 * j]));
			br[k + j + 2] = r - ((r1 * tr[64 * j]) - (i1 * ti[64 * j]));
			bi[k + j + 2] = i - ((r1 * ti[64 * j]) + (i1 * tr[64 * j]));

		}
	}
	//third stage
	for (k = 0; k < 256; k += 8)
	{
		for (j = 0; j < 4; j++)
		{
			r = br[k + j];
			i = bi[k + j];
			r1 = br[k + j + 4];
			i1 = bi[k + j + 4];
			br[k + j] = r + ((r1 * tr[32 * j]) - (i1 * ti[32 * j]));
			bi[k + j] = i + ((r1 * ti[32 * j]) + (i1 * tr[32 * j]));
			br[k + j + 4] = r - ((r1 * tr[32 * j]) - (i1 * ti[32 * j]));
			bi[k + j + 4] = i - ((r1 * ti[32 * j]) + (i1 * tr[32 * j]));

		}
	}
	//fourth stage
	for (k = 0; k < 256; k += 16)
	{
		for (j = 0; j < 8; j++)
		{
			r = br[k + j];
			i = bi[k + j];
			r1 = br[k + j + 8];
			i1 = bi[k + j + 8];
			br[k + j] = r + ((r1 * tr[16 * j]) - (i1 * ti[16 * j]));
			bi[k + j] = i + ((r1 * ti[16 * j]) + (i1 * tr[16 * j]));
			br[k + j + 8] = r - ((r1 * tr[16 * j]) - (i1 * ti[16 * j]));
			bi[k + j + 8] = i - ((r1 * ti[16 * j]) + (i1 * tr[16 * j]));

		}
	}
	//fifth stage
	for (k = 0; k < 256; k += 32)
	{
		for (j = 0; j < 16; j++)
		{
			r = br[k + j];
			i = bi[k + j];
			r1 = br[k + j + 16];
			i1 = bi[k + j + 16];
			br[k + j] = r + ((r1 * tr[8 * j]) - (i1 * ti[8 * j]));
			bi[k + j] = i + ((r1 * ti[8 * j]) + (i1 * tr[8 * j]));
			br[k + j + 16] = r - ((r1 * tr[8 * j]) - (i1 * ti[8 * j]));
			bi[k + j + 16] = i - ((r1 * ti[8 * j]) + (i1 * tr[8 * j]));

		}
	}
	//sixth stage
	for (k = 0; k < 256; k += 64)
	{
		for (j = 0; j < 32; j++)
		{
			r = br[k + j];
			i = bi[k + j];
			r1 = br[k + j + 32];
			i1 = bi[k + j + 32];
			br[k + j] = r + ((r1 * tr[4 * j]) - (i1 * ti[4 * j]));
			bi[k + j] = i + ((r1 * ti[4 * j]) + (i1 * tr[4 * j]));
			br[k + j + 32] = r - ((r1 * tr[4 * j]) - (i1 * ti[4 * j]));
			bi[k + j + 32] = i - ((r1 * ti[4 * j]) + (i1 * tr[4 * j]));

		}
	}
	//seventh stage
	for (k = 0; k < 256; k += 128)
	{
		for (j = 0; j < 64; j++)
		{
			r = br[k + j];
			i = bi[k + j];
			r1 = br[k + j + 64];
			i1 = bi[k + j + 64];
			br[k + j] = r + ((r1 * tr[2 * j]) - (i1 * ti[2 * j]));
			bi[k + j] = i + ((r1 * ti[2 * j]) + (i1 * tr[2 * j]));
			br[k + j + 64] = r - ((r1 * tr[2 * j]) - (i1 * ti[2 * j]));
			bi[k + j + 64] = i - ((r1 * ti[2 * j]) + (i1 * tr[2 * j]));
		}
	}
	//eight stage
	for (k = 0; k < 256; k += 256)
	{
		for (j = 0; j < 128; j++)
		{
			r = br[k + j];
			i = bi[k + j];
			r1 = br[k + j + 128];
			i1 = bi[k + j + 128];
			br[k + j] = r + ((r1 * tr[j]) - (i1 * ti[j]));
			bi[k + j] = i + ((r1 * ti[j]) + (i1 * tr[j]));
			br[k + j + 128] = r - ((r1 * tr[j]) - (i1 * ti[j]));
			bi[k + j + 128] = i - ((r1 * ti[j]) + (i1 * tr[j]));

		}
	}
	for (k = 0; k < 256; k++)
	{
		*(op + k) = sqrt((br[k] * br[k]) + (bi[k] * bi[k]));
	}
}
unsigned char buff2[2];
unsigned char buff4[4];
wavheader h;
int main()
{
	FILE* fp1, * fp2, * fp3, * fp4;

	const char c[200] = "as3.wav";//Change the path to your input .wav file
	fp1 = fopen(c, "rb");//Opening the .wav file


	/*Check whether the input directory is valid or not
	if invalid, exit the program*/


	if (fp1 == NULL)
	{
		errorDisplay("No such directory..");
		exit(0);
	}
	int read_count = 0;//Number of Bytes read form the .wav file
	read_count = fread(h.riff, sizeof(h.riff), 1, fp1);

	if (strcmp((const char*)h.riff, "RIFF"))
	{
		errorDisplay("Invalid Format.");
		exit(0);
	}
	//reading the overall file size
	read_count = fread(buff4, sizeof(buff4), 1, fp1);
	h.size = buff4[0] | (buff4[1] << 8) | (buff4[2] << 16) | (buff4[3] << 24);
	//reading the next 4 bytes 9-12
	read_count = fread(h.wave, sizeof(h.wave), 1, fp1);
	if (strcmpi((const char*)h.wave, "WAVE"))
	{
		errorDisplay("Not a WAVE file.");
		exit(0);
	}
	//reading the FMT 13-16
	read_count = fread(h.fmt, sizeof(h.fmt), 1, fp1);
	read_count = fread(buff4, sizeof(buff4), 1, fp1);
	//converting to big endian notation
	h.length = buff4[0] | (buff4[1] << 8) | (buff4[2] << 16) | (buff4[3] << 24);
	//checking whether the wave is PCM or not
	read_count = fread(buff2, sizeof(buff2), 1, fp1);
	//converting to big endian notation
	h.format = buff2[0] | (buff2[1] << 8);
	//reading the number of channels
	read_count = fread(buff2, sizeof(buff2), 1, fp1);
	//converting to big endian notation
	h.no_of_channels = buff2[0] | (buff2[1] << 8);
	//reading the sample rate
	read_count = fread(buff4, sizeof(buff4), 1, fp1);
	//converting to big endian notation
	h.sample_rate = buff4[0] | (buff4[1] << 8) | (buff4[2] << 16) | (buff4[3] << 24);
	//reading the byte rate
	read_count = fread(buff4, sizeof(buff4), 1, fp1);
	//converting to big endian notation
	h.byte_rate = buff4[0] | (buff4[1] << 8) | (buff4[2] << 16) | (buff4[3] << 24);
	//reading the alignment
	read_count = fread(buff2, sizeof(buff2), 1, fp1);
	//converting to big endian notation
	h.align = buff2[0] | (buff2[1] << 8);
	//reading the bits per sample
	read_count = fread(buff2, sizeof(buff2), 1, fp1);
	//converting to big endian notation
	h.bits_per_sample = buff2[0] | (buff2[1] << 8);
	//reading the data chunk marker
	read_count = fread(h.data_chunk_marker, sizeof(h.data_chunk_marker), 1, fp1);
	//converting into big endian notation
	//reading the data size
	read_count = fread(buff4, sizeof(buff4), 1, fp1);
	//converting into big endian notation
	h.size_of_data = buff4[0] | (buff4[1] << 8) | (buff4[2] << 16) | (buff4[3] << 24);
	printf("\n \n");
	printf("1 - 4: %s\n", h.riff);
	printf("5 - 8: Size = %d\n", h.size);
	printf("9 - 12: Format : %s\n", h.wave);
	printf("13 - 16: fmt : %s\n", h.fmt);
	printf("17 - 20: Length of format data : %d\n", h.length);
	switch (h.format)
	{
	case 1: printf("PCM\n");
		break;
	case 6: printf("A-LAW\n");
		break;
	case 7: printf("Mu-law\n");
		break;
	default: {printf("Invalid format\n");
		exit(0);
	}

	}
	//Print the header information of the .wav file
	printf("23 - 24: Number of channels : %d channels\n", h.no_of_channels);
	printf("25 - 28: Sample rate : %d hz\n", h.sample_rate);
	printf("29 - 32: Byte rate : %d\n", h.byte_rate);
	printf("33 - 34: Align number : %d\n", h.align);
	printf("35 - 36: Bit per entry : %d bits\n", h.bits_per_sample);
	printf("37 - 40: Data chunk header : %s\n", h.data_chunk_marker);
	printf("40 - 44: Size of the data section : %d\n", h.size_of_data);
	float duration = (float)h.size_of_data / h.byte_rate;
	printf("Duration of track: %0.3f seconds\n", duration);
	long sample_num = (8 * h.size_of_data) / (h.no_of_channels * h.bits_per_sample);
	printf("Number of samples : %ld samples\n", sample_num);
	long sample_size = (h.no_of_channels * h.bits_per_sample) / 8;
	printf("Size of each sample : %ld bytes\n", sample_size);
	long bytes_per_sample = sample_size / h.no_of_channels;
	if (sample_size != (bytes_per_sample * h.no_of_channels))
	{
		errorDisplay("Error. Invalid value.");
		exit(0);
	}
	long low_limit = 01, high_limit = 01; //Setting the high and low limits for the samples
	switch (h.bits_per_sample)
	{
	case 8:low_limit = -128;
		high_limit = 127;
		break;
	case 16:
		low_limit = -32768;
		high_limit = 32767;
		break;
	case 32:
		low_limit = -2147483648;
		high_limit = 2147483647;
		break;

	}
	printf("Lower and higher limits are %ld and %ld respectively.\n", low_limit, high_limit);
	printf("Proceed for FFT y/n ?\n");
	char flag = 'y';
	scanf("%c", &flag);

	if ((flag == 'n') || (flag == 'N'))
		exit(0);

	int i, j, k;
	unsigned char data_buff[2];//data buffer to read the samples
	int data;
	int pass_num = 0;
	printf("Enter the nummber of points desired\n");
	scanf("%d", &i);
	if (i == 256)
	{
		fp3 = fopen("cosine256.txt", "r");//Opening the twiddle factors stored in a seperate text file
		fp4 = fopen("sinee256.txt", "r");//Opening the twiddle factors stored in a seperate text file
		fp2 = fopen("Output_spectrogram256.txt", "w");//opening the output file
		float datatowrite[256], tr[128], ti[128];//input array, twiddle array
		double output[256];
		pass_num = sample_num / 256; //Calculating the Number of Windows
		printf("%d\n", pass_num);
		k = 0;
		while ((fscanf(fp3, "%f", &tr[k]) != EOF) && (fscanf(fp4, "%f", &ti[k]) != EOF))//loading the twiddle factor values
		{
			k++;
		}
		for (k = 0; k < pass_num; k++)
		{

			for (i = 0; i < 256; i++)
			{

				read_count = fread(data_buff, sizeof(data_buff), 1, fp1);
				data = 0;
				if (read_count == 1)
				{
					data = data_buff[0] | (data_buff[1] << 8);
					datatowrite[i] = (float)data / high_limit;
					if (datatowrite[i] > 1.0)
					{
						datatowrite[i] = datatowrite[i] - 2;
					}
				}

			}
			fprintf(fp2, "%d.\n", k);
			fft256(datatowrite, tr, ti, output);//passing the input to the fft subroutine
			for (j = 0; j < 256; j++)
			{
				fprintf(fp2, "k = %d-%0.5e\n", j + 1, output[j]);
			}
			fprintf(fp2, "\n\n new Pass \n\n");


		}
	}
	else
	{
		fp3 = fopen("cosine128.txt", "r");//Opening the twiddle factors stored in a seperate text file
		fp4 = fopen("sinee128.txt", "r");//Opening the twiddle factors stored in a seperate text file
		fp2 = fopen("Output_spectrogram128.txt", "w");//opening the output file
		float datatowrite[128], tr[64], ti[64];//input array, twiddle array
		double output[128];
		pass_num = sample_num / 128; //Calculating the Number of Windows
		printf("%d\n", pass_num);
		k = 0;
		while ((fscanf(fp3, "%f", &tr[k]) != EOF) && (fscanf(fp4, "%f", &ti[k]) != EOF))//loading the twiddle factor values
		{
			k++;
		}
		for (k = 0; k < pass_num; k++)
		{

			for (i = 0; i < 128; i++)
			{

				read_count = fread(data_buff, sizeof(data_buff), 1, fp1);
				data = 0;
				if (read_count == 1)
				{
					data = data_buff[0] | (data_buff[1] << 8);
					datatowrite[i] = (float)data / high_limit;
					if (datatowrite[i] > 1.0)
					{
						datatowrite[i] = datatowrite[i] - 2;
					}
				}

			}
			fprintf(fp2, "%d.\n", k);
			fft(datatowrite, tr, ti, output);//passing the input to the fft subroutine
			for (j = 0; j < 128; j++)
			{
				fprintf(fp2, "k = %d-%0.5e\n", j + 1, output[j]);
			}
			fprintf(fp2, "\n\n new Pass \n\n");


		}
	}
	printf("Closing file..\n");//Closing all the opened files
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);




	return 0;
}