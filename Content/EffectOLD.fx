//#if OPENGL
//	#define SV_POSITION POSITION
//	#define VS_SHADERMODEL vs_3_0
//	#define PS_SHADERMODEL ps_3_0
//#else
#define VS_SHADERMODEL vs_4_0//_level_9_1
#define PS_SHADERMODEL ps_4_0//_level_9_1
//#endif

matrix WorldViewProjection;

float volume[3];
float low;
float currentVolume;
int sampleRate;
int volumeDelay = 1000;

Texture2D tex;
sampler2D texture_sampler = sampler_state
{
	Texture = <tex>;
};

struct VertexShaderInput
{
	float4 Position : POSITION0;
	float4 Color : COLOR0;
	float2 tex : TEXCOORD0;
};
struct VertexShaderOutput
{
	float4 Position : SV_POSITION;
	float4 Color : COLOR0;
	float2 tex : TEXCOORD0;
};
VertexShaderOutput MainVS(in VertexShaderInput input)
{
	VertexShaderOutput output = (VertexShaderOutput)0;

	output.Position = mul(input.Position, WorldViewProjection);
	output.Color = input.Color;
	output.tex = input.tex;

	return output;
}
float4 MainPS(VertexShaderOutput input) : COLOR
{
	float4 color = tex2D(texture_sampler, input.tex);
	/*color.r *= volume[0];
	color.g *= volume[1];
	color.b *= volume[2];*/
	//color.a *= (volume[0] + volume[1] + volume[2])/3;
	return color;
}

VertexShaderOutput AudioVS(in VertexShaderInput input)
{
	VertexShaderOutput output = (VertexShaderOutput)0;

	output.Position = mul(input.Position, WorldViewProjection);
	output.Color = input.Color;
	output.tex = input.tex;

	return output;
}

/*
float4 sample = tex2D(texture_sampler, input.tex);
	//return sample;
	uint leftChannel =  uint(sample.r * 255.0) | (uint(sample.g * 255.0) << 8);
	uint rightChannel = uint(sample.b * 255.0) | (uint(sample.a * 255.0) << 8);

	float v = min(max(volume, 0.0), 1.0);

	//leftChannel  *= v; //  uint(float(leftChannel) * v);
	//rightChannel *= v; // uint(float(rightChannel) * v);

	sample = float4((leftChannel & 0xff), (leftChannel & 0xff00) >> 8,(rightChannel & 0xff), (rightChannel & 0xff00) >> 8) / 255.0;

	//color *= (input.tex.x);
	return sample;
*/
/*	float4 AudioPS(VertexShaderOutput input) : COLOR
	{


		float4 sample = tex2D(texture_sampler, input.tex);
		//return sample;
		int inLChannel = (int(ceil(sample.r * 255.0))) | ((int(ceil(sample.g * 255.0))) << 8);
		int inRChannel = (int(ceil(sample.b * 255.0))) | ((int(ceil(sample.a * 255.0))) << 8);

		//leftChannel  = leftChannel  - (leftChannel  >> 16) * 256;
		//rightChannel = rightChannel - (rightChannel >> 16) * 256;

		float v = min(max(volume, 0.0), 1.0);

		//int L = inLeftChannel  - 10 * (((inLeftChannel  >> 15) * -2) + 1);
		//int R = inRightChannel - 10 * (((inRightChannel >> 15) * -2) + 1);

	//	int L = ((inLeftChannel  - ((0xffff + 1) * (((inLeftChannel   >> 15) * -2) + 1))) & 0x0ffff) / 2;
	//	int R = ((inRightChannel - ((0xffff + 1) * (((inRightChannel  >> 15) * -2) + 1))) & 0x0ffff) / 2;

	//	int outLeftChannel  = (L + ((0xffff + 1) * (((inLeftChannel   >> 15) * -2) + 1))) & 0x0ffff;
	//	int outRightChannel = (R + ((0xffff + 1) * (((inLeftChannel   >> 15) * -2) + 1)))& 0x0ffff;

		int L = ((~(inLChannel - 1)) & 0xffff) * ((inLChannel & 0xffff) >> 15) + inLChannel * ((~inLChannel & 0xffff) >> 15);
		int R = ((~(inRChannel - 1)) & 0xffff) * ((inRChannel & 0xffff) >> 15) + inRChannel * ((~inRChannel & 0xffff) >> 15);

		L *= v;
		R *= v;

		int outLChannel = ((~L + 1) & 0xffff) * ((inLChannel & 0xffff) >> 15) + L * ((~inLChannel & 0xffff) >> 15);
		int outRChannel = ((~R + 1) & 0xffff) * ((inRChannel & 0xffff) >> 15) + R * ((~inRChannel & 0xffff) >> 15);

		sample = float4((outLChannel & 0x00ff) / 255.0, ((outLChannel >> 8) & 0xff) / 255.0, (outRChannel & 0x00ff) / 255.0, ((outRChannel >> 8) & 0xff) / 255.0);
	//  sample = float4((outLeftChannel & 0x00ff)*0, ((outLeftChannel >> 8) & 0xff)*0 + 1, (outRightChannel & 0x00ff)*0 + 2, ((outRightChannel >> 8) & 0xff)*0 + 3) / 255.0;

		//color *= (input.tex.x);
		return sample;
	}*/

int2 getPCMData(float4 sample)
{
	int inLChannel = (int(ceil(sample.r * 255.0))) | ((int(ceil(sample.g * 255.0))) << 8);
	int inRChannel = (int(ceil(sample.b * 255.0))) | ((int(ceil(sample.a * 255.0))) << 8);

	return int2(inLChannel, inRChannel);
}
int2 absPCM(int2 inChannels)
{
	int L = ((~(inChannels.x - 1)) & 0xffff) * ((inChannels.x & 0xffff) >> 15) + inChannels.x * ((~inChannels.x & 0xffff) >> 15);
	int R = ((~(inChannels.y - 1)) & 0xffff) * ((inChannels.y & 0xffff) >> 15) + inChannels.y * ((~inChannels.y & 0xffff) >> 15);

	return int2(L, R);
}
int2 resignPCM(int2 channels, int2 origChannels)
{
	int outLChannel = ((~channels.x + 1) & 0xffff) * ((origChannels.x & 0xffff) >> 15) + channels.x * ((~origChannels.x & 0xffff) >> 15);
	int outRChannel = ((~channels.y + 1) & 0xffff) * ((origChannels.y & 0xffff) >> 15) + channels.y * ((~origChannels.y & 0xffff) >> 15);

	return int2(outLChannel, outRChannel);
}

int lowPass(int input, int lastInput, float a)
{
	return int(a * input + (1 - a) * lastInput);
}

//float4 getAudio(float2 pos, int index)
//{
//	float4 sample = tex2D(texture_sampler, float2(pos.x / 1.0, (pos.y / 3.0) + (index / 3.0)));
//
//	int2 inChannels = getPCMData(sample);
//	int2 channels = absPCM(inChannels);
//
//	float v = min(max(volume, 0.0), 1.0);
//	float l = min(max(low, 0.0), 1.0);
//
//	int2 outChannels = resignPCM(channels, inChannels);
//	sample = float4((outChannels.x & 0x00ff) / 255.0, ((outChannels.x >> 8) & 0xff) / 255.0, (outChannels.y & 0x00ff) / 255.0, ((outChannels.y >> 8) & 0xff) / 255.0);
//
//	return sample;
//}

float4 getAudio(float2 pos)
{
	//float[] v; = min(max(volume, 0.0), 1.0);
	
	float4 sample;
	int index;

	//1
	index = 0;
	sample = tex2D(texture_sampler, float2(pos.x / 1.0, (pos.y / 3.0) + (index / 3.0)));

	int2 inChannels1 = getPCMData(sample);
	int2 channels1 = absPCM(inChannels1);	
	//2
	index = 1;
	sample = tex2D(texture_sampler, float2(pos.x / 1.0, (pos.y / 3.0) + (index / 3.0)));

	int2 inChannels2 = getPCMData(sample);
	int2 channels2 = absPCM(inChannels2);
	//3
	index = 2;
	sample = tex2D(texture_sampler, float2(pos.x / 1.0, (pos.y / 3.0) + (index / 3.0)));

	int2 inChannels3 = getPCMData(sample);
	int2 channels3 = absPCM(inChannels3);

	channels1 *= volume[0];
	channels2 *= volume[1];
	channels3 *= volume[2];

	int2 outChannels1 = resignPCM(channels1, inChannels1);
	int2 outChannels2 = resignPCM(channels2, inChannels2);
	int2 outChannels3 = resignPCM(channels3, inChannels3);
	int2 outChannels = outChannels1 + outChannels2 + outChannels3;

	sample = float4((outChannels.x & 0x00ff) / 255.0, ((outChannels.x >> 8) & 0xff) / 255.0, (outChannels.y & 0x00ff) / 255.0, ((outChannels.y >> 8) & 0xff) / 255.0);

	return sample;
}

float4 AudioPS(VertexShaderOutput input) : COLOR
{	
	//float4 sample = tex2D(texture_sampler, float2(input.tex.x / 1.0, (input.tex.y / 3.0) + (1.0 / 3)));

	//int2 inChannels = getPCMData(sample);
	//int2 channels = absPCM(inChannels);
	//
	//float v = min(max(volume, 0.0), 1.0);
	//float l = min(max(low, 0.0), 1.0);

	//////--------------------------------------
	////float w = 60.0;
	////float h = 40.0;
	////
	////float passX =  (1.0 / w);
	////float passY =  (1.0 / h);
	////
	////float2 lastPos = input.tex;
	////lastPos.x -= passX;
	////
	////if (lastPos.x < 0.0)
	////{
	////	lastPos.y -= passY;
	////	lastPos.x = 1.0;
	////}
	////
	////int2 inLastChannel = getPCMData(tex2D(texture_sampler, lastPos));
	////int2 lastChannel = inLastChannel;//absPCM(inLastChannel);
	////
	////channels.x = (channels.x + lastChannel.x) / 2.0; //lowPass(channels.x, lastChannel.x, l);
	////channels.y = (channels.y + lastChannel.y) / 2.0; //lowPass(channels.y, lastChannel.y, l);
	////
	////inChannels = inLastChannel;
	///*inChannels.x = lowPass(inChannels.x, inLastChannel.x, l);
	//inChannels.y = lowPass(inChannels.y, inLastChannel.y, l);*/
	////_______________________________________

	////channels1.xy *= v;

	////int2 outChannels = channels;// resignPCM(channels, inChannels);
	////int2 outChannels = resignPCM(channels1.xy, channels1.ba);
	//int2 outChannels = resignPCM(channels, inChannels);
	//sample = float4((outChannels.x & 0x00ff) / 255.0, ((outChannels.x >> 8) & 0xff) / 255.0, (outChannels.y & 0x00ff) / 255.0, ((outChannels.y >> 8) & 0xff) / 255.0);

	

	return getAudio(input.tex);
}

technique BasicColorDrawing
{
	pass P0
	{
		VertexShader = compile VS_SHADERMODEL AudioVS();
		PixelShader = compile PS_SHADERMODEL AudioPS();
	}

	pass P1
	{
		VertexShader = compile VS_SHADERMODEL MainVS();
		PixelShader = compile PS_SHADERMODEL MainPS();
	}
};