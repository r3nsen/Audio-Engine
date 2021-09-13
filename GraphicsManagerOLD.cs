using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Audio;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Input;

using NVorbis;

using System;
using System.Diagnostics;
using System.Numerics;

using Vector2 = Microsoft.Xna.Framework.Vector2;
using Vector3 = Microsoft.Xna.Framework.Vector3;

namespace Audio_Engine
{
    internal class GraphicsManager
    {
        private GraphicsDevice gm;

        Matrix view;
        Matrix projection;

        public Effect fx;
        public Texture2D cor;
        RenderTarget2D render;

        RenderTarget2D audioBufferTexture;
        public Texture2D bufferTexture;

        private Texture2D blank;

        int samplerPointer;

        public DynamicSoundEffectInstance dSong;
        public struct Pieces
        {
            public byte[] audio;
        }
        public VorbisReader[] vorbis = new VorbisReader[10];

        public Pieces[] pieces = new Pieces[4];
        int currentPiece = 0;
        int sampleRate = 48000;

        public float[] volume = { 1.0f, 1.0f, 1.0f };
        float low = 1.0f;

        VertexPositionColorTexture[] vertex = new VertexPositionColorTexture[12000];
        short[] index = new short[(12000 * 3) / 2];
        int vertexCount, indexCount;

        public GraphicsManager(GraphicsDevice graphicsDevice)
        {
            this.gm = graphicsDevice;

            cor = new Texture2D(gm, 480, 480);
            render = new RenderTarget2D(gm, 480, 480);

            bufferTexture = new Texture2D(gm, 60, 40 * (3));
            audioBufferTexture = new RenderTarget2D(gm, 60, 40);

            dSong = new DynamicSoundEffectInstance(sampleRate, AudioChannels.Stereo);
            blank = new Texture2D(graphicsDevice, 1, 1);
            blank.SetData(new Color[] { new Color(220, 220, 220) });
        }

        public void Begin(int x, int y, int w, int h)
        {
            Vector3 pos = new Vector3(0, 0, 3);
            Vector3 target = new Vector3(0, 0, 0);
            Vector3 up = Vector3.Up;

            Matrix.CreateLookAt(ref pos, ref target, ref up, out view);
            Matrix.CreateOrthographicOffCenter(x, x + w, y + h, y, -500, 500, out projection);
        }

        public void Draw(int x, int y, int w, int h)
        {
            index[indexCount++] = (short)(vertexCount + 0);
            index[indexCount++] = (short)(vertexCount + 1);
            index[indexCount++] = (short)(vertexCount + 2);

            index[indexCount++] = (short)(vertexCount + 0);
            index[indexCount++] = (short)(vertexCount + 3);
            index[indexCount++] = (short)(vertexCount + 2);


            vertex[vertexCount++] = new VertexPositionColorTexture(new Vector3(x + 0, y + 0, 0), Color.Transparent, new Vector2(0, 0));
            vertex[vertexCount++] = new VertexPositionColorTexture(new Vector3(x + w, y + 0, 0), Color.Transparent, new Vector2(1.0f, 0));
            vertex[vertexCount++] = new VertexPositionColorTexture(new Vector3(x + w, y + h, 0), Color.Transparent, new Vector2(1.0f, 1.0f));
            vertex[vertexCount++] = new VertexPositionColorTexture(new Vector3(x + 0, y + h, 0), Color.Transparent, new Vector2(0, 1.0f));
        }

        Span<Complex> FFT(Span<Complex> array)
        {
            double PI = Math.Acos(-1);
            int n = array.Length;

            if (n == 1) return array;

            Span<Complex> par = stackalloc Complex[n / 2];
            Span<Complex> impar = stackalloc Complex[n / 2];

            for (int i = 0; 2 * i < n; i++)
            {
                par[i] = array[2 * i];
                impar[i] = array[2 * i + 1];
            }

            Span<Complex> outPar = FFT(par);
            Span<Complex> outImpar = FFT(impar);

            double angle = (2 * PI / n);
            float freqMul = 1;// 20000f / (n);
            Span<Complex> outArray = new Complex[n];

            for (int i = 0; 2 * i < n; i++)
            {
                outArray[i] = (outPar[i] + outImpar[i] * (new Complex(Math.Cos(angle * (i * freqMul)), Math.Sin(angle * (i * (freqMul))))));
                outArray[i + n / 2] = (outPar[i] - outImpar[i] * (new Complex(Math.Cos(angle * i), Math.Sin(angle * i))));
            }
            return outArray;
        }
        Complex[] IFFT(Complex[] array)
        {
            double PI = Math.Acos(-1);
            int n = array.Length;

            if (n == 1) return array;

            Complex[] par = new Complex[n / 2];
            Complex[] impar = new Complex[n / 2];

            for (int i = 0; 2 * i < n; i++)
            {
                par[i] = array[2 * i];
                impar[i] = array[2 * i + 1];
            }

            Complex[] outPar = IFFT(par);
            Complex[] outImpar = IFFT(impar);

            double angle = -2 * PI / n;

            Complex[] outArray = new Complex[n];

            for (int i = 0; 2 * i < n; i++)
            {
                outArray[i] = (outPar[i] + outImpar[i] * (new Complex(Math.Cos(angle * i), Math.Sin(angle * i))));
                outArray[i + n / 2] = (outPar[i] - outImpar[i] * (new Complex(Math.Cos(angle * i), Math.Sin(angle * i))));
            }
            return outArray;
        }

        public short byteToShort(byte[] b)
        {
            return (short)(b[0] | (b[1] << 8));
        }
        public (byte a, byte b) shortToByte(short s)
        {
            return ((byte)(s & 0xff), (byte)((s >> 8) & 0xff));
        }

        public void DrawAudio(int x, int y, int r, byte[] buffer, float vol)
        {
            vertex[vertexCount++] = new VertexPositionColorTexture(new Vector3(x, y, 0), Color.White, Vector2.Zero);
            int len = 500;// buffer.Length / 20;
            int currRadius = 0;
            for (int i = 0; i < buffer.Length / 2; i++)
                currRadius += Math.Abs(Math.Max(byteToShort(buffer[(2 * i)..(2 * i + 2)]), (short)-32767));

            currRadius /= buffer.Length / 2;
            currRadius /= 30;
            currRadius = (int)(currRadius * vol);
            for (int i = 0; i < len; i++)
            {
                index[indexCount++] = (short)(vertexCount + i + 0);
                index[indexCount++] = (short)(vertexCount + (i + 1) % len);
                index[indexCount++] = (short)(vertexCount - 1 + 0);



                //int mediaNum = buffer.Length / len;
                //for (int j = 0; j < mediaNum && i * mediaNum + j < buffer.Length; j++)
                //{
                //    currRadius += buffer[i * mediaNum + j];
                //}
                //currRadius /= mediaNum;

                vertex[vertexCount + i] = new VertexPositionColorTexture(
                    new Vector3(x, y, 0) + new Vector3(
                        (float)((currRadius + r) * Math.Sin(((float)i / len) * Math.PI * 2)),
                        (float)((currRadius + r) * Math.Cos(((float)i / len) * Math.PI * 2)),
                        0),
                    Color.White, Vector2.Zero);
            }
            vertexCount += len;
        }

        public void DrawAudioFreq(int x, int y, int r, byte[] buffer, float vol)
        {
            Span<Complex> buff = new Complex[buffer.Length / 4];
            for (int i = 0; i < buff.Length; i++)
            {
                buff[i] = new Complex((byteToShort(buffer[(4 * i)..(4 * i + 2)]) + byteToShort(buffer[(4 * i + 2)..(4 * i + 4)])) / 2.0, 0);
            }
            buff = FFT(buff[0..2048]);

            vertex[vertexCount++] = new VertexPositionColorTexture(new Vector3(x, y, 0), Color.White, Vector2.Zero);

            int len = 60;
            int[] radius = new int[len];
            int pointsBetween = 19;

            {
                int currRadius;

                for (int i = 0; i < len / 2; i++)
                {
                    double log20k = Math.Log(20000, 2) - 5; // l

                    double passo = log20k / (len / 2);

                    double init = (passo * (i - 1));
                    double end = (passo * (i));

                    init = Math.Pow(2, init + 5) / 10;
                    end = Math.Pow(2, end + 5) / 10;

                    float re = 0;
                    float im = 0;

                    for (int j = (int)init; j < end && j < buff.Length; j++)
                    {
                        re += (float)buff[j].Real;
                        im += (float)buff[j].Imaginary;
                    }

                    currRadius = (int)Math.Sqrt(re * re + im * im);
                    currRadius /= 20000;
                    radius[i] = radius[^(i + 1)] = currRadius;
                }
            }
            for (int i = 0; i < len + len * pointsBetween; i++)
            {
                int pts = pointsBetween;
                int origVert = i / (pointsBetween + 1);

                int currVert = vertexCount + i;

                index[indexCount++] = (short)(vertexCount + i + 0);
                index[indexCount++] = (short)(vertexCount + (i + 1) % (len + len * pointsBetween));
                index[indexCount++] = (short)(vertexCount - 1 + 0);

                vertex[vertexCount + i] = new VertexPositionColorTexture(
                    new Vector3(x, y, 0) + new Vector3(
                        (float)((radius[origVert] + r) * Math.Sin(((float)origVert / len) * Math.PI * 2)),
                        (float)((radius[origVert] + r) * Math.Cos(((float)origVert / len) * Math.PI * 2)),
                        0), Color.White, Vector2.Zero);

                (float, float) bezier4((float, float) init, (float, float) end, (float, float) center, float st)
                {
                    (float, float) B = ((init.Item1 + end.Item1) * .5f, (init.Item2 + end.Item2) * .5f);
                    (float, float) D = B;

                    float a, b, c, d, e;
                    float et = (1 - st);
                    a = 1 * st * st * st * st;
                    b = 4 * st * st * st * et;
                    c = 6 * st * st * et * et;
                    d = 4 * st * et * et * et;
                    e = 1 * et * et * et * et;

                    return (
                        (a * init.Item1 + b * B.Item1 + c * center.Item1 + d * D.Item1 + e * end.Item1),
                        (a * init.Item2 + b * B.Item2 + c * center.Item2 + d * D.Item2 + e * end.Item2)
                        );
                }

                while (pts > 0)
                {
                    ++i; --pts;

                    index[indexCount++] = (short)(vertexCount + i + 0);
                    index[indexCount++] = (short)(vertexCount + (i + 1) % (len + len * pointsBetween));
                    index[indexCount++] = (short)(vertexCount - 1 + 0);

                    int nextR = (origVert + 1) % radius.Length;
                    float f = (1.0f * (pts + 1)) / (pointsBetween + 1);

                    (float x, float y) pos = bezier4(
                        (vertex[currVert].Position.X, vertex[currVert].Position.Y),
                        ((float)(x + (radius[nextR] + r) * Math.Sin(((float)nextR / len) * Math.PI * 2)),
                         (float)(y + (radius[nextR] + r) * Math.Cos(((float)nextR / len) * Math.PI * 2))),
                        ((float)(x + r * Math.Sin(((2f * origVert + 1) / (len * 2f)) * Math.PI * 2)),
                         (float)(y + r * Math.Cos(((2f * origVert + 1) / (len * 2f)) * Math.PI * 2))
                        ), f
                    );

                    vertex[vertexCount + i] = new VertexPositionColorTexture(
                    new Vector3(pos.x, pos.y, 0),
                    Color.White, Vector2.Zero);
                }
            }
            vertexCount += len + len * pointsBetween;
        }

        public void SubmitVoidBuffer()
        {
            byte[] data = new byte[2400 * 4];
            dSong.SubmitBuffer(data, 0, 2400 * 4);
        }

        public void AudioEndOLD()
        {
            if (dSong.PendingBufferCount < pieces.Length)
            {
                cor.SetData(pieces[currentPiece].audio);

                gm.BlendState = BlendState.AlphaBlend;
                gm.RasterizerState = RasterizerState.CullNone;
                gm.SamplerStates[0] = SamplerState.PointClamp;

                fx.Parameters["WorldViewProjection"].SetValue(view * projection);
                fx.Parameters["tex"].SetValue(cor);
                fx.Parameters["sampleRate"].SetValue(sampleRate);
                fx.Parameters["volume"].SetValue(0.5f);

                gm.SetRenderTarget(render);
                gm.Clear(Color.Transparent);

                fx.CurrentTechnique.Passes[0].Apply();
                gm.DrawUserIndexedPrimitives(PrimitiveType.TriangleList, vertex, 0, vertexCount, index, 0, indexCount / 3);

                byte[] bytes = new byte[960 * 960];

                render.GetData(bytes, 0, 960 * 960);

                dSong.SubmitBuffer(bytes);
                ++currentPiece;
                currentPiece %= pieces.Length;
            }
        }

        float[] H = { 1.0f };
        int lastFA, lastFB, lastatt;
        private void CalcCoeff(int fa, int fb, int M = 21, int att = 60)
        {
            if (lastFA == fa && lastFB == fb && lastatt == att)
                return;

            lastFA = fa;
            lastFB = fb;
            lastatt = att;

            float a = 0;
            int fs = 48000;
            //int M = 21;

            if ((M % 2) == 0) M++;

            if (att > 50) a = 0.1102f * (att - 8.7f);
            else if (att >= 21) a = (float)(0.5842f * Math.Pow((att - 21), 0.4)) + 0.07886f * (att - 21);

            //int M = (int)((att - 8) / ((14.36 * df) / fs));// lenght

            Span<float> w = stackalloc float[M];
            Span<float> A = stackalloc float[M];

            if (H == null || H.Length != (M + 1) / 2)
                H = new float[(M + 1) / 2];

            float Ino(float x)
            {
                float d = 0;
                float ds = 1;
                float s = 1;

                do
                {
                    d += 2;
                    ds *= (x * x) / (d * d);
                    s += ds;
                } while (ds > (s / 1000000f));

                return s;
            }

            int np = (M - 1) / 2;
            float pi = 3.141592653589793238462643383279f;

            for (int i = 0; i <= np; i++)
            {
                float b = (float)(i * i) / (np * np);//(float)(i - np) / np;
                float up = Ino((float)(a * Math.Sqrt(1 - b /** b*/)));
                float down = Ino(a);
                w[i] = (up / down);

                if (i == 0) A[0] = (2f * ((fb - fa) / (float)fs));
                else
                    A[i] = (float)((Math.Sin(2 * pi * i * ((float)fb / fs)) - Math.Sin(2 * pi * i * ((float)fa / fs))) / (pi * i));

                H[i] = A[i] * w[i];
            }
        }

        byte[] bbuffer = new byte[60 * 40 * 2 * 2 * (3)];
        public byte[] bbbuffer = new byte[60 * 40 * 2 * 2];
        int freq = 22;//440;
        int coeffNum = 10;
        float[] coeficients = { 1f };

        int _fa = 0, _fb = 22000, _att = 0;

        public void AudioEnd()
        {
            if (Keyboard.GetState().IsKeyDown(Keys.Q)) volume[0] = Math.Min(1.0f, volume[0] + 0.01f);
            if (Keyboard.GetState().IsKeyDown(Keys.W)) volume[0] = Math.Max(0.0f, volume[0] - 0.01f);

            if (Keyboard.GetState().IsKeyDown(Keys.A)) volume[1] = Math.Min(1.0f, volume[1] + 0.01f);
            if (Keyboard.GetState().IsKeyDown(Keys.S)) volume[1] = Math.Max(0.0f, volume[1] - 0.01f);

            if (Keyboard.GetState().IsKeyDown(Keys.Z)) volume[2] = Math.Min(1.0f, volume[2] + 0.01f);
            if (Keyboard.GetState().IsKeyDown(Keys.X)) volume[2] = Math.Max(0.0f, volume[2] - 0.01f);

            //if (Keyboard.GetState().IsKeyDown(Keys.Left)) low = Math.Min(1.0f, low + 0.01f);
            //if (Keyboard.GetState().IsKeyDown(Keys.Right)) low = Math.Max(0.0f, low - 0.01f);

            if (dSong.PendingBufferCount < 2)
            {
                gm.BlendState = BlendState.AlphaBlend;
                gm.RasterizerState = RasterizerState.CullNone;
                gm.SamplerStates[0] = SamplerState.PointClamp;

                //  int count = 0;
                float[] fbuffer = new float[60 * 40 * 2 * (3)];

            a:
                int count1 = vorbis[0].ReadSamples(fbuffer, (fbuffer.Length / 3) * 0, fbuffer.Length / 3);
                int count2 = vorbis[1].ReadSamples(fbuffer, (fbuffer.Length / 3) * 1, fbuffer.Length / 3);
                int count3 = vorbis[2].ReadSamples(fbuffer, (fbuffer.Length / 3) * 2, fbuffer.Length / 3);

                if (vorbis[0].IsEndOfStream)
                {
                    vorbis[0].SamplePosition = 0;
                    vorbis[1].SamplePosition = 0; //
                    vorbis[2].SamplePosition = 0;
                    goto a;
                }

                //for (int i = 0; i < fbuffer.Length; i++)
                //    fbuffer[i] = (float)Math.Sin(((i/2f)*Math.PI*2)/5f);
                //float last = 0;
                //float a = low;
                //for (int i = 0; i < fbuffer.Length; i++)
                //{
                //    //float temp = fbuffer[i];
                //    //fbuffer[i] = (float)(a * last + (1.0f - a)*fbuffer[i]);
                //    //last = fbuffer[i];
                //    last = last - (a * (last - fbuffer[i]));
                //    fbuffer[i] = last;
                //}

                //bbuffer = new byte[60 * 40 * 2 * 2 * (3)];
                for (int i = 0; i < count1; i++)
                {
                    bbuffer[2 * i + 0 + (60 * 40 * 2 * 2) * 0] = (byte)(((short)(fbuffer[i + (60 * 40 * 2) * 0] * short.MaxValue) & 0x00ff));
                    bbuffer[2 * i + 1 + (60 * 40 * 2 * 2) * 0] = (byte)(((short)(fbuffer[i + (60 * 40 * 2) * 0] * short.MaxValue) & 0xff00) >> 8);

                    bbuffer[2 * i + 0 + (60 * 40 * 2 * 2) * 1] = (byte)(((short)(fbuffer[i + (60 * 40 * 2) * 1] * short.MaxValue) & 0x00ff));
                    bbuffer[2 * i + 1 + (60 * 40 * 2 * 2) * 1] = (byte)(((short)(fbuffer[i + (60 * 40 * 2) * 1] * short.MaxValue) & 0xff00) >> 8);

                    bbuffer[2 * i + 0 + (60 * 40 * 2 * 2) * 2] = (byte)(((short)(fbuffer[i + (60 * 40 * 2) * 2] * short.MaxValue) & 0x00ff));
                    bbuffer[2 * i + 1 + (60 * 40 * 2 * 2) * 2] = (byte)(((short)(fbuffer[i + (60 * 40 * 2) * 2] * short.MaxValue) & 0xff00) >> 8);
                }

                bufferTexture.SetData(bbuffer, 0, 60 * 40 * 2 * 2 * (3));

                fx.Parameters["WorldViewProjection"].SetValue(view * projection);
                fx.Parameters["tex"].SetValue(bufferTexture);
                fx.Parameters["sampleRate"].SetValue(sampleRate);
                fx.Parameters["volume"].SetValue(volume);
                fx.Parameters["low"].SetValue(low);

                gm.SetRenderTarget(audioBufferTexture);
                gm.Clear(Color.Transparent);

                fx.CurrentTechnique.Passes[0].Apply();
                gm.DrawUserIndexedPrimitives(PrimitiveType.TriangleList, vertex, 0, vertexCount, index, 0, indexCount / 3);

                audioBufferTexture.GetData(bbbuffer, 0, 60 * 40 * 2 * 2);

                {
                    //if (Keyboard.GetState().IsKeyDown(Keys.Right)) ++coeffNum;
                    //if (Keyboard.GetState().IsKeyDown(Keys.Left)) --coeffNum;
                    //
                    //coeffNum = Math.Max(1, Math.Min(coeffNum, 100));
                    //
                    //float[] highPass = new float[] {
                    //    +0.400000f,
                    //    -0.289226f,
                    //    +0.077837f,
                    //    +0.041025f,
                    //    -0.035457f,
                    //    +0.000000f,
                    //    +0.008267f,
                    //    -0.002034f,
                    //    -0.000650f,
                    //    +0.000229f,
                    //    -0.000000f
                    //};
                    //float[] lowPass = new float[] {
                    //    0.083333f,
                    //    0.078709f,
                    //    0.066212f,
                    //    0.049353f,
                    //    0.032287f,
                    //    0.018251f,
                    //    0.008693f,
                    //    0.003343f,
                    //    0.000958f,
                    //    0.000170f,
                    //    0.000008f
                    // };
                    //float[] normalPass = new float[] { 1f };
                    //
                    //if (Keyboard.GetState().IsKeyDown(Keys.I))
                    //    coeficients = normalPass;
                    //if (Keyboard.GetState().IsKeyDown(Keys.O)) coeficients = highPass;
                    //if (Keyboard.GetState().IsKeyDown(Keys.P)) coeficients = lowPass;
                } // coment

                if (Keyboard.GetState().IsKeyDown(Keys.Down)) _att = (int)Math.Min(3000, _att + 1);
                if (Keyboard.GetState().IsKeyDown(Keys.Up)) _att = (int)Math.Max(21, _att - 1);

                if (Keyboard.GetState().IsKeyDown(Keys.J)) _fa = (int)Math.Min(22000, _fa + 100);
                if (Keyboard.GetState().IsKeyDown(Keys.H)) _fa = (int)Math.Max(0, _fa - 100);

                if (Keyboard.GetState().IsKeyDown(Keys.M)) _fb = (int)Math.Min(22000, _fb + 100);
                if (Keyboard.GetState().IsKeyDown(Keys.N)) _fb = (int)Math.Max(0, _fb - 100);

                CalcCoeff(_fa, _fb, 51, _att);
                coeficients = H;

                for (int i = 0; i < bbbuffer.Length / 4; i++)
                {
                    float left = 0;
                    float right = 0;

                    for (int j = 0; j < coeficients.Length; j++)
                    {
                        if (i - j >= 0)
                        {
                            left += (byteToShort(bbbuffer[(4 * (i - j) + 0)..(4 * (i - j) + 2)]) * coeficients[j] * 1);
                            right += (byteToShort(bbbuffer[(4 * (i - j) + 2)..(4 * (i - j) + 4)]) * coeficients[j] * 1);
                        }
                        if (i + j < bbbuffer.Length / 4 && j != 0)
                        {
                            left += (byteToShort(bbbuffer[(4 * (i + j) + 0)..(4 * (i + j) + 2)]) * coeficients[j] * 1);
                            right += (byteToShort(bbbuffer[(4 * (i + j) + 2)..(4 * (i + j) + 4)]) * coeficients[j] * 1);
                        }
                    }

                    (bbbuffer[4 * i + 0], bbbuffer[4 * i + 1]) = shortToByte((short)left);
                    (bbbuffer[4 * i + 2], bbbuffer[4 * i + 3]) = shortToByte((short)right);
                }
                dSong.SubmitBuffer(bbbuffer, 0, count1 * 2);
            }
            vertexCount = indexCount = 0;
        }

        public void End()
        {
            gm.SetRenderTarget(null);
            gm.Clear(new Color(20, 20, 20));


            if (vertexCount == 0) return;
            gm.BlendState = BlendState.NonPremultiplied;
            RasterizerState rast = new RasterizerState();
            rast.FillMode = FillMode.Solid;
            rast.CullMode = CullMode.None;
            gm.RasterizerState = rast;
            gm.SamplerStates[0] = SamplerState.PointClamp;

            fx.Parameters["WorldViewProjection"].SetValue(view * projection);
            fx.Parameters["tex"].SetValue(blank);

            fx.CurrentTechnique.Passes[1].Apply();
            gm.DrawUserIndexedPrimitives(PrimitiveType.TriangleList, vertex, 0, vertexCount, index, 0, indexCount / 3);

            vertexCount = indexCount = 0;
        }
    }
}