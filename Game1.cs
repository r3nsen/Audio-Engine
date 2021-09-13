using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Input;

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace Audio_Engine
{
    public class Game1 : Game
    {
        private GraphicsDeviceManager _graphics;
        private GraphicsManager gm;

        public Game1()
        {
            _graphics = new GraphicsDeviceManager(this);
            Content.RootDirectory = "Content";
            IsMouseVisible = true;

            _graphics.PreferredBackBufferWidth = 960;
            _graphics.PreferredBackBufferHeight = 540;

            _graphics.PreferredDepthStencilFormat = DepthFormat.Depth24Stencil8;
            _graphics.GraphicsProfile = GraphicsProfile.HiDef;
            _graphics.PreferMultiSampling = false;
        }

        protected override void Initialize()
        {
            // TODO: Add your initialization logic here

            base.Initialize();
        }

        struct sample
        {
            public byte lowL;
            public byte highL;
            public byte lowR;
            public byte highR;
        }
        struct fourrier
        {
            public double real;
            public double imag;
            public double freq;
            public double amp;
            public double phase;

            public static fourrier operator *(fourrier a, fourrier b)
            {
                fourrier retorno = new fourrier();

                retorno.real = a.real * b.real - a.imag * b.imag;
                retorno.imag = a.real * b.imag + a.imag * b.real;

                return retorno;
            }

            public static fourrier operator /(fourrier a, int b)
            {
                fourrier retorno = new fourrier();

                retorno.real = a.real / b;
                retorno.imag = a.imag / b;

                return retorno;
            }

            public static fourrier operator +(fourrier a, fourrier b)
            {
                fourrier retorno = new fourrier();

                retorno.real = a.real + b.real;
                retorno.imag = a.imag + b.imag;

                return retorno;
            }

            public static fourrier operator -(fourrier a, fourrier b)
            {
                fourrier retorno = new fourrier();

                retorno.real = a.real - b.real;
                retorno.imag = a.imag - b.imag;

                return retorno;
            }

            internal void getData()
            {
                //   freq = k;
                amp = Math.Sqrt(real * real + imag * imag);
                if (real != 0)
                    phase = Math.Atan2(imag, real);
                else
                    phase = double.NaN;
            }
        }
        Complex[] FFT(Complex[] array)
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

            Complex[] outPar = FFT(par);
            Complex[] outImpar = FFT(impar);

            double angle = 2 * PI / n;

            Complex[] outArray = new Complex[n];

            for (int i = 0; 2 * i < n; i++)
            {
                outArray[i] = (outPar[i] + outImpar[i] * (new Complex(Math.Cos(angle * i), Math.Sin(angle * i))));
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
        void FFTOLD(ref fourrier[] array, bool invert = false)
        {
            double PI = Math.Acos(-1);
            int size = array.Length;
            if (size == 1) return;

            fourrier[] par = new fourrier[size / 2];
            fourrier[] impar = new fourrier[size / 2];

            for (int i = 0; 2 * i < size; i++)
            {
                par[i] = array[2 * i];
                impar[i] = array[2 * i + 1];
            }
            FFTOLD(ref par, invert);
            FFTOLD(ref impar, invert);

            double angle = 2 * PI / size * (invert ? -1 : 1);

            fourrier w, wn;
            w = new fourrier() { real = 1, imag = 1 };
            wn = new fourrier() { real = Math.Cos(angle), imag = Math.Sin(angle) };

            for (int i = 0; 2 * i < size; i++)
            {
                array[i] = par[i] + w * impar[i];
                array[i + size / 2] = par[i] - w * impar[i];

                if (invert)
                {
                    array[i] /= 2;
                    array[i + size / 2] /= 2;
                }
                w *= wn;
            }

        }
        fourrier[] DFT1(double[] array, int sampleRate = 48000)
        {
            int N = array.Length;
            fourrier[] X = new fourrier[N];

            for (int k = 1; k < (N + 0) / 1; k++)
            {
                double real = 0;
                double imag = 0;
                for (int i = 0; i < N; i++)
                {
                    double angle = (2 * Math.PI * k * i) / N;
                    real += (+array[i] * Math.Cos(angle));
                    imag += (-array[i] * Math.Sin(angle));
                }
                X[k].real = Math.Round((real * 2) / N, 4);
                X[k].imag = Math.Round((imag * 2) / N, 4);
                X[k].freq = k;
                X[k].amp = Math.Sqrt(X[k].real * X[k].real + X[k].imag * X[k].imag);
                X[k].phase = Math.Atan2(X[k].imag, X[k].real);
            }

            return X;
        }

        public double[] MDCT(double[] input)
        {
            int size = input.Length;
            double PI = 3.141592653589793238462643383279;

            double[] output = new double[size / 2];

            for (int k = 0; k < output.Length; k++)
            {
                for (int n = 0; n < input.Length; n++)
                {
                    output[k] += input[n] * Math.Cos((PI / output.Length) * (n + 0.5 + output.Length / 2.0) * (k + 0.5));
                }
            }

            return output;
        }

        public double[] IMDCT(double[] input)
        {
            int size = input.Length;
            double PI = 3.141592653589793238462643383279;

            double[] output = new double[size * 2];

            for (int n = 0; n < output.Length; n++)
            {
                for (int k = 0; k < input.Length; k++)
                {
                    output[n] += input[k] * Math.Cos((PI / output.Length) * (n + 0.5 + output.Length / 2.0) * (k + 0.5));
                }
                output[n] /= output.Length;
            }

            return output;
        }

        //fourrier[] FFT(short[] array, int sampleRate = 48000)
        //{
        //    int N = array.Length;            
        //    if (N == 1)
        //    {
        //        double angle = (2 * Math.PI * k * 0) / N;
        //        fourrier[] temp = new fourrier[1];
        //        temp[0].real = array[0] * Math.Cos(angle);
        //        temp[0].imag = array[0] * Math.Sin(angle);
        //        return array[0];
        //    }

        //    double W = (2 * Math.PI) / N;
        //    short[] par = array.Where((x, i) => i % 2 == 0).ToArray();
        //    short[] impar = array.Where((x, i) => i % 2 != 0).ToArray();

        //    fourrier[] return1 = FFT(par);
        //    fourrier[] return2 = FFT(impar);

        //    fourrier[] Y = new fourrier[N];

        //    for (int i = 0; i < N / 2; i++)
        //    {
        //        Y[i].real = return1[i].real + 0;
        //        Y[i + N / 2].imag = return1[i].real + 0;
        //    }

        //    return null;
        //}

        short IDFT(fourrier[] input)
        {

            return 0;
        }

        //struct MyRand
        //{
        //    static uint u, v;
        //    public MyRand(uint _u)
        //    {
        //        u = _u;
        //        v = 124;
        //    }
        //    public uint Next()
        //    {
        //        v = 36969 * (v & 65535) + (v >> 16);
        //        u = 18000 * (u & 65535) + (u >> 16);
        //        return (v << 16) + (u & 65535);
        //    }
        //}

        //struct HalfSquareRand
        //{
        //    static int NumOfDigits = 4;
        //    static int trim = NumOfDigits / 2;
        //    static int[] range = { 1, 10, 100, 1000, 10000, 100000, 1000000 };
        //    static int next = 0;
        //    public HalfSquareRand(int v)
        //    {
        //        next = v;
        //    }
        //    public int Next()
        //    {
        //        int square = next * next;
        //        square = square / range[trim];
        //        for (int i = 0; i < NumOfDigits; i++)
        //        {
        //            next += (square & range[trim]) * range[i];
        //            square /= 10;
        //        }
        //        return next;
        //    }
        //}

        protected override void LoadContent()
        {        
            {
                ////int X0, X1, N, P1, P2, init;
                ////init = 79;
                ////X0 = init;
                ////N = 100;
                ////P1 = 263;
                ////P2 = 71;
                ////int[] count = new int[100];
                //Random rand1 = new Random(101);
                //Random rand2 = new Random(320);
                //Random rand3 = new Random(570);
                //MyRand rand4 = new MyRand(467);
                //HalfSquareRand rand5 = new HalfSquareRand(5790986);
                //int p = 387;
                //while (true)
                //{
                //    Debug.WriteLine("------------------------------------------");
                //    for (int i = 0; i < 200; i++)
                //    {
                //        //X0 = i;
                //        //X1 = (X0 * P1 + P2) % N;
                //        //Debug.WriteLine(X1);
                //        ////X0 = X1;
                //        //++count[X1];
                //        //if (count[X1] > 1)
                //        //    ;
                //        //    //if (X0 == init)
                //        //    //;
                //        //rand5 = new MyRand((uint)i);
                //        //rand4.Next();
                //        //rand4.Next();
                //        // Debug.WriteLine("{0} - {1} - {2}",rand1.Next(), rand2.Next(), rand3.Next());
                //        //p = (i * i) % 10000;
                //        p = (p * p * p) % 100000;
                //        Debug.WriteLine("({0}, {1})", i, (p / 10000d).ToString().Replace(',', '.'));

                //        //Debug.WriteLine("({0}, {1})", i, (rand1.Next()/1000000000d + 4).ToString().Replace(',', '.'));


                //    }
                //}
            }
            gm = new GraphicsManager(GraphicsDevice);
            gm.fx = Content.Load<Effect>("Effect");
            // gm.cor = Content.Load<Texture2D>("zerinho");

            for (int i = 0; i < 4; i++)
            {
                byte[] bytes = File.ReadAllBytes("Content/piano teste_Seq0" + (i + 1) + ".wav");
                gm.pieces[i].audio = bytes[44..];
                // gm.dSong.SubmitBuffer(gm.pieces[i].audio,0, gm.pieces[i].audio.Length);
            }
            //gm.vorbis[0] = new NVorbis.VorbisReader("Content/piano teste.ogg");

            //gm.vorbis[0] = new NVorbis.VorbisReader("Content/music test/ritmado 1.ogg");
            //gm.vorbis[1] = new NVorbis.VorbisReader("Content/music test/ritmado 2.ogg");
            //gm.vorbis[2] = new NVorbis.VorbisReader("Content/music test/ritmado 3.ogg");

            gm.vorbis[0] = new NVorbis.VorbisReader("Content/layers/music_Seq01_melodia.ogg");
            gm.vorbis[1] = new NVorbis.VorbisReader("Content/layers/music_Seq01_ornamento.ogg");
            gm.vorbis[2] = new NVorbis.VorbisReader("Content/layers/music_Seq01_bass.ogg");

            gm.vorbis[3] = new NVorbis.VorbisReader("Content/layers/music_Seq02_melodia.ogg");
            gm.vorbis[4] = new NVorbis.VorbisReader("Content/layers/music_Seq02_ornamento.ogg");
            gm.vorbis[5] = new NVorbis.VorbisReader("Content/layers/music_Seq02_bass.ogg");

            gm.chunkSize = gm.vorbis[0].TotalSamples;

            gm.SubmitVoidBuffer();
            gm.dSong.Play();
        }

        protected override void Update(GameTime gameTime)
        {
            if (GamePad.GetState(PlayerIndex.One).Buttons.Back == ButtonState.Pressed || Keyboard.GetState().IsKeyDown(Keys.Escape))
                Exit();           

            // TODO: Add your update logic here

            base.Update(gameTime);
        }
        float freq = 1;
        protected override void Draw(GameTime gameTime)
        {

            int w = _graphics.PreferredBackBufferWidth;
            int h = _graphics.PreferredBackBufferHeight;

            gm.Begin(0, 0, 60, 40);
            gm.Draw(0, 0, 60, 40);
            gm.AudioEnd();

            gm.Begin(0, 0, w, h);
            //float pi = 3.141592653589793238462643383279f;

            //byte[] buff = gm.bbbuffer;

            //if (Keyboard.GetState().IsKeyDown(Keys.Left))
            //    freq = Math.Max(1.0f, freq - 10.0f);
            //if (Keyboard.GetState().IsKeyDown(Keys.Right))
            //    freq = Math.Min(20000, freq + 10.0f);

            //float i = 0;

            //float wave = (float)Math.Sin(2f * pi * (i / 2400f) * freq);
            //float wawe = (float)Math.Sin(2f * pi * (i / 48000f) * freq);            

            //gm.DrawAudioFreq(w / 2, (h / 2), 100, buff, gm.volume[0]);

            //gm.DrawAudioFreq(w / 2 - (w / 3), (h / 3) * 2, 50, buff, gm.volume[0]);
            //gm.DrawAudioFreq(w / 2, h / 2, 80, gm.bbuffer[((gm.bbuffer.Length / 3) * 1)..((gm.bbuffer.Length / 3) * 2)], gm.volume[1]);
            //gm.DrawAudioFreq(w / 2 + (w / 3), (h / 3) * 2, 50, gm.bbuffer[((gm.bbuffer.Length / 3) * 2)..((gm.bbuffer.Length / 3) * 3)], gm.volume[2]);
            gm.DrawTrack(0, 0, 0);
            gm.End();

            base.Draw(gameTime);
        }
    }
    static class BoolExt
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static byte AsByte(this bool value) =>
            Unsafe.As<bool, byte>(ref value);
    }
}
