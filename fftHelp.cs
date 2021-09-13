using System;

namespace Extreme.Numerics.QuickStart.CSharp
{
    // We'll need real and complex vectors.
    using Extreme.Mathematics;
    // The FFT classes reside in the Extreme.Mathematics.SignalProcessing
    // namespace.
    using Extreme.Mathematics.SignalProcessing;

    /// <summary>
    /// Illustrates the use of the FftProvider and Fft classes for computing 
    /// the Fourier transform of real and complex signals.
    /// </summary>
    class FourierTransforms
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main(string[] args)
        {
            // This QuickStart sample shows how to compute the Fouier
            // transform of real and complex signals.

            // Some vectors to play with:
            var r1 = Vector.Create(1000, i => 1.0 / (1 + i));
            var c1 = Vector.Create(1000,
                i => new Complex<double>(Math.Sin(0.03 * i), Math.Cos(0.07 * i)));

            var r2 = Vector.Create(new double[] { 1.0, 2.0, 3.0, 4.0 });
            var c2 = Vector.Create(new Complex<double>[] {
                new Complex<double>(1, 2), new Complex<double>(3, 4),
                new Complex<double>(5, 6), new Complex<double>(7, 8)
            });

            //
            // One-time FFT's
            //

            // The var and ComplexVector classes have static methods to compute FFT's:
            ComplexConjugateSignalVector<double> c3 = Vector.FourierTransform(r2);
            var r3 = Vector.InverseFourierTransform(c3);
            Console.WriteLine("fft(r2) = {0:F3}", c3);
            Console.WriteLine("ifft(fft(r2)) = {0:F3}", r3);
            // The ComplexConjugateSignalVector type represents a complex vector
            // that is the Fourier transform of a real signal. 
            // It enforces certain symmetry properties:
            Console.WriteLine("c3[i] == conj(c3[N-i]): {0} == conj({1})", c3[1], c3[3]);

            //
            // FFT objects
            //

            // FFT's require a fair bit of pre-computation. When performing
            // many transforms of the same length, it is more efficient 
            // to use an Fft object that caches these computations.

            // Here, we create an FFT implementation for a real signal:
            var realFft = Fft<double>.CreateReal(r1.Length);
            // For a complex to complex transform:
            var complexFft = Fft<double>.CreateComplex(c1.Length);

            // You can set the scale factor for the forward transform.
            // The default is 1/N.
            realFft.ForwardScaleFactor = 1.0 / Math.Sqrt(c1.Length);
            // and the backward transform, with default 1:
            realFft.BackwardScaleFactor = realFft.ForwardScaleFactor;

            // The ForwardTransform method performs a forward transform:
            var c4 = realFft.ForwardTransform(r1);
            Console.WriteLine("First 5 terms of fft(r1):");
            for (int i = 0; i < 5; i++)
                Console.WriteLine("   {0}: {1}", i, c4[i]);
            var c5 = complexFft.ForwardTransform(c1);
            Console.WriteLine("First 5 terms of fft(c1):");
            for (int i = 0; i < 5; i++)
                Console.WriteLine("   {0}: {1}", i, c5[i]);

            // ForwardTransform has many overloads for real to complex and
            // complex to complex transforms.

            // A one-sided transform returns only the first half of the FFT of
            // a real signal. The rest can be deduced from the symmetry properties.
            // Here's how to compute a one-sided FFT:
            var c6 = Vector.Create<Complex<double>>(r1.Length / 2 + 1);
            realFft.ForwardTransform(r1, c6, RealFftFormat.OneSided);

            // The BackwardTransform method has a similar set of overloads:
            var r4 = Vector.Create<double>(r1.Length);
            realFft.BackwardTransform(c6, r4, RealFftFormat.OneSided);

            // As the last step, we need to dispose the two FFT implementations.
            realFft.Dispose();
            complexFft.Dispose();

            //
            // 2D transforms
            //

            // 2D transforms are handled in a completely analogous way.
            var m = Matrix.Create(36, 56, (i, j) =>
                Math.Exp(-0.1 * i) * Math.Sin(0.01 * (i * i + j * j - i * j)));
            var mFft = Matrix.Create<Complex<double>>(m.RowCount, m.ColumnCount);

            using (var fft2 = Fft2D<double>.CreateReal(m.RowCount, m.ColumnCount))
            {
                fft2.ForwardTransform(m, mFft);

                Console.WriteLine("First few terms of fft(m):");
                for (int i = 0; i < 4; i++)
                {
                    string comma = string.Empty;
                    for (int j = 0; j < 4; j++)
                    {
                        Console.Write(comma);
                        Console.Write("{0}", mFft[i, j].ToString("F4"));
                        comma = ", ";
                    }
                    Console.WriteLine();
                }

                // and the backward transform:
                fft2.BackwardTransform(mFft, m);

                // Dispose is called automatically.
            }

            Console.Write("Press Enter key to exit...");
            Console.ReadLine();
        }
    }
}