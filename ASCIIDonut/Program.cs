using System;
using System.IO;
using System.Threading;
using System.Threading.Tasks;

// https://www.a1k0n.net/2011/07/20/donut-math.html
// Simple version?
// https://youtu.be/DEqXNfs_HhY?t=53

namespace ASCIIDonut {
    public static class Program {
        private static byte[] conBuffer;
        private static Stream stdOut;

        private const double thetaSpacing = 0.07;
        private const double phiSpacing = 0.02;

        private static int conWidth;
        private static int conHeight;

        private const double R1 = 1.0;
        private const double R2 = 2.0;
        private const double K2 = 5.0;
        
        private static double K1;

        private static char[,] output;
        private static double[,] zBuffer;

        public static void Main(string[] args) {
            Console.CursorVisible = false;
            stdOut = Console.OpenStandardOutput();

            Task.Run(() => {
                double A = 1.0;
                double B = 1.0;

                while(!Console.KeyAvailable) {
                    if(conWidth  != Console.WindowWidth ||
                       conHeight != Console.WindowHeight) {
                        conWidth = Console.WindowWidth;
                        conHeight = Console.WindowHeight;

                        // Calculate K1 based on screen size: the maximum x-distance occurs
                        // roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
                        // want that to be displaced 3/8ths of the width of the screen, which
                        // is 3/4th of the way from the center to the side of the screen.
                        // conWidth*3/8 = K1*(R1+R2)/(K2+0)
                        // conWidth*K2*3/(8*(R1+R2)) = K1
                        K1 = Math.Min(conWidth, conHeight) * K2 * 3.0 / (8.0 * (R1 + R2));

                        output = new char[conWidth, conHeight]; // = ' ';
                        zBuffer = new double[conWidth, conHeight]; // = 0;
                        conBuffer = new byte[conWidth * conHeight];

                        for(int y = 0; y < conHeight; y++)
                            for(int x = 0; x < conWidth; x++)
                                output[x, y] = ' ';
                    }

                    A += 0.07;
                    B += 0.03;
                    RenderFrame(A, B);

                    Thread.Sleep(25);
                }
            }).Wait();
            stdOut.Close();
            
            Console.ReadKey(true);
            Console.CursorVisible = true;
        }

        private static void RenderFrame(double A, double B) {
            // precompute sines and cosines of A and B
            double cosA = Math.Cos(A), sinA = Math.Sin(A);
            double cosB = Math.Cos(B), sinB = Math.Sin(B);

            // theta goes around the cross-sectional circle of a torus
            for(double theta = 0; theta < 2 * Math.PI; theta += thetaSpacing) {
                // precompute sines and cosines of theta
                double cosTheta = Math.Cos(theta), sinTheta = Math.Sin(theta);

                // phi goes around the center of revolution of a torus
                for(double phi = 0; phi < 2 * Math.PI; phi += phiSpacing) {
                    // precompute sines and cosines of phi
                    double cosPhi = Math.Cos(phi), sinPhi = Math.Sin(phi);

                    // the x,y coordinate of the circle, before revolving (factored
                    // out of the above equations)
                    double circleX = R2 + R1 * cosTheta;
                    double circleY = R1 * sinTheta;

                    // final 3D (x,y,z) coordinate after rotations, directly from
                    // our math above
                    double x = circleX * (cosB * cosPhi + sinA * sinB * sinPhi)
                             - circleY * cosA * sinB;
                    double y = circleX * (sinB * cosPhi - sinA * cosB * sinPhi)
                             + circleY * cosA * cosB;
                    double z = K2 + cosA * circleX * sinPhi + circleY * sinA;
                    double ooz = 1.0 / z;  // "one over z"

                    // calculate luminance. ugly, but correct.
                    double L = cosPhi * cosTheta * sinB - cosA * cosTheta * sinPhi -
                               sinA * sinTheta + cosB * (cosA * sinTheta - 
                                                         cosTheta * sinA * sinPhi);
                    // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
                    // is pointing away from us, so we won't bother trying to plot it.

                    //if(L >= 0) {
                        // x and y projection. Note that y is negated here, because y
                        // goes up in 3D space but down on 2D displays.
                        int xp = (int)(conWidth / 2.0 + K1 * ooz * x);
                        int yp = (int)(conHeight / 2.0 - K1 * ooz * y);

                        // test against the z-buffer.  larger 1/z means the pixel is
                        // closer to the viewer than what's already plotted.
                        if(xp >= 0 && xp < conWidth && 
                           yp >= 0 && yp < conHeight &&
                           ooz > zBuffer[xp, yp]) {

                            zBuffer[xp, yp] = ooz;
                            int luminanceIndex = (int)Math.Max(0, L * 8);
                            // luminanceIndex is now in the range 0..11 (8*sqrt(2) = 11.3)
                            // now we lookup the character corresponding to the
                            // luminance and plot it in our output:
                            output[xp, yp] = ".,-~:;=!*#$@"[luminanceIndex];
                        }
                    //}
                }
            }

            // now, dump output[] to the screen.
            // bring cursor to "home" location, in just about any currently-used
            // terminal emulation mode
            Console.SetCursorPosition(0, 0);
            for(int y = 0; y < conHeight; y++) {
                for(int x = 0; x < conWidth; x++) {
                    conBuffer[x + y * conWidth] = (byte)output[x, y];

                    output[x, y] = ' ';
                    zBuffer[x, y] = 0;
                }
            }
            stdOut.Write(conBuffer, 0, conBuffer.Length - 1);
        }
    }
}