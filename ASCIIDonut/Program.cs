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
        private static Stream stdout;

        private const double theta_spacing = 0.07;
        private const double phi_spacing = 0.02;

        private static int screen_width;
        private static int screen_height;

        private const double R1 = 1.0;
        private const double R2 = 2.0;
        private const double K2 = 5.0;
        // Calculate K1 based on screen size: the maximum x-distance occurs
        // roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
        // want that to be displaced 3/8ths of the width of the screen, which
        // is 3/4th of the way from the center to the side of the screen.
        // screen_width*3/8 = K1*(R1+R2)/(K2+0)
        // screen_width*K2*3/(8*(R1+R2)) = K1
        private static double K1;

        private static char[,] output;
        private static double[,] zbuffer;

        public static void Main(string[] args) {
            stdout = Console.OpenStandardOutput();

            Task.Run(() => {
                double A = 1.0;
                double B = 1.0;

                while(!Console.KeyAvailable) {
                    if(screen_width != Console.WindowWidth || screen_height != Console.WindowHeight) {
                        Console.CursorVisible = false;

                        screen_width = Console.WindowWidth;
                        screen_height = Console.WindowHeight;
                        K1 = Math.Min(screen_width, screen_height) * K2 * 3.0 / (8.0 * (R1 + R2));

                        output = new char[screen_width, screen_height]; // = ' ';
                        zbuffer = new double[screen_width, screen_height]; // = 0;

                        conBuffer = new byte[screen_width * screen_height];
                    }

                    A += 0.07;
                    B += 0.03;
                    RenderFrame(A, B);

                    Thread.Sleep(25);
                }
            }).Wait();
            Console.ReadKey(true);
            Console.CursorVisible = true;

            stdout.Close();
        }

        private static void RenderFrame(double A, double B) {
            // precompute sines and cosines of A and B
            double cosA = Math.Cos(A), sinA = Math.Sin(A);
            double cosB = Math.Cos(B), sinB = Math.Sin(B);

            for(int x = 0; x < screen_width; x++) {
                for(int y = 0; y < screen_height; y++) {
                    output[x, y] = ' ';
                    zbuffer[x, y] = 0;
                }
            }

            // theta goes around the cross-sectional circle of a torus
            for(double theta = 0; theta < 2 * Math.PI; theta += theta_spacing) {
                // precompute sines and cosines of theta
                double costheta = Math.Cos(theta), sintheta = Math.Sin(theta);

                // phi goes around the center of revolution of a torus
                for(double phi = 0; phi < 2 * Math.PI; phi += phi_spacing) {
                    // precompute sines and cosines of phi
                    double cosphi = Math.Cos(phi), sinphi = Math.Sin(phi);

                    // the x,y coordinate of the circle, before revolving (factored
                    // out of the above equations)
                    double circlex = R2 + R1 * costheta;
                    double circley = R1 * sintheta;

                    // final 3D (x,y,z) coordinate after rotations, directly from
                    // our math above
                    double x = circlex * (cosB * cosphi + sinA * sinB * sinphi)
                             - circley * cosA * sinB;
                    double y = circlex * (sinB * cosphi - sinA * cosB * sinphi)
                             + circley * cosA * cosB;
                    double z = K2 + cosA * circlex * sinphi + circley * sinA;
                    double ooz = 1.0 / z;  // "one over z"

                    // calculate luminance.  ugly, but correct.
                    double L = cosphi * costheta * sinB - cosA * costheta * sinphi -
                               sinA * sintheta + cosB * (cosA * sintheta - 
                                                         costheta * sinA * sinphi);
                    // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
                    // is pointing away from us, so we won't bother trying to plot it.

                    //if(L >= 0) {
                        // x and y projection. Note that y is negated here, because y
                        // goes up in 3D space but down on 2D displays.
                        int xp = (int)(screen_width / 2.0 + K1 * ooz * x);
                        int yp = (int)(screen_height / 2.0 - K1 * ooz * y);

                        // test against the z-buffer.  larger 1/z means the pixel is
                        // closer to the viewer than what's already plotted.
                        if(xp >= 0 && xp < screen_width && yp >= 0 && yp < screen_height &&
                            ooz > zbuffer[xp, yp]) {

                            zbuffer[xp, yp] = ooz;
                            int luminance_index = (int)Math.Max(0, L * 8);
                            // luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
                            // now we lookup the character corresponding to the
                            // luminance and plot it in our output:
                            output[xp, yp] = ".,-~:;=!*#$@"[luminance_index];
                        }
                    //}
                }
            }

            // now, dump output[] to the screen.
            // bring cursor to "home" location, in just about any currently-used
            // terminal emulation mode
            Console.SetCursorPosition(0, 0);
            for(int j = 0; j < screen_height; j++) {
                for(int i = 0; i < screen_width; i++) {
                    //Console.Write(output[i, j]);
                    conBuffer[i + j * screen_width] = (byte)output[i, j];
                }
            }
            stdout.Write(conBuffer, 0, conBuffer.Length - 1);
        }
    }
}