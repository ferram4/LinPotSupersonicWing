/*
Copyright (c) 2014, Michael Ferrara
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using System;

namespace SupersonicWing
{
    class Program
    {
        static void Main()
        {
            Console.Write("Choose Wing Geometry File: ");
            WingGeometry geometry = GeometryReader.ReadWingGeometryFromFile(Console.ReadLine());

            Console.Write("Choose Mach Number: ");
            double beta = double.Parse(Console.ReadLine());
            beta = Math.Sqrt(beta * beta - 1);

            Console.Write("Choose Angle Of Attack: ");
            double angleOfAttack = double.Parse(Console.ReadLine());

            Console.Write("Choose Sideslip Angle: ");
            double sideslip = double.Parse(Console.ReadLine());

            Console.Write("Choose Horizontal Grid Size: ");
            int gridSize = int.Parse(Console.ReadLine());

            Console.Write("Num Sims: ");
            int numSims = int.Parse(Console.ReadLine());
            System.Diagnostics.Stopwatch watch = new System.Diagnostics.Stopwatch();

            int i = 0;
            geometry.ScaleWingBased(sideslip * Math.PI / 180, angleOfAttack * Math.PI / 180, beta);
            WingSim sim = new WingSim(geometry, gridSize, beta);

            watch.Start();
            for (i = 0; i < numSims; i++)
            {
                geometry.ScaleWingBased(sideslip * Math.PI / 180, angleOfAttack * Math.PI / 180, beta);


                sim = new WingSim(geometry, gridSize, beta);
                sim.RunSim();
            }
            watch.Stop();

            Console.WriteLine("Sim Cells: " + sim.countCells);
            Console.WriteLine("Total elapsed time: " + watch.ElapsedMilliseconds + " ms");
            Console.WriteLine("Elapsed Time Per Sim: " + (double)watch.ElapsedMilliseconds / numSims + " ms");

            sim.DumpToFile();

            sim.IntegrateAndPrintCoefficients(angleOfAttack * Math.PI / 180);
            Console.ReadKey();
        }
    }
}
