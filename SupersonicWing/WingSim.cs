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
using System.IO;
using System.Diagnostics;

namespace SupersonicWing
{
    class WingSim
    {
        WingGeometry wing;

        int jGridSize, iGridSize;
        double gridScale, invGridScale;
        double beta, twoDivBeta;

        const double INV_PI = 0.31830988618379067153776752674503;

        VortexCell[,] vortexGrid;
        Stopwatch watch;

        public WingSim(WingGeometry wing, int horizontalGridSize, double beta)
        {
            this.wing = wing;
            this.jGridSize = horizontalGridSize;
            gridScale = (double)horizontalGridSize / (wing.maxY - wing.minY);
            invGridScale = 1 / gridScale;
            iGridSize = (int)Math.Round(gridScale * (wing.maxX - wing.minX));

            this.beta = beta;
            twoDivBeta = 2 / beta;

            watch = new Stopwatch();
            InitializeGrid();
        }

        private void InitializeGrid()
        {
            watch.Start();
            vortexGrid = new VortexCell[this.iGridSize, this.jGridSize];

            int countCells = 0;

            Console.WriteLine(iGridSize + " " + jGridSize);
            for (int i = 0; i < iGridSize; i++)
                for (int j = 0; j < jGridSize; j++)
                {
                    VortexCell vortCell = new VortexCell();

                    Vector3 test = new Vector3(-(i + 0.5) * invGridScale + wing.maxX, (j + 0.5) * invGridScale + wing.minY, 0);

                    double tmp;
                    Triangle tri;
                    vortCell.exists = wing.TryFindSlopeAtPoint(test, out tmp, out tri);
                    vortCell.localSlope = tmp;

                    if (vortCell.exists)
                    {
                        countCells++;
                        vortCell.influenceFactor = CalculateInfluenceFactor(tri, test);
                    }

                    vortexGrid[i, j] = vortCell;

                }

            watch.Stop();
            Console.WriteLine("Grid init time elapsed: " + watch.ElapsedMilliseconds.ToString() + " ms");
            watch.Reset();
            Console.WriteLine(countCells + " vortex cells for sim.");
        }

        private double CalculateInfluenceFactor(Triangle tri, Vector3 test)
        {
            if (wing.TryFindSlopeAtPoint(test + new Vector3(1, 0, 0)))
                return 1;

            double forwardOffset = 0;

            if (tri.p0.x > test.x)
            {
                double x1;
                if ((tri.p0.y - test.y) * (tri.p1.y - test.y) < 0)
                    x1 = (tri.p0.x - tri.p1.x) / (tri.p0.y - tri.p1.y) * (test.y - tri.p1.y) + tri.p1.x;
                else
                    x1 = (tri.p0.x - tri.p2.x) / (tri.p0.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;

                forwardOffset = x1 - test.x;
            }
            else if (tri.p1.x > test.x)
            {
                double x1;
                if ((tri.p0.y - test.y) * (tri.p1.y - test.y) < 0)
                    x1 = (tri.p0.x - tri.p1.x) / (tri.p0.y - tri.p1.y) * (test.y - tri.p1.y) + tri.p1.x;
                else
                    x1 = (tri.p1.x - tri.p2.x) / (tri.p1.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;

                forwardOffset = x1 - test.x;
            }
            else
            {
                double x1;
                if ((tri.p0.y - test.y) * (tri.p2.y - test.y) < 0)
                    x1 = (tri.p0.x - tri.p2.x) / (tri.p0.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;
                else
                    x1 = (tri.p1.x - tri.p2.x) / (tri.p1.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;

                forwardOffset = x1 - test.x;
            }

            Triangle oldTri = tri;
            if (wing.TryFindSlopeAtPoint(test + new Vector3(forwardOffset + 1e-6, 0, 0), out tri) && !oldTri.Equals(tri))
                return CalculateInfluenceFactor(tri, test);
            else
                return forwardOffset;
        }

        public void RunSim()
        {
            watch.Start();
            for (int i = 0; i < this.iGridSize; i++)
            {
                for (int j = 0; j < this.jGridSize; j++)
                {
                    VortexCell vortCell = vortexGrid[i, j];
                    if (!vortCell.exists)
                        continue;

                    vortCell.deltaU = -twoDivBeta * vortCell.localSlope;
                    vortCell.deltaU += CalculateInfluence(i, j);
                    vortexGrid[i, j] = vortCell;
                }
                if (i - 1 > 0)
                {
                    for (int j = 0; j < this.jGridSize; j++)
                    {
                        VortexCell vortCell = vortexGrid[i - 1, j];
                        if (!vortCell.exists)
                            continue;

                        vortCell.deltaU = vortCell.deltaU * 0.5 * (1 + vortCell.influenceFactor / (1 + vortCell.influenceFactor)) + vortexGrid[i, j].deltaU * 0.5 / (1 + vortCell.influenceFactor);
                        vortexGrid[i - 1, j] = vortCell;
                    }
                    for (int j = 0; j < this.jGridSize; j++)
                    {
                        VortexCell vortCell = vortexGrid[i, j];
                        if (!vortCell.exists)
                            continue;

                        vortCell.deltaU = -twoDivBeta * vortCell.localSlope;
                        vortCell.deltaU += CalculateInfluence(i, j);
                        vortexGrid[i, j] = vortCell;
                    }
                }
            }
            watch.Stop();
            Console.WriteLine("Time elapsed: " + watch.ElapsedMilliseconds.ToString() + " ms");
            DumpToFile();
        }

        private void DumpToFile()
        {
            StreamWriter writer = new StreamWriter(new FileStream("sim.csv", FileMode.Create));

            for (int i = 0; i < this.iGridSize; i++)
            {
                for (int j = 0; j < this.jGridSize; j++)
                {
                    writer.Write(vortexGrid[i, j].deltaU * -2);
                    writer.Write(",");
                }
                writer.WriteLine();
            }

            writer.Close();
        }

        public void IntegrateAndPrintCoefficients(double angleOfAttackRad)
        {
            double cL, cD, area;
            cL = 0;
            area = 0;

            for (int i = 0; i < this.iGridSize; i++)
                for (int j = 0; j < this.jGridSize; j++)
                {
                    VortexCell vortCell = vortexGrid[i, j];
                    if (!vortCell.exists)
                        continue;

                    cL += vortCell.deltaU * -2 * vortCell.influenceFactor;
                    area += vortCell.influenceFactor;
                }

            cL /= area;
            cD = cL * Math.Tan(angleOfAttackRad);

            Console.WriteLine("C_L: " + cL);
            Console.WriteLine("C_D: " + cD);
            Console.WriteLine("L/D: " + cL / cD);
        }

        private double CalculateInfluence(int iIndex, int jIndex)
        {
            double influence = 0;
            for (int i = iIndex - 1; i >= 0; i--)         //Start from the points right in front of it and work forward to find all influences
            {
                int forwardIndex = iIndex - i;
                for (int j = jIndex - forwardIndex; j <= jIndex + forwardIndex; j++)        //Only get points in a triangle in front of the point in question
                {
                    if (j < 0 || j >= jGridSize)
                        continue;

                    VortexCell otherCell = vortexGrid[i, j];
                    if (!otherCell.exists)
                        continue;

                    double forwardBack = forwardIndex + 0.5;
                    double forwardBackSqr = forwardBack * forwardBack;

                    double sidePos = jIndex - j + 0.5;
                    double sideNeg = sidePos - 1;

                    double curInfluence = Math.Sqrt((forwardBackSqr - sideNeg * sideNeg)) / (forwardBack * sideNeg);
                    curInfluence -= Math.Sqrt((forwardBackSqr - sidePos * sidePos)) / (forwardBack * sidePos);

                    curInfluence *= otherCell.deltaU * otherCell.influenceFactor;

                    influence += curInfluence;
                }
            }

            return influence * INV_PI;
        }
    }

    struct VortexCell
    {
        public bool exists;
        public double influenceFactor;
        public double deltaU;
        public double localSlope;
    }
}
