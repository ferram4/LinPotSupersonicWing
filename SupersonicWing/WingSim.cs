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
        static double[,] influenceGrid;

        public int countCells = 0;

        public WingSim(WingGeometry wing, int horizontalGridSize, double beta)
        {
            this.wing = wing;
            this.jGridSize = horizontalGridSize;
            gridScale = (double)horizontalGridSize / (wing.maxY - wing.minY);
            invGridScale = 1 / gridScale;
            iGridSize = (int)Math.Round(gridScale * (wing.maxX - wing.minX));

            this.beta = beta;
            twoDivBeta = 2 / beta;

            InitializeVortexGrid();
            CheckInfluenceGrid();
        }

        private void CheckInfluenceGrid()
        {
            if (influenceGrid == null)
            {
                influenceGrid = new double[iGridSize, jGridSize];
                RecalculateInfluenceGrid();
            }
            else
            {
                int iLength, jLength;
                iLength = influenceGrid.GetLength(0);
                jLength = influenceGrid.GetLength(1);

                if (iLength < iGridSize || jLength < jGridSize)
                {
                    influenceGrid = new double[Math.Max(iGridSize, iLength), Math.Max(jGridSize, jLength)];
                    RecalculateInfluenceGrid();
                }
            }
        }

        private void RecalculateInfluenceGrid()
        {
            int iLength, jLength;
            iLength = influenceGrid.GetLength(0);
            jLength = influenceGrid.GetLength(1);

            Console.WriteLine("Recalcing Influences");

            for (int i = 1; i < iLength; i++)         //Start from the points right in front of it and work forward to find all influences
            {
                for (int j = 0; j <= i; j++)        //Only get points in a triangle in front of the point in question
                {
                    if (j >= jLength)
                        continue;

                    double forwardBack = i + 0.5;
                    double forwardBackSqr = forwardBack * forwardBack;

                    double sidePos = j + 0.5;
                    double sideNeg = sidePos - 1;

                    double curInfluence = Math.Sqrt((forwardBackSqr - sideNeg * sideNeg)) / (forwardBack * sideNeg);
                    curInfluence -= Math.Sqrt((forwardBackSqr - sidePos * sidePos)) / (forwardBack * sidePos);
                    curInfluence *= INV_PI;

                    influenceGrid[i, j] = curInfluence;
                }
            }
        }

        private void InitializeVortexGrid()
        {
            vortexGrid = new VortexCell[this.iGridSize, this.jGridSize];

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
                        vortCell.forwardVelInfluence = 0.5 * (1 + vortCell.influenceFactor / (1 + vortCell.influenceFactor));
                        vortCell.backwardVelInfluence = 0.5 / (1 + vortCell.influenceFactor);
                    }

                    vortexGrid[i, j] = vortCell;

                }
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
            double[] farInfluence = new double[jGridSize];      //used to store far influences, since they do not change in the relaxation procedure

            for (int i = 0; i < this.iGridSize; i++)
            {
                for (int j = 0; j < this.jGridSize; j++)        //Calculating the initial value for the cell
                {
                    VortexCell vortCell = vortexGrid[i, j];
                    if (!vortCell.exists)
                        continue;

                    vortCell.deltaU = -twoDivBeta * vortCell.localSlope;
                    farInfluence[j] = CalculateFarInfluence(i, j);
                    vortCell.deltaU += farInfluence[j] + CalculateNearInfluence(i, j);
                }
                if (i - 1 > 0)      //If this isn't the frontmost row, we want to go and adjust the velocities of the i-1 row in order to smooth out some oscillations
                {
                    for (int j = 0; j <= this.jGridSize; j++)
                    {
                        VortexCell vortCell;
                        if (j < this.jGridSize)
                        {
                            vortCell = vortexGrid[i - 1, j];
                            if (!vortCell.exists)
                                continue;

                            vortCell.deltaU = vortCell.deltaU * vortCell.forwardVelInfluence + vortexGrid[i, j].deltaU * vortCell.backwardVelInfluence;
                        }
                        int jmin1 = j - 1;
                        if (jmin1 >= 0)
                        {
                            vortCell = vortexGrid[i, jmin1];
                            if (!vortCell.exists)
                                continue;

                            vortCell.deltaU = -twoDivBeta * vortCell.localSlope;
                            vortCell.deltaU += farInfluence[jmin1] + CalculateNearInfluence(i, jmin1);
                        }
                    }
                }
            }
        }

        public void DumpToFile()
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

        //This only accounts for the influence of the closest row ahead, which is the only one affected by the relaxation factor
        private double CalculateNearInfluence(int iIndex, int jIndex)
        {
            double influence = 0;
            int i = iIndex - 1;
            int forwardIndex = iIndex - i;

            if (i < 0)
                return 0;

            for (int j = jIndex - forwardIndex; j <= jIndex + forwardIndex; j++)        //Only get points in a triangle in front of the point in question
            {
                if (j < 0 || j >= jGridSize)
                    continue;

                VortexCell otherCell = vortexGrid[i, j];
                if (!otherCell.exists)
                    continue;

                double curInfluence = influenceGrid[forwardIndex, Math.Abs(jIndex - j)];
                curInfluence *= otherCell.deltaU * otherCell.influenceFactor;

                influence += curInfluence;
            }
            
            return influence;
        }
        
        //Only calculates the effect of influences beyond the row right in front of the cell
        private double CalculateFarInfluence(int iIndex, int jIndex)
        {
            double influence = 0;
            for (int i = iIndex - 2; i >= 0; i--)         //Start from the points right in front of it and work forward to find all influences
            {
                int forwardIndex = iIndex - i;
                for (int j = jIndex - forwardIndex; j <= jIndex + forwardIndex; j++)        //Only get points in a triangle in front of the point in question
                {
                    if (j < 0 || j >= jGridSize)
                        continue;

                    VortexCell otherCell = vortexGrid[i, j];
                    if (!otherCell.exists)
                        continue;

                    double curInfluence = influenceGrid[forwardIndex, Math.Abs(jIndex - j)];
                    curInfluence *= otherCell.deltaU * otherCell.influenceFactor;

                    influence += curInfluence;
                }
            }

            return influence;
        }
    }

    class VortexCell
    {
        public bool exists;
        public double influenceFactor;
        public double forwardVelInfluence;
        public double backwardVelInfluence;
        public double deltaU;
        public double localSlope;
    }
}
