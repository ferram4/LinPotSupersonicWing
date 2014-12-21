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
using System.Collections.Generic;
using System.Diagnostics;

namespace SupersonicWing
{
    class WingSim
    {
        WingGeometry wing;

        int jGridSize, iGridSize;
        double gridScale, invGridScale;
        double beta, negTwoDivBeta;

        const double INV_PI = 0.31830988618379067153776752674503;

        VortexCell[][] vortexGrid;
        List<int>[] colBlockLims;

        static double[][] influenceGrid;

        public int countCells = 0;
        public static Stopwatch gridWatch = new Stopwatch();
        public static Stopwatch simWatch = new Stopwatch();


        public void Init(WingGeometry wing, int horizontalGridSize, double beta)
        {
            gridWatch.Start();
            this.wing = wing;
            this.jGridSize = horizontalGridSize;
            gridScale = (double)horizontalGridSize / (wing.maxY - wing.minY);
            invGridScale = 1 / gridScale;
            iGridSize = (int)Math.Round(gridScale * (wing.maxX - wing.minX));

            this.beta = beta;
            negTwoDivBeta = -2 / beta;

            InitializeVortexGrid();
            CheckInfluenceGrid();
            gridWatch.Stop();
        }
        private void CheckInfluenceGrid()
        {
            if (influenceGrid == null)
            {
                influenceGrid = new double[iGridSize][];
                RecalculateInfluenceGrid();
            }
            else
            {
                if (influenceGrid.Length < iGridSize)
                {
                    influenceGrid = new double[Math.Max(iGridSize, influenceGrid.Length)][];
                    RecalculateInfluenceGrid();
                }
            }
        }

        private void RecalculateInfluenceGrid()
        {
            Console.WriteLine("Recalcing Influences");

            for (int i = 0; i < influenceGrid.Length; i++)         //Start from the points right in front of it and work forward to find all influences
            {
                double[] influenceRow = new double[i + 1];
                for (int j = 0; j <= i; j++)        //Only get points in a triangle in front of the point in question
                {
                    double forwardBack = i + 0.5;
                    double forwardBackSqr = forwardBack * forwardBack;

                    double sidePos = j + 0.5;
                    double sideNeg = sidePos - 1;

                    double curInfluence = Math.Sqrt((forwardBackSqr - sideNeg * sideNeg)) / (forwardBack * sideNeg);
                    curInfluence -= Math.Sqrt((forwardBackSqr - sidePos * sidePos)) / (forwardBack * sidePos);
                    curInfluence *= INV_PI;

                    influenceRow[j] = curInfluence;
                }
                influenceGrid[i] = influenceRow;
            }
        }

        private void InitializeVortexGrid()
        {
            if (vortexGrid == null || vortexGrid.Length < iGridSize)
            {
                vortexGrid = new VortexCell[this.iGridSize][];
                colBlockLims = new List<int>[iGridSize];
            }

            bool toggle = false;
            countCells = 0;
            for (int i = 0; i < iGridSize; i++)
            {
                VortexCell[] rowCell = new VortexCell[jGridSize];
                List<int> newList = new List<int>();
                for (int j = 0; j < jGridSize; j++)
                {
                    VortexCell vortCell = new VortexCell();

                    Vector3 test = new Vector3(-(i + 0.5) * invGridScale + wing.maxX, (j + 0.5) * invGridScale + wing.minY, 0);

                    double tmp;
                    Triangle tri;
                    bool exists = wing.TryFindSlopeAtPoint(ref test, out tmp, out tri);
                    if (exists)
                    {
                        countCells++;
                        vortCell.influenceFactor = CalculateInfluenceFactor(ref tri, ref test);
                        vortCell.forwardVelInfluence = 0.5 * (1 + vortCell.influenceFactor / (1 + vortCell.influenceFactor));
                        vortCell.backwardVelInfluence = 0.5 / (1 + vortCell.influenceFactor);
                        vortCell.localSlope = tmp;
                    }
                    if(exists != toggle)
                    {
                        toggle = !toggle;
                        newList.Add(j);
                    }

                    rowCell[j] = vortCell;
                }
                vortexGrid[i] = rowCell;
                colBlockLims[i] = newList;
                toggle = false;
            }
        }

        private double CalculateInfluenceFactor(ref Triangle tri, ref Vector3 test)
        {
            if (wing.TryFindSlopeAtPoint(test + new Vector3(1, 0, 0)))
                return 1;

            double forwardOffset = 0;

            if (tri.p0.x > test.x)
            {
                double x1;
                if ((tri.p0.y > test.y) ^ (tri.p1.y > test.y))
                    x1 = (tri.p0.x - tri.p1.x) / (tri.p0.y - tri.p1.y) * (test.y - tri.p1.y) + tri.p1.x;
                else
                    x1 = (tri.p0.x - tri.p2.x) / (tri.p0.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;

                forwardOffset = x1 - test.x;
            }
            else if (tri.p1.x > test.x)
            {
                double x1;
                if ((tri.p0.y > test.y) ^ (tri.p1.y > test.y))
                    x1 = (tri.p0.x - tri.p1.x) / (tri.p0.y - tri.p1.y) * (test.y - tri.p1.y) + tri.p1.x;
                else
                    x1 = (tri.p1.x - tri.p2.x) / (tri.p1.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;

                forwardOffset = x1 - test.x;
            }
            else
            {
                double x1;
                if ((tri.p0.y > test.y) ^ (tri.p2.y > test.y))
                    x1 = (tri.p0.x - tri.p2.x) / (tri.p0.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;
                else
                    x1 = (tri.p1.x - tri.p2.x) / (tri.p1.y - tri.p2.y) * (test.y - tri.p2.y) + tri.p2.x;

                forwardOffset = x1 - test.x;
            }


            Triangle oldTri = tri;
            if (wing.TryFindSlopeAtPoint(test + new Vector3(forwardOffset + 1e-6, 0, 0), out tri) && !oldTri.Equals(tri))
                return CalculateInfluenceFactor(ref tri, ref test);
            else
                return forwardOffset;
        }

        public void RunSim()
        {
            simWatch.Start();
            double[] farInfluence = new double[jGridSize];      //used to store far influences, since they do not change in the relaxation procedure

            //VortexCell vortCell;
            List<int> rowLims = colBlockLims[0];
            int k = 0;
            int lowerLim, upperLim;

            VortexCell[] rowCell = vortexGrid[0];

            while (k < rowLims.Count)
            {
                lowerLim = rowLims[k];
                k++;
                if (k == rowLims.Count)
                    upperLim = jGridSize;
                else
                {
                    upperLim = rowLims[k];
                    k++;
                }
                for (int j = 0; j < this.jGridSize; j++)
                {
                    VortexCell vortCell = rowCell[j];

                    vortCell.deltaU = negTwoDivBeta * vortCell.localSlope;
                }
            }
            for (int i = 1; i < this.iGridSize; i++)
            {
                rowLims = colBlockLims[i];
                k = 0;
                rowCell = vortexGrid[i];
                VortexCell[] forwardRow = vortexGrid[i - 1];

                while (k < rowLims.Count)
                {
                    lowerLim = rowLims[k];
                    k++;
                    if (k == rowLims.Count)
                        upperLim = jGridSize;
                    else
                    {
                        upperLim = rowLims[k];
                        k++;
                    }
                    for (int j = lowerLim; j < upperLim; j++)        //Calculating the initial value for the cell
                    {

                        VortexCell vortCell = rowCell[j];

                        vortCell.deltaU = negTwoDivBeta * vortCell.localSlope;      //calculate the effects on this panel
                        farInfluence[j] = CalculateFarInfluence(i, j);
                        vortCell.deltaU += farInfluence[j] + CalculateNearInfluence(i, j);
                    }

                    for (int j = lowerLim; j < upperLim; j++)
                    {

                        VortexCell vortCell = forwardRow[j];            //get the one in the row in front...
                        //and adjust its pert vel to remove oscillations
                        vortCell.deltaU = vortCell.deltaU * vortCell.forwardVelInfluence + rowCell[j].deltaU * vortCell.backwardVelInfluence;
                    }
                    for (int j = lowerLim; j < upperLim; j++)
                    {
                        VortexCell vortCell = rowCell[j];

                        vortCell.deltaU = negTwoDivBeta * vortCell.localSlope;
                        vortCell.deltaU += farInfluence[j] + CalculateNearInfluence(i, j);
                    }
                }
            }
            simWatch.Stop();
        }

        public void DumpToFile()
        {
            StreamWriter writer = new StreamWriter(new FileStream("sim.csv", FileMode.Create));

            for (int i = 0; i < this.iGridSize; i++)
            {
                VortexCell[] rowCell = vortexGrid[i];
                for (int j = 0; j < this.jGridSize; j++)
                {
                    writer.Write(rowCell[j].deltaU * -2);
                    writer.Write(",");
                }
                writer.WriteLine();
            }

            writer.Close();
        }

        public void DumpToInfluenceCoeffs()
        {
            StreamWriter writer = new StreamWriter(new FileStream("influence.csv", FileMode.Create));

            for (int i = 0; i < influenceGrid.Length; i++)
            {
                double[] influenceRow = influenceGrid[i];
                for (int j = 0; j < influenceRow.Length; j++)
                {
                    writer.Write(influenceRow[j]);
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
            {
                VortexCell[] rowCell = vortexGrid[i];
                for (int j = 0; j < this.jGridSize; j++)
                {
                    VortexCell vortCell = rowCell[j];
                    //if (!vortCell.exists)
                    //    continue;

                    cL += vortCell.deltaU * -2 * vortCell.influenceFactor;
                    area += vortCell.influenceFactor;
                }
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

            VortexCell[] rowCell = vortexGrid[i];
            double[] influenceRow = influenceGrid[1];
            List<int> rowLims = colBlockLims[i];
            int k = 0;
            int lowerLim, upperLim;


            while (k < rowLims.Count)
            {
                lowerLim = rowLims[k];
                k++;
                if (k == rowLims.Count)
                    upperLim = jGridSize - 1;
                else
                {
                    upperLim = rowLims[k];
                    k++;
                }

                for (int j = Math.Max(jIndex - 1, lowerLim); j <= Math.Min(jIndex + 1, upperLim); j++)        //Only get points in a triangle in front of the point in question
                {
                    VortexCell otherCell = rowCell[j];
                    //if (!otherCell.exists)
                    //    continue;

                    double curInfluence = influenceRow[Math.Abs(jIndex - j)];
                    curInfluence *= otherCell.deltaU * otherCell.influenceFactor;

                    influence += curInfluence;
                }
            }

            return influence;
        }
        
        //Only calculates the effect of influences beyond the row right in front of the cell
        //Currently the most expensive function in this sim
        private double CalculateFarInfluence(int iIndex, int jIndex)
        {
            double influence = 0;
            for (int i = 0; i < iIndex - 1; i++)         //Start from the points right in front of it and work forward to find all influences
            {
                int forwardIndex = iIndex - i;

                VortexCell[] rowCell = vortexGrid[i];
                double[] influenceRow = influenceGrid[forwardIndex];
                List<int> rowLims = colBlockLims[i];
                int k = 0;
                int lowerLim, upperLim;

                while (k < rowLims.Count)
                {
                    lowerLim = rowLims[k];
                    k++;
                    if (k == rowLims.Count)
                        upperLim = jGridSize - 1;
                    else
                    {
                        upperLim = rowLims[k];
                        k++;
                    }

                    for (int j = Math.Max(jIndex - forwardIndex, lowerLim); j <= Math.Min(jIndex + forwardIndex, upperLim); j++)        //Only get points in a triangle in front of the point in question
                    {
                        VortexCell otherCell = rowCell[j];

                        double curInfluence = influenceRow[Math.Abs(jIndex - j)];
                        curInfluence *= otherCell.deltaU * otherCell.influenceFactor;

                        influence += curInfluence;
                    }
                }
            }

            return influence;
        }
    }

    class VortexCell
    {
        public double influenceFactor;
        public double forwardVelInfluence;
        public double backwardVelInfluence;
        public double deltaU;
        public double localSlope;
    }
}
