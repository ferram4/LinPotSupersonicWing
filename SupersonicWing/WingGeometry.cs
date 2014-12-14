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
    class WingGeometry
    {
        Triangle[] wingTriangles;   //triangles making up wing geometry

        double[] dz_dx;             //slopes associated with each triangle
        double[] dz_dy;

        Triangle[] transformedTriangles;    //triangles rotated and scaled to account for mach number and flow direction
        double[] currentTriAngle;           //AoA and actual slope combined

        public double maxX, minX, maxY, minY;      //extents of the transformed wing for later gridding

        public WingGeometry(Triangle[] unflattenedTris)
        {
            dz_dx = new double[unflattenedTris.Length];         //create new arrays for the slopes
            dz_dy = new double[unflattenedTris.Length];

            //Calculate slopes from tris, then flatten them out

            wingTriangles = unflattenedTris;

            transformedTriangles = new Triangle[wingTriangles.Length];      //Initialize these arrays now
            currentTriAngle = new double[unflattenedTris.Length];           

        }

        public void ScaleWingBased(double sideslipAngleRad, double angleOfAttackRad, double beta)
        {
            double cosSideslip = Math.Cos(sideslipAngleRad);
            double sinSideslip = Math.Sin(sideslipAngleRad);
            double tanAoA = Math.Tan(angleOfAttackRad);

            maxX = maxY = minX = minY = 0;

            for (int i = 0; i < wingTriangles.Length; i++)
            {
                Triangle tri = wingTriangles[i];
                for (int j = 0; j < 3; j++)
                {
                    Vector3 tmp;
                    tmp.x = tri[j].x * cosSideslip + tri[j].y * sinSideslip;    //rotate the triangles
                    tmp.y = -tri[j].x * sinSideslip + tri[j].y * cosSideslip;
                    tmp.z = 0;

                    tmp.y *= beta;      //and squish in the appropriate direction

                    maxX = Math.Max(maxX, tmp.x);
                    minX = Math.Min(minX, tmp.x);       //Update the extents of the wing

                    maxY = Math.Max(maxY, tmp.y);
                    minY = Math.Min(minY, tmp.y);

                    tri[j] = tmp;
                }
                transformedTriangles[i] = tri;

                currentTriAngle[i] = dz_dx[i] * cosSideslip + dz_dy[i] * sinSideslip + tanAoA;
            }
        }

        //Returns true if a triangle is found there
        public bool TryFindSlopeAtPoint(Vector3 test, out double slope, out Triangle tri)
        {
            slope = 0;
            tri = new Triangle();
            for(int i = 0; i < transformedTriangles.Length; i++)
            {
                if(transformedTriangles[i].Contains(test))
                {
                    slope = currentTriAngle[i];
                    tri = transformedTriangles[i];
                    return true;
                }
            }
            return false;
        }
        public bool TryFindSlopeAtPoint(Vector3 test, out Triangle tri)
        {
            tri = new Triangle();
            for (int i = 0; i < transformedTriangles.Length; i++)
            {
                if (transformedTriangles[i].Contains(test))
                {
                    tri = transformedTriangles[i];
                    return true;
                }
            }
            return false;
        }
        public bool TryFindSlopeAtPoint(Vector3 test)
        {
            for (int i = 0; i < transformedTriangles.Length; i++)
            {
                if (transformedTriangles[i].Contains(test))
                {
                    return true;
                }
            }
            return false;
        }
    }

    struct Triangle
    {
        public Vector3 p0, p1, p2;

        public Vector3 this[int index]
        {
            get
            {
                if (index == 0)
                    return p0;
                else if (index == 1)
                    return p1;
                else
                    return p2;
            }
            set
            {
                if (index == 0)
                    p0 = value;
                else if (index == 1)
                    p1 = value;
                else
                    p2 = value;

            }
        }

        public bool Contains(Vector3 test)
        {
            Vector3 v0, v1, v2;
            v0 = p2 - p0;
            v1 = p1 - p0;
            v2 = test - p0;

            double dot00, dot01, dot11, dot02, dot12;
            dot00 = Vector3.Dot(v0, v0);
            dot01 = Vector3.Dot(v0, v1);
            dot02 = Vector3.Dot(v0, v2);
            dot11 = Vector3.Dot(v1, v1);
            dot12 = Vector3.Dot(v1, v2);

            double recipDenom = 1 / (dot00 * dot11 - dot01 * dot01);
            double u = (dot11 * dot02 - dot01 * dot12) * recipDenom;
            double v = (dot00 * dot12 - dot01 * dot02) * recipDenom;


            return (u >= 0) && (v >= 0) && (u + v < 1);
        }
    }

    struct Vector3
    {
        public Vector3(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public double x, y, z;
        public double this[int index]
        {
            get
            {
                if (index == 0)
                    return x;
                else if (index == 1)
                    return y;
                else
                    return z;
            }
            set
            {
                if (index == 0)
                    x = value;
                else if (index == 1)
                    y = value;
                else
                    z = value;

            }
        }

        public static Vector3 operator + (Vector3 a, Vector3 b)
        {
            return new Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static Vector3 operator - (Vector3 a, Vector3 b)
        {
            return new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static double Dot(Vector3 a, Vector3 b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
    }
}
