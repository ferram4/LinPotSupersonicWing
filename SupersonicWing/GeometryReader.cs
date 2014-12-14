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
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace SupersonicWing
{
    static class GeometryReader
    {
        public static WingGeometry ReadWingGeometryFromFile(string file)
        {
            StreamReader reader = new StreamReader(new FileStream(file, FileMode.Open));
            List<Triangle> triList = new List<Triangle>();

            int i = 0;
            while(!reader.EndOfStream)
            {
                if(reader.ReadLine().ToLowerInvariant() == "triangle")
                {
                    Triangle tri = new Triangle();
                    string[] pointStr;
                    for (int j = 0; j < 3; j++)
                    {
                        pointStr = reader.ReadLine().Split(new char[] { ',', ' ', ';' });

                        if (pointStr.Length > 3)
                            throw new Exception("Triangle " + i + " has an improperly defined point");
                        Vector3 point = new Vector3();

                        for(int k = 0; k < pointStr.Length; k++)
                        {
                            double tmp;
                            if (!double.TryParse(pointStr[k], out tmp))
                                throw new Exception("Point in triangle " + i + " has a point definition that is not parsable as a double");

                            point[k] = tmp;
                        }
                        tri[j] = point;
                    }
                    triList.Add(tri);
                }
            }
            reader.Close();
            return new WingGeometry(triList.ToArray());
        }
    }
}
