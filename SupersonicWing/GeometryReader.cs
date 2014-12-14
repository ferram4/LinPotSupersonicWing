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
