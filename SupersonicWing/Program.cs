using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

            geometry.ScaleWingBased(sideslip * Math.PI / 180, angleOfAttack * Math.PI / 180, beta);

            Console.Write("Choose Horizontal Grid Size: ");
            int gridSize = int.Parse(Console.ReadLine());

            WingSim sim = new WingSim(geometry, gridSize, beta);
            sim.RunSim();

            Console.ReadKey();
        }
    }
}
