using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;
using MathNet.Numerics.Statistics.Mcmc;
using System.Diagnostics;
using MathNet.Spatial.Euclidean;
using System.IO;

namespace FluidMath
{
    // See https://aka.ms/new-console-template for more information
    //Console.WriteLine("Hello, World!");

    public class P2G1
    {
        // Note: There will probably be a class overlaid on top of the double array. 
        public int[] particlePositionToCellPosition(double[] particlePosition)
        {
            return Array.ConvertAll<double, int>(particlePosition, x => (int)x);
        }

        public double[] computeDistanceFromParticleToCell(double[] particlePosition, int[] correspondingCellPosition)
        {
            double[] cellPosition = Array.ConvertAll<int, double>(correspondingCellPosition, x => x);
            //return (particlePosition - cellPosition) - 0.5;
            double[] distance = { particlePosition[0] - cellPosition[0] - 0.5, particlePosition[1] - cellPosition[1] - 0.5 };
            return distance;
        }

        public double[] computeWeight0(double[] distanceFromParticleToCell)
        {
            double x = 0.5 * Math.Pow((0.5 - distanceFromParticleToCell[0]), 2);
            double y = 0.5 * Math.Pow((0.5 - distanceFromParticleToCell[1]), 2);
            double[] weight0 = { x, y };
            return weight0;
        }

        public double[] computeWeight1(double[] distanceFromParticleToCell)
        {
            double x = 0.75 - Math.Pow(distanceFromParticleToCell[0], 2);
            double y = 0.75 - Math.Pow(distanceFromParticleToCell[1], 2);
            double[] weight1 = { x, y };
            return weight1;
        }

        public double[] computeWeight2(double[] distanceFromParticleToCell)
        {
            double x = 0.5 * Math.Pow(0.5 + distanceFromParticleToCell[0], 2);
            double y = 0.5 * Math.Pow(0.5 + distanceFromParticleToCell[1], 2);
            double[] weight2 = { x, y };
            return weight2;
        }

        public double[][] computeWeights(double[] distanceFromParticleToCell)
        {
            double[] weight0 = computeWeight0(distanceFromParticleToCell);
            double[] weight1 = computeWeight1(distanceFromParticleToCell);
            double[] weight2 = computeWeight2(distanceFromParticleToCell);
            double[][] weights = { weight0, weight1, weight2 };
            return weights;
        }

        public double computeWeight(double[][] weights, int nx, int ny)
        {
            return weights[nx][0] * weights[ny][1];
        }

        public double[] computeNeighborPosition(int[] correspondingCellPosition, int nx, int ny)
        {
            double x = correspondingCellPosition[0] + nx - 1;
            double y = correspondingCellPosition[1] + ny - 1;
            double[] neighborPosition = { x, y };
            return neighborPosition;
        }

        public double[] computeDistanceFromCurrentParticleToCurrentNeighbor(int[] currentNeighborPosition, double[] particlePosition)
        {
            double x = (currentNeighborPosition[0] - particlePosition[0]) + 0.5;
            double y = (currentNeighborPosition[1] - particlePosition[1]) + 0.5;
            double[] distanceFromCurrentParticleToCurrentNeighbor = { x, y };
            return distanceFromCurrentParticleToCurrentNeighbor;
        }

        public double[,] computeQ(double[,] C, double[,] distanceFromCurrentParticleToCurrentNeighbor)
        {
            Matrix<double> a = DenseMatrix.OfArray(C);
            Matrix<double> b = DenseMatrix.OfArray(distanceFromCurrentParticleToCurrentNeighbor);
            Matrix<double> Q = a * b;
            return Q.ToArray();
        }

        public double computeMassContribution(double weight, double particleMass)
        {
            return weight * particleMass;
        }

        public double recomputeCellMassAndReturnIt(double initialCellMass, double massContribution)
        {
            return initialCellMass + massContribution;
        }

        public double[] recomputeCellVelocityAndReturnIt(double massContribution, double[] particleVelocity, double[,] Q)
        {
            double x = massContribution * (particleVelocity[0] + Q[0, 0]);
            double y = massContribution * (particleVelocity[1] + Q[1, 0]);
            double[] newCellVelocity = { particleVelocity[0] + x, particleVelocity[1] + y };
            return newCellVelocity;
        }

        static void Main()
        {

        }
    }

    public class P2G2
    {
        public int[] findNearestGridCellToParticle(double[] particlePosition)
        {
            // Do NOT cast! C# will simply truncate the decimal!
            int x = Convert.ToInt32(particlePosition[0]);
            int y = Convert.ToInt32(particlePosition[1]);
            int[] nearestGridCellToParticle = { x, y };
            return nearestGridCellToParticle;
        }

        public double computeUpdatedDensity(double weight, double gridCellMass, double initialDensity)
        {
            return initialDensity + (gridCellMass * weight);
        }

        public double computeVolume(double particleMass, double density)
        {
            return particleMass * density;
        }

        public double computePressure(double eosStiffness, double density, double restDensity, double eosPower)
        {
            // Note: eosStiffness is applied to the term after it's raised to the power and 1 is subtracted from it.
            return Math.Max(-0.1, eosStiffness * (Math.Pow((density / restDensity), eosPower) - 1));
        }

        public double[,] createStressMatrix(double pressure)
        {
            double[,] stressMatrix = { { -pressure, 0 }, { 0, -pressure } };
            return stressMatrix;
        }

        public double[,] initializeStrainMatrix(double[,] C)
        {
            double[,] strainMatrix = C;
            return strainMatrix;
        }

        public double computeTrace(double[,] strain)
        {
            return strain[1, 0] + strain[0, 1];
        }

        public double[,] updateStrain(double[,] initialStrain, double trace)
        {
            initialStrain[0, 1] = trace;
            initialStrain[1, 0] = trace;
            return initialStrain;
        }

        public double[,] computeVisoscity(double[,] strain, double dynamicViscosity)
        {
            Matrix<double> a = DenseMatrix.OfArray(strain);
            a *= dynamicViscosity;
            double[,] viscosity = a.ToArray();
            return viscosity;
        }

        public double[,] updateStress(double[,] initialStress, double[,] viscosity)
        {
            Matrix<double> a = DenseMatrix.OfArray(initialStress);
            Matrix<double> b = DenseMatrix.OfArray(viscosity);
            Matrix<double> updatedStress = a + b;
            return updatedStress.ToArray();
        }

        public double[,] computeEquation16Term0(double[,] stress, double volume, double dt)
        {
            Matrix<double> stressMatrix = DenseMatrix.OfArray(stress);
            Matrix<double> term0 = -volume * 4 * stressMatrix * dt;
            return term0.ToArray();
        }

        public double[] computeDistanceFromCellToNeighbor(double[] neighborCellPosition, double[] currentCellPosition)
        {
            double x = (neighborCellPosition[0] - currentCellPosition[0]) + 0.5;
            double y = (neighborCellPosition[1] - currentCellPosition[1]) + 0.5;
            double[] distanceFromCellToNeighbor = { x, y };
            return distanceFromCellToNeighbor;
        }

        public double[,] computeMomentum(double[,] equation16Term0, double weight, double[,] distanceFromCellToNeighbor)
        {
            Matrix<double> equation16Term0Matrix = DenseMatrix.OfArray(equation16Term0);
            Matrix<double> firstFactorMatrix = equation16Term0Matrix * weight;
            Matrix<double> secondFactorMatrix = DenseMatrix.OfArray(distanceFromCellToNeighbor);
            // Unity will have its own math library
            // Matrix multiplication here
            return (firstFactorMatrix * secondFactorMatrix).ToArray();
        }

        public double[,] updateCellVelocity(double[,] momentum, double initialCellVelocity)
        {
            return (DenseMatrix.OfArray(momentum) + initialCellVelocity).ToArray();
        }
    }

}