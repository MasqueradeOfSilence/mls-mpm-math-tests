using System;
using FluidMath;
using NUnit.Framework;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
namespace FluidMathTests
{
	[TestClass]
	public class G2PTest
	{
		private G2P g2p = new G2P();

        [Test]
		public void TestComputeWeightedVelocity()
		{
			double[] neighborCellVelocity = { 0, 1 };
			double weight = 0.16;
			double[] expectedWeightedVelocity = { 0, 0.16 };
			double[] actualWeightedVelocity = g2p.computeWeightedVelocity(neighborCellVelocity, weight);
            NUnit.Framework.Assert.That(actualWeightedVelocity, Is.EqualTo(expectedWeightedVelocity).Within(0.01));
        }

		[Test]
		public void TestComputeTerm()
		{
			// Term should be 2x2! 
			double[] weightedVelocity = { 0, 0.16 };
			double[] distanceFromCellToNeighbor = { -0.7, 0.4 };
			double[,] expectedTerm = { { 0, -0.112 }, { 0, 0.064 } };
			double[,] actualTerm = g2p.computeTerm(weightedVelocity, distanceFromCellToNeighbor);
			NUnit.Framework.Assert.That(actualTerm, Is.EqualTo(expectedTerm).Within(0.01));
		}

		[Test]
		public void TestUpdateB()
		{
			double[,] initialB = { { 0, 0 }, { 0, 0 } };
			double[,] term = { { 0, -0.112 }, { 0, 0.064 } };
			double[,] expectedUpdatedB = { { 0, -0.112 }, { 0, 0.064 } };
			double[,] actualUpdatedB = g2p.computeUpdatedB(initialB, term);
			NUnit.Framework.Assert.That(actualUpdatedB, Is.EqualTo(expectedUpdatedB).Within(0.01));
        }

		[Test]
		public void TestUpdateParticleVelocity()
		{
			double[] initialParticleVelocity = { 2, 2 };
			double[] weightedVelocity = { 0, 0.16 };
			double[] expectedUpdatedParticleVelocity = { 2, 2.16 };
			double[] actualUpdatedParticleVelocity = g2p.computeUpdatedParticleVelocity(initialParticleVelocity, weightedVelocity);
			NUnit.Framework.Assert.That(actualUpdatedParticleVelocity, Is.EqualTo(expectedUpdatedParticleVelocity).Within(0.01));
        }

		[Test]
		public void TestRecomputeCMatrix()
		{
			double[,] B = { { 0, -0.112 }, { 0, 0.064 } };
			double[,] expectedUpdatedC = { { 0, -0.448 }, { 0, 0.256 } };
			double[,] actualUpdatedC = g2p.recomputeCMatrix(B);
			NUnit.Framework.Assert.That(actualUpdatedC, Is.EqualTo(expectedUpdatedC).Within(0.01));
        }

		[Test]
		public void testParticleAdvection()
		{
			double[] initialParticlePosition = { 3.2, 6 };
			double[] particleVelocity = { 2, 2.16 };
			double dt = 0.2;
			double[] expectedParticlePositionAfterAdvection = { 3.6, 6.432 };
			// The real function will probably initially take a Particle class object as its parameter, but will call an operation like this under the hood
			double[] actualParticlePositionAfterAdvection = g2p.advectParticle(initialParticlePosition, particleVelocity, dt);
			NUnit.Framework.Assert.That(actualParticlePositionAfterAdvection, Is.EqualTo(expectedParticlePositionAfterAdvection).Within(0.01));
        }
	}
}

