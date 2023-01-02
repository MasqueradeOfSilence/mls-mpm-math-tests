using System;
using FluidMath;
using NUnit.Framework;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace FluidMathTests;

[TestClass]
public class P2G2Test
{
	private P2G2 p2g2 = new P2G2();

	[Test]
	public void TestFindNearestGridCellToParticle()
	{
		double[] particlePosition = { 0.76, 1.61 };
		int[] expectedGridCell = { 1, 2 };
		int[] actualGridCell = p2g2.findNearestGridCellToParticle(particlePosition);
        NUnit.Framework.Assert.That(expectedGridCell, Is.EqualTo(actualGridCell).Within(0.01));
	}
	// next few computations are already tested in the first one.
	// in the actual code, consider keeping these in a different class outside of P2G1
	[Test]
	public void TestUpdateDensity()
	{
		// Not to be confused with particle mass. Remember, data transfer
		double gridCellMass = 3;
		double weight = 1.5376;
		double initialDensity = 0;
		double expectedUpdatedDensity = 4.6128;
		double actualUpdatedDensity = p2g2.computeUpdatedDensity(weight, gridCellMass, initialDensity);
		NUnit.Framework.Assert.That(actualUpdatedDensity, Is.EqualTo(expectedUpdatedDensity).Within(0.01));
	}

	[Test]
	public void TestComputeVolume()
	{
		double particleMass = 2;
		double density = 4.6128;
		double expectedVolume = 9.2256;
		double actualVolume = p2g2.computeVolume(particleMass, density);
		NUnit.Framework.Assert.That(actualVolume, Is.EqualTo(expectedVolume).Within(0.01));
    }

	[Test]
	public void TestComputePressure()
	{
		double expectedPressure = 7.6855;
		double eosStiffness = 10;
		// updated density from previous step
		double density = 4.6128;
		double restDensity = 4;
		double eosPower = 4;
		double actualPressure = p2g2.computePressure(eosStiffness, density, restDensity, eosPower);
		NUnit.Framework.Assert.That(actualPressure, Is.EqualTo(expectedPressure).Within(0.01));
	}

	[Test]
	public void TestCreateStressMatrix()
	{
        double pressure = 7.6855;
        double[,] expectedStress = { { -pressure, 0 }, { 0, -pressure } };
		//NUnit.Framework.TestContext.WriteLine(System.String.Join(" ", expectedStress.Cast<double>()));
		double[,] actualStress = p2g2.createStressMatrix(pressure);
		NUnit.Framework.Assert.That(actualStress, Is.EqualTo(expectedStress).Within(0.01));
	}

	[Test]
	public void TestInitializeStrainMatrix()
	{
		double[,] C = { { 2, 3 }, { 4, 1 } };
		// It's just initialized to the particle's C value at first.
		// But we're not building the particle right now since we're just isolating the math here.
		double[,] expectedStrain = { { 2, 3 }, { 4, 1 } };
		double[,] actualStrain = p2g2.initializeStrainMatrix(C);
		NUnit.Framework.Assert.That(actualStrain, Is.EqualTo(expectedStrain));
    }

	[Test]
	public void TestComputeTrace()
	{
		double[,] initialStrain = { { 2, 3 }, { 4, 1 } };
		double expectedTrace = 7;
		double actualTrace = p2g2.computeTrace(initialStrain);
		NUnit.Framework.Assert.That(actualTrace, Is.EqualTo(expectedTrace));
    }

	[Test]
	public void TestUpdateStrain()
	{
		double[,] expectedUpdatedStrain = { { 2, 7 }, { 7, 1 } };
        double[,] initialStrain = { { 2, 3 }, { 4, 1 } };
        double trace = 7;
		double[,] actualUpdatedStrain = p2g2.updateStrain(initialStrain, trace);
		NUnit.Framework.Assert.That(actualUpdatedStrain, Is.EqualTo(expectedUpdatedStrain));
	}

	[Test]
	public void TestComputeViscosity()
	{
		double dynamicViscosity = 0.1;
		double[,] strain = { { 2, 7 }, { 7, 1 } };
		double[,] expectedViscosity = { { 0.2, 0.7 }, { 0.7, 0.1 } };
		double[,] actualViscosity = p2g2.computeVisoscity(strain, dynamicViscosity);
		NUnit.Framework.Assert.That(actualViscosity, Is.EqualTo(expectedViscosity).Within(0.01));
    }

	[Test]
	public void TestUpdateStress()
	{
		double[,] viscosity = { { 0.2, 0.7 }, { 0.7, 0.1 } };
        double[,] initialStress = { { -7.6855, 0 }, { 0, -7.6855 } };
        // stress += viscosity, with viscosity also being a 2x2 matrix
        double[,] expectedUpdatedStress = { { -7.4855 , 0.7 }, { 0.7, -7.5855 } };
		double[,] actualUpdatedStress = p2g2.updateStress(initialStress, viscosity);
		NUnit.Framework.Assert.That(actualUpdatedStress, Is.EqualTo(expectedUpdatedStress).Within(0.01));
    }

	[Test]
	public void TestComputeTerm0()
	{
		// -volume * 4 * stress * dt
		double[,] stress = { { -7.4855, 0.7 }, { 0.7, -7.5855 } };
		double volume = 1;
		double dt = 0.2;
		// Based on wolfram alpha calculator
		double[,] expectedTerm0 = { { 5.9884, -0.56 }, { -0.56, 6.0684 } };
		double[,] actualTerm0 = p2g2.computeEquation16Term0(stress, volume, dt);
		NUnit.Framework.Assert.That(actualTerm0, Is.EqualTo(expectedTerm0).Within(0.01));
    }

	[Test]
	public void TestComputeDistanceFromCellToNeighbor()
	{
		double[] neighborCellPosition = { 0, 1 };
		double[] currentCellPosition = { 1, 2 };
		double[] expectedDistance = { -0.5, -0.5 };
		double[] actualDistance = p2g2.computeDistanceFromCellToNeighbor(neighborCellPosition, currentCellPosition);
		NUnit.Framework.Assert.That(actualDistance, Is.EqualTo(expectedDistance).Within(0.01));
	}

	[Test]
	public void TestComputeMomentum()
	{
		// fused force + momentum update
		double[,] equation16Term0 = { { 5.9884, -0.56 }, { -0.56, 6.0684 } };
		double weight = 1.5376;
		double[,] distanceFromCellToNeighbor = { { -0.5 }, { -0.5 } };
		double[,] expectedMomentum = { { -4.17335 }, { -4.23486 } };
		double[,] actualMomentum = p2g2.computeMomentum(equation16Term0, weight, distanceFromCellToNeighbor);
		NUnit.Framework.Assert.That(actualMomentum, Is.EqualTo(expectedMomentum).Within(0.01));
    }

	[Test]
	public void TestUpdateCellVelocity()
	{
		// initial state has this velocity at 0
		double initialCellVelocity = 0;
        double[,] momentum = { { -4.17335 }, { -4.23486 } };
		double[,] expectedUpdatedCellVelocity = { { -4.17335 }, { -4.23486 } };
		double[,] actualUpdatedCellVelocity = p2g2.updateCellVelocity(momentum, initialCellVelocity);
		NUnit.Framework.Assert.That(actualUpdatedCellVelocity, Is.EqualTo(expectedUpdatedCellVelocity).Within(0.01));
    }

}

