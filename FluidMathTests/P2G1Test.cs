namespace FluidMathTests;
using FluidMath;
using NUnit.Framework;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

[TestClass]
public class P2G1Test
{
    private P2G1 p2g1 = new P2G1();
    private TestContext testContextInstance;
    // Use TextContext.WriteLine to print to unit test results
    public TestContext TestContext
    {
        get { return testContextInstance; }
        set { testContextInstance = value; }
    }

    [Test]
    public void TestCastParticlePositiontoCell()
    {
        double[] position = { 3.2, 6 };
        int[] expectedCorrespondingCellPosition = { 3, 6 };
        int[] correspondingCellPosition = p2g1.particlePositionToCellPosition(position);
        CollectionAssert.AreEqual(correspondingCellPosition, expectedCorrespondingCellPosition);
        // sanity check
        int[] sanityCheck1 = { 0, 2 };
        int[] sanityCheck2 = { 1, 2 };
        CollectionAssert.AreNotEqual(sanityCheck1, sanityCheck2);
    }

    [Test]
    public void TestComputeDistanceFromParticleToCell()
    {
        double[] position = { 3.2, 6 };
        int[] correspondingCellPosition = { 3, 6 };
        double[] expectedDistance = { -0.3, -0.5 };
        double[] distanceFromParticleToCell = p2g1.computeDistanceFromParticleToCell(position, correspondingCellPosition);
        Assert.That(expectedDistance, Is.EqualTo(distanceFromParticleToCell).Within(0.01));
    }

    [Test]
    public void TestComputeWeight0()
    {
        double[] expectedWeight0 = { 0.32, 0.5 };
        double[] distanceFromParticleToCell = { -0.3, -0.5 };
        double[] actualWeight0 = p2g1.computeWeight0(distanceFromParticleToCell);
        Assert.That(expectedWeight0, Is.EqualTo(actualWeight0).Within(0.01));
    }

    [Test]
    public void TestComputeWeight1()
    {
        double[] expectedWeight1 = { 0.66, 0.5 };
        double[] distanceFromParticleToCell = { -0.3, -0.5 };
        double[] actualWeight1 = p2g1.computeWeight1(distanceFromParticleToCell);
        Assert.That(expectedWeight1, Is.EqualTo(actualWeight1).Within(0.01));
    }

    [Test]
    public void TestComputeWeight2()
    {
        double[] expectedWeight2 = { 0.02, 0 };
        double[] distanceFromParticleToCell = { -0.3, -0.5 };
        double[] actualWeight2 = p2g1.computeWeight2(distanceFromParticleToCell);
        Assert.That(expectedWeight2, Is.EqualTo(actualWeight2).Within(0.01));
    }

    [Test]
    public void TestComputeAllWeights()
    {
        double[] weight0 = { 0.32, 0.5 };
        double[] weight1 = { 0.66, 0.5 };
        double[] weight2 = { 0.02, 0 };
        double[][] expectedWeights = { weight0, weight1, weight2 };
        double[] distanceFromParticleToCell = { -0.3, -0.5 };
        double[][] actualWeights = p2g1.computeWeights(distanceFromParticleToCell);
        Assert.That(expectedWeights, Is.EqualTo(actualWeights).Within(0.01));
    }

    [Test]
    public void TestComputeWeight()
    {
        double expectedWeight = 0.16;
        double[] weight0 = { 0.32, 0.5 };
        double[] weight1 = { 0.66, 0.5 };
        double[] weight2 = { 0.02, 0 };
        double[][] weights = { weight0, weight1, weight2 };
        int nx = 0;
        int ny = 0;
        double actualWeight = p2g1.computeWeight(weights, nx, ny);
        Assert.That(expectedWeight, Is.EqualTo(actualWeight).Within(0.01));
    }

    [Test]
    public void TestComputeCurrentNeighborPosition()
    {
        int nx = 0;
        int ny = 0;
        int[] correspondingCellPosition = { 3, 6 };
        double[] expectedNeighborPosition = { 2, 5 };
        double[] actualNeighborPosition = p2g1.computeNeighborPosition(correspondingCellPosition, nx, ny);
        Assert.That(expectedNeighborPosition, Is.EqualTo(actualNeighborPosition).Within(0.01));
    }

    [Test]
    public void TestDistanceFromCurrentParticleToCurrentNeighbor()
    {
        double[] expectedDistance = { -0.7, -0.5 };
        int[] neighborPosition = { 2, 5 };
        double[] particlePosition = { 3.2, 6 };
        double[] actualDistance = p2g1.computeDistanceFromCurrentParticleToCurrentNeighbor(neighborPosition, particlePosition);
        Assert.That(expectedDistance, Is.EqualTo(actualDistance).Within(0.01));
    }

    [Test]
    public void TestComputeQ()
    {
        // note the slight difference in format, may need to be accounted for
        double[,] C = { { 0.1, 1 }, { 0.2, 2 } };
        double[,] distanceFromCurrentParticleToCurrentNeighbor = { { -0.7 }, { -0.5 } };
        double[,] Q = p2g1.computeQ(C, distanceFromCurrentParticleToCurrentNeighbor);
        double[,] expectedQ = { { -0.57 }, { -1.14 } };
        Assert.That(Q, Is.EqualTo(expectedQ).Within(0.01));
    }

    [Test]
    public void TestComputeMassContribution()
    {
        double expectedMassContribution = 0.16;
        double particleMass = 1;
        double weight = 0.16;
        double actualMassContribution = p2g1.computeMassContribution(weight, particleMass);
        Assert.That(expectedMassContribution, Is.EqualTo(actualMassContribution));
    }

    [Test]
    public void TestCellMassUpdate()
    {
        double initialCellMass = 1;
        double massContribution = 0.16;
        double expectedUpdatedCellMass = 1.16;
        // recompute cell mass because it just returns the updated value, which we would then update inside of a hypothetical cell object
        double actualUpdatedCellMass = p2g1.recomputeCellMassAndReturnIt(initialCellMass, massContribution);
        Assert.That(expectedUpdatedCellMass, Is.EqualTo(actualUpdatedCellMass).Within(0.01));
    }

    [Test]
    public void TestCellVelocityUpdate()
    {
        // This will be for the neighbor cell in the big picture. 
        // How do we update cell velocity? By addding momentum to it.
        // Momentum is computed as mass contribution * (particle velocity + Q)
        double massContribution = 0.16;
        double[] particleVelocity = { 2, 2 };
        double[] expectedNewCellVelocity = { 2.2288, 2.1376 };
        double[,] Q = { { -0.57 }, { -1.14 } };
        double[] actualNewCellVelocity = p2g1.recomputeCellVelocityAndReturnIt(massContribution, particleVelocity, Q);
        Assert.That(expectedNewCellVelocity, Is.EqualTo(actualNewCellVelocity).Within(0.01));
    }

    // Make sure to also test the synthesis of the two, when you put the whole code together.
    // And commit this to GitHub
}
