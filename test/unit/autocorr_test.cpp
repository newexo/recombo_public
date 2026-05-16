#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "autocorr.h"

using namespace std;

// Test data from legacy autocorr tests
const double autocorrTestData[] = {
    1.2190079223,
    -8.9262668602,
    7.7765452396,
    -1.3448207639,
    -7.7495624404,
    -1.7062141187,
    -7.2927681822,
    3.8728410844,
    -6.2263957225,
    -6.1935190763,
    7.8775635548,
    4.2208929267,
    -1.3419686165,
    8.3451695926,
    6.3461351953,
    3.8499901164,
    7.0610383339,
    -9.9842282571,
    0.0627100654,
    -8.4068057034,
    -1.7178643402,
    -4.5115231071,
    -4.5541757066,
    0.8381974231,
    2.4940358102,
    -9.5345366467,
    -1.4188661613,
    -1.5643732436,
    2.2443082184,
    0.4945876356
};

const int autocorrTestDataSize = 30;
const double EPSILON = 1.0e-10;

class AutocorrTest : public ::testing::Test
{
protected:
    autocorr ac;
    vector<double> data;

    void SetUp() override
    {
        // Load test data into vector
        for (int i = 0; i < autocorrTestDataSize; i++)
            data.push_back(autocorrTestData[i]);
    }
};

TEST_F(AutocorrTest, ComputeMeanFromArray)
{
    double mean = ac.computeMean(autocorrTestData, autocorrTestDataSize);
    double expectedMean = -0.8590288609;
    EXPECT_NEAR(expectedMean, mean, EPSILON);
}

TEST_F(AutocorrTest, ComputeMeanFromVector)
{
    double mean = ac.computeMean(data);
    double expectedMean = -0.8590288609;
    EXPECT_NEAR(expectedMean, mean, EPSILON);
}

TEST_F(AutocorrTest, ComputeMeanBothMethodsMatch)
{
    double meanArray = ac.computeMean(autocorrTestData, autocorrTestDataSize);
    double meanVector = ac.computeMean(data);
    EXPECT_NEAR(meanArray, meanVector, EPSILON);
}

TEST_F(AutocorrTest, ComputeMeanSquareFromArray)
{
    double meanSq = ac.computeMeanSquare(autocorrTestData, autocorrTestDataSize);
    double expectedMeanSq = 30.9991617759;
    EXPECT_NEAR(expectedMeanSq, meanSq, EPSILON);
}

TEST_F(AutocorrTest, ComputeMeanSquareFromVector)
{
    double meanSq = ac.computeMeanSquare(data);
    double expectedMeanSq = 30.9991617759;
    EXPECT_NEAR(expectedMeanSq, meanSq, EPSILON);
}

TEST_F(AutocorrTest, ComputeMeanSquareBothMethodsMatch)
{
    double meanSqArray = ac.computeMeanSquare(autocorrTestData, autocorrTestDataSize);
    double meanSqVector = ac.computeMeanSquare(data);
    EXPECT_NEAR(meanSqArray, meanSqVector, EPSILON);
}

TEST_F(AutocorrTest, ComputeVarianceFromArray)
{
    double var = ac.computeVariance(autocorrTestData, autocorrTestDataSize);
    double expectedVariance = 30.261231192;
    EXPECT_NEAR(expectedVariance, var, EPSILON);
}

TEST_F(AutocorrTest, ComputeVarianceFromVector)
{
    double var = ac.computeVariance(data);
    double expectedVariance = 30.261231192;
    EXPECT_NEAR(expectedVariance, var, EPSILON);
}

TEST_F(AutocorrTest, ComputeVarianceBothMethodsMatch)
{
    double varArray = ac.computeVariance(autocorrTestData, autocorrTestDataSize);
    double varVector = ac.computeVariance(data);
    EXPECT_NEAR(varArray, varVector, EPSILON);
}

TEST_F(AutocorrTest, ComputeVarianceMatchesMeanSquareMinusMeanSquared)
{
    double mean = ac.computeMean(autocorrTestData, autocorrTestDataSize);
    double meanSq = ac.computeMeanSquare(autocorrTestData, autocorrTestDataSize);
    double var = ac.computeVariance(autocorrTestData, autocorrTestDataSize);

    double expectedVar = meanSq - mean * mean;
    EXPECT_NEAR(expectedVar, var, EPSILON);
}

TEST_F(AutocorrTest, AutocorrInfoConstruction)
{
    autocorrInfo aci(1.5, 2.0, 3.5, 4.0, 5);
    EXPECT_DOUBLE_EQ(1.5, aci.autocorr);
    EXPECT_DOUBLE_EQ(2.0, aci.error_autocorr);
    EXPECT_DOUBLE_EQ(3.5, aci.mean);
    EXPECT_DOUBLE_EQ(4.0, aci.error_mean);
    EXPECT_EQ(5, aci.window_size);
}

TEST_F(AutocorrTest, AutocorrInfoDefaultConstruction)
{
    autocorrInfo aci;
    EXPECT_FALSE(aci.computationError);
}

TEST_F(AutocorrTest, AutocorrInfoAssignment)
{
    autocorrInfo aci1(1.5, 2.0, 3.5, 4.0, 5);
    autocorrInfo aci2;
    aci2 = aci1;
    EXPECT_EQ(aci1, aci2);
}

TEST_F(AutocorrTest, AutocorrInfoEquality)
{
    autocorrInfo aci1(1.5, 2.0, 3.5, 4.0, 5);
    autocorrInfo aci2(1.5, 2.0, 3.5, 4.0, 5);
    EXPECT_EQ(aci1, aci2);
}

TEST_F(AutocorrTest, AutocorrInfoInequality)
{
    autocorrInfo aci1(1.5, 2.0, 3.5, 4.0, 5);
    autocorrInfo aci2(1.0, 2.0, 3.5, 4.0, 5);
    EXPECT_NE(aci1, aci2);
}

TEST_F(AutocorrTest, AutocorrelationFromArray)
{
    autocorrInfo infoArray = ac.autocorrelation(autocorrTestData, autocorrTestDataSize, false);
    EXPECT_FALSE(infoArray.computationError);
}

TEST_F(AutocorrTest, AutocorrelationFromVector)
{
    autocorrInfo infoVector = ac.autocorrelation(data, false);
    EXPECT_FALSE(infoVector.computationError);
}

TEST_F(AutocorrTest, AutocorrelationBothMethodsMatch)
{
    autocorrInfo infoArray = ac.autocorrelation(autocorrTestData, autocorrTestDataSize, false);
    autocorrInfo infoVector = ac.autocorrelation(data, false);
    EXPECT_EQ(infoArray, infoVector);
}
