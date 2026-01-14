#include "pch.h"
#include "CppUnitTest.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

#include "TestHelpers.h"

static bool LogError(const std::string& msg)
{
	std::string errorMessage = "[Error] " + msg;
    Logger::WriteMessage(errorMessage.c_str());
    return false;
}


namespace TestHelpers
{
    bool AreMatricesEqual(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B, double eps)
    {
        if (A.size() != B.size())
        {
            LogError("Matrix row size is different!");
            return false;
        }

        for (size_t i = 0; i < A.size(); i++)
        {
            if (A[i].size() != B[i].size())
            {
                LogError("Matrix column size is different!");
                return false;
            }

            for (size_t j = 0; j < A[i].size(); j++)
            {
				const double diff = std::abs(A[i][j] - B[i][j]);
                if (diff > eps)
                {
                    std::ostringstream msg;
                    msg << "Values in matrix are not equal!\n" <<
                        "Diff in row : " << i << ", column : " << j + "\n" <<
						"First value: " << A[i][j] << ", Second value: " << B[i][j] << ", Diff: " << diff << "\n";

                    LogError(msg.str());
                    return false;
                }
            }
        }

        return true;
    }

    bool AreVectorsEqual(const std::vector<double>& A, const std::vector<double>& B, double eps)
    {
        if (A.size() != B.size())
        {
            LogError("Vector size is different!");
            return false;
        }

        for (size_t i = 0; i < A.size(); i++)
        {
			double diff = std::abs(A[i] - B[i]);
            if (diff > eps)
            {
                std::ostringstream msg;
                msg << "Values in vector are not equal!\n" <<
                    "Diff at index : " + i << "\n" <<
                    "First value: " << A[i] << ", Second value: " << B[i] << ", Diff: " << diff << "\n";

                LogError(msg.str());
                return false;
            }
        }

        return true;
    }
}
