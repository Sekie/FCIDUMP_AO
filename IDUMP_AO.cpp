#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip> 

int main()
{
    std::string BaseName = "h2h2h2";
    int NumAO = 6;
    int NumOps = 8;
    int NumFiles = 1;
    for (int r = 0; r < NumFiles; r++)
    {
        std::ifstream idump(BaseName + "_" + std::to_string(r));
        std::map< std::string, double > Integrals;

        int tmpInt1, tmpInt2, tmpInt3, tmpInt4;
        double tmpDouble;

        std::vector< std::vector< int > > TwoElectronIndexes;
        std::vector< std::vector< int > > OneElectronIndexes;
        std::vector< std::string > OptionLines;
        for (int s = 0; s < NumOps; s++)
        {
            std::string tmpstring;
            getline(idump, tmpstring);
            OptionLines.push_back(tmpstring);
        }

        while(!idump.eof())
        {
            /* We have to include all 8-fold permuation symmetries. This holds each integral in chemistry notation. We represent
            (ij|kl) as "i j k l". h_ij is "i j 0 0", as given in QChem. */
            idump >> tmpDouble >> tmpInt1 >> tmpInt2 >> tmpInt3 >> tmpInt4;
            // IntegralsFile >> tmpInt1 >> tmpInt2 >> tmpInt3 >> tmpInt4 >> tmpDouble;
            Integrals[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
            Integrals[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
            Integrals[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
            Integrals[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
            Integrals[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
            Integrals[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
            Integrals[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
            Integrals[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;

            std::vector< int > tmpVec;
            tmpVec.push_back(tmpInt1); tmpVec.push_back(tmpInt2); tmpVec.push_back(tmpInt3); tmpVec.push_back(tmpInt4);
            if(tmpInt1 == 0)
            {
                continue;
            }
            else if(tmpInt2 == 0)
            {
                continue;
            }
            else if(tmpInt4 == 0)
            {
                OneElectronIndexes.push_back(tmpVec);
            }
            else
            {
                TwoElectronIndexes.push_back(tmpVec);
            }
        }
        
        Eigen::MatrixXd CoeffMatrix(NumAO, NumAO);
        Eigen::MatrixXd InvCoeff(NumAO, NumAO);
        Eigen::MatrixXd OverlapMatrix(NumAO, NumAO);

        // InvCoeff(0, 0) = 1 / sqrt(2);
        // InvCoeff(1, 0) = 1 / sqrt(2);
        // InvCoeff(0, 1) = 1 / sqrt(2);
        // InvCoeff(1, 1) = -1 / sqrt(2);

        // CoeffMatrix << 
        // 1 / sqrt(2) , 1 / sqrt(2), 
        // 1 / sqrt(2) , -1 / sqrt(2);

        // H2O
    //     CoeffMatrix << 
    //     0.99422,   0.23375 , -0.00000  , 0.10406 , -0.00000  , 0.12588, 0.00000 ,
    //     0.02586 , -0.84430 ,  0.00000,  -0.53830,   0.00000,  -0.82084, 0,
    //     0.00000  , 0.00000  , 0.61268  , 0.00000 , -0.00000 , -0.00000, 0.96005,
    //     -0.00000 , -0.00000 ,  0.00000 ,  0.00000  , 1.00000 ,  0.00000, 0,
    //     -0.00417  , 0.12297 ,  0.00000,  -0.75600 ,  0.00000  , 0.76353, 0,
    // -0.00559  ,-0.15563 , -0.44920 ,  0.29494  ,-0.00000  , 0.76945, 0.81497,
    //     -0.00559 , -0.15563  , 0.44920 ,  0.29494 , -0.00000  , 0.76945, -0.81497;


    // CoeffMatrix << 
    // 1 / sqrt(2) , 0, 1 / sqrt(2), 0,
    // 1 / sqrt(2) , 0, -1/ sqrt(2), 0,
    // 0, 1 / sqrt(2) , 0, 1 / sqrt(2),
    // 0, 1 / sqrt(2) , 0, -1/ sqrt(2);
    // CoeffMatrix << 
    // 0.5, 0.5, 0.5, 0.5,
    // 0.5, 0.5, -0.5, -0.5,
    // 0.5, -0.5, 0.5, -0.5,
    // 0.5, -0.5, -0.5, 0.5;

        //Hydrogen Ring

        std::ifstream CoeffFile(BaseName + "_coeff_" + std::to_string(r));
        std::ifstream OverlapFile(BaseName + "_ovlp_" + std::to_string(r));
        int CoeffCount = 0;
        while(CoeffFile >> tmpDouble)
        {
            CoeffMatrix(CoeffCount / NumAO, CoeffCount % NumAO) = tmpDouble;
            CoeffCount++;
        }
        CoeffCount = 0;
        while(OverlapFile >> tmpDouble)
        {
            OverlapMatrix(CoeffCount / NumAO, CoeffCount % NumAO) = tmpDouble;
            CoeffCount++;
        }

        // for(int i = 0; i < CoeffMatrix.cols(); i++)
        // {
        //     double norm = 0;
        //     for(int j = 0; j < CoeffMatrix.rows(); j++)
        //     {
        //         norm += CoeffMatrix(j, i) * CoeffMatrix(j, i);
        //     }
        //     norm = sqrt(norm);
        //     for(int j = 0; j < CoeffMatrix.rows(); j++)
        //     {
        //         CoeffMatrix(j, i) /= norm;
        //     }
        // }

        InvCoeff = CoeffMatrix.transpose();

        // for(int i = 0; i < InvCoeff.cols(); i++)
        // {
        //     double norm = 0;
        //     for(int j = 0; j < InvCoeff.rows(); j++)
        //     {
        //         norm += InvCoeff(j, i) * InvCoeff(j, i);
        //     }
        //     norm = sqrt(norm);
        //     for(int j = 0; j < InvCoeff.rows(); j++)
        //     {
        //         InvCoeff(j, i) /= norm;
        //     }
        // }

        Eigen::MatrixXd SCS = OverlapMatrix * CoeffMatrix * OverlapMatrix.transpose();
        
        Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > SolveS(OverlapMatrix);
        Eigen::MatrixXd SOrtho = Eigen::MatrixXd::Zero(NumAO, NumAO);
        for (int ii = 0; ii < NumAO; ii++)
        {
            SOrtho(ii, ii) = sqrt(SolveS.eigenvalues()[ii]);
        }
        SOrtho = SolveS.eigenvectors() * SOrtho * SolveS.eigenvectors().transpose();
        Eigen::MatrixXd Cp = SOrtho * CoeffMatrix;

        std::cout << Cp << "\n*****************\n" << SCS << "\n*************\n" << Cp * Cp.transpose() << std::endl;



        std::ofstream idumpAO(BaseName + "_loc_" + std::to_string(r));
        idumpAO << std::setprecision(12);

        for (int s = 0; s < OptionLines.size(); s++)
        {
            idumpAO << OptionLines[s] << std::endl;
        }

        for(int i = 0; i < NumAO; i++)
        {
            for (int j = 0; j < NumAO; j++)
            {
                for (int k = 0; k < NumAO; k++)
                {
                    for (int l = 0; l < NumAO; l++)
                    {
                        tmpDouble = 0;
                        for (int n = 0; n < NumAO; n++)
                        {
                            for (int m = 0; m < NumAO; m++)
                            {
                                for (int o = 0; o < NumAO; o++)
                                {
                                    for (int p = 0; p < NumAO; p++)
                                    {
                                        for(int a = 0; a < NumAO; a++)
                                        {
                                            for(int b = 0; b < NumAO; b++)
                                            {
                                                for(int c = 0; c < NumAO; c++)
                                                {
                                                    for(int d = 0; d < NumAO; d++)
                                                    {
                                                        tmpDouble +=  OverlapMatrix(n, i) * OverlapMatrix(m, j) * OverlapMatrix(k, o) * OverlapMatrix(l, p)
                                                                    * InvCoeff(a, n) * InvCoeff(b, m) * InvCoeff(c, o) * InvCoeff(d, p)
                                                                    * Integrals[std::to_string(a + 1) + " " + std::to_string(b + 1) + " " + std::to_string(c + 1) + " " + std::to_string(d + 1)];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        idumpAO << tmpDouble << "\t" << i + 1 << "\t" << j + 1 << "\t" << k + 1 << "\t" << l + 1 << std::endl;
                        // idumpAO << tmpDouble << "\t" << TwoElectronIndexes[i][0] << "\t" << TwoElectronIndexes[i][1] << "\t" << TwoElectronIndexes[i][2] << "\t" << TwoElectronIndexes[i][3] << std::endl;
                    }
                }
            }
        }
        for (int i = 0; i < NumAO; i++)
        {
            for (int j = 0; j < NumAO; j++)
            {
                tmpDouble = 0;
                for (int n = 0; n < NumAO; n++)
                {
                    for (int m = 0; m < NumAO; m++)
                    {
                        for(int a = 0; a < NumAO; a++)
                        {
                            for(int b = 0; b < NumAO; b++)
                            {
                                tmpDouble +=  OverlapMatrix(n, i) * OverlapMatrix(j, m) * InvCoeff(a, n) * InvCoeff(b, m)
                                            * Integrals[std::to_string(a + 1) + " " + std::to_string(b + 1) + " 0 0"];
                            }
                        }
                    }
                }
                idumpAO << tmpDouble << "\t" << i + 1 << "\t" << j + 1 << "\t0\t0" << std::endl;
                // idumpAO << tmpDouble << "\t" << OneElectronIndexes[i][0] << "\t" << OneElectronIndexes[i][1] << "\t0\t0" << std::endl;
            }
        }

        // Now write the orbital energies and nuclear repulsion
        for (int i = 0; i < NumAO + 1; i++)
        {
            idumpAO << Integrals[std::to_string(i) + " 0 0 0"] << "\t" << i << "\t0\t0\t0" << std::endl;
        }

        idumpAO.close();
        // std::ifstream idumpnew("INTDUMP_AO");
        // std::map< std::string, double > IntegralsNew;
        // while(!idumpnew.eof())
        // {
        //     /* We have to include all 8-fold permuation symmetries. This holds each integral in chemistry notation. We represent
        //     (ij|kl) as "i j k l". h_ij is "i j 0 0", as given in QChem. */
        //     idumpnew >> tmpDouble >> tmpInt1 >> tmpInt2 >> tmpInt3 >> tmpInt4;
        //     // IntegralsFile >> tmpInt1 >> tmpInt2 >> tmpInt3 >> tmpInt4 >> tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
        //     IntegralsNew[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
        // }

        // std::string tmpstring;
        // std::cout << IntegralsNew["1 1 1 1"];
        // while(tmpstring != "0")
        // {
        //     std::getline(std::cin, tmpstring);
        //     std::cout << IntegralsNew[tmpstring] << std::endl;
        // }
    }
    return 0;
}

