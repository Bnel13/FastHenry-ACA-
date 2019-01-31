#include <Eigen/Dense>
extern "C"
{
#include "induct.h"
}

//#include <cmath>


/*Starts by applying QR decomposition using the householder transform and then applies 
SVD to obtain the recompressed U and V matrices.
The Eigen Library is used to assist in computing QR and SVD.*/
int SVD_QR(double **U, double **V, int m, int n, int rank, double tol)
{
	using namespace Eigen;

	MatrixXd U2(m, rank);
	MatrixXd V2(n, rank);
	MatrixXd Q_u(m, rank);
	MatrixXd R_u(rank, rank);
	MatrixXd Q_v(n, rank);
	MatrixXd R_v(rank, rank);

	int new_rank = 0;

	for (int i = 0; i < rank; i++)
	{
		for (int j = 0; j < m; j++)
		{
			U2(j, i) = U[i][j];
		}
		for (int j = 0; j < n; j++)
		{
			V2(j, i) = V[i][j];
		}
	}

	HouseholderQR<MatrixXd> qr(U2);
	Q_u = qr.householderQ()* MatrixXd::Identity(m, rank);

	R_u = Q_u.transpose()*U2;

	HouseholderQR<MatrixXd> qr2(V2);
	//Q_v = V2.colPivHouseholderQr().householderQ().setLength(V2.colPivHouseholderQr().nonzeroPivots());
	Q_v = qr2.householderQ()*MatrixXd::Identity(n, rank);

	//R_v = V2.colPivHouseholderQr().matrixR().topLeftCorner(rank, rank).template triangularView<Upper>();
	R_v = Q_v.transpose()*V2;

	JacobiSVD<MatrixXd> svd(R_u*R_v.transpose(), ComputeFullU | ComputeFullV);
	//BDCSVD<MatrixXd> svd(R_u.block(0, 0, rank, rank)*R_v.transpose().block(0, 0, rank, rank), ComputeFullU | ComputeFullV);


	for (int i = 0; i < rank; i++)
	{
		if (svd.singularValues()(i) >= tol*svd.singularValues()(0))
		{
			++new_rank;
		}
		else
		{
			for (int j = 0; j < m; j++)
			{
				U[i][j] = 0;
			}
			for (int j = 0; j < n; j++)
			{
				V[i][j] = 0;
			}
		}
	}

	double temp_U = 0;
	double temp_V = 0;
	MatrixXd U_temp(new_rank, new_rank);


	for (int i = 0; i < new_rank; i++)
	{
		for (int j = 0; j < new_rank; j++)
		{
			U_temp(i, j) = svd.matrixU()(j, i)*svd.singularValues()(i);//check
		}
	}



	for (int i = 0; i < m; i++)
	{
		for (int k = 0; k < new_rank; k++)
		{

			for (int j = 0; j < new_rank; j++)
			{
				temp_U += Q_u(i, j)*U_temp(k, j);
			}
			U[k][i] = temp_U;
			temp_U = 0;
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < new_rank; k++)
		{

			for (int j = 0; j < new_rank; j++)
			{
				temp_V += Q_v(i, j)*svd.matrixV()(j, k);
			}
			V[k][i] = temp_V;
			temp_V = 0;
		}
	}
	return new_rank;
}
